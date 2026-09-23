"""Offline, auditable SNP QC and GWAS-oriented pairing; never runs coloc."""
import csv
import gzip
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict
from contextlib import ExitStack
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ACQ = ROOT / 'data/processed/acquisition'
OUT = ROOT / 'data/processed/harmonization'
COMMON = ['locus', 'chr', 'pos', 'effect_allele', 'other_allele', 'beta', 'se', 'p', 'build']
PAIR_COLS = ['locus', 'gene_id', 'tissue', 'key', 'chr', 'pos', 'build',
             'effect_allele', 'other_allele', 'beta_gwas', 'se_gwas', 'p_gwas',
             'effect_allele_eqtl_original', 'other_allele_eqtl_original',
             'beta_eqtl_original', 'beta_eqtl_aligned', 'se_eqtl', 'p_eqtl', 'sdY',
             'allele_flipped', 'gwas_source_row', 'eqtl_source_row',
             'gwas_rsids', 'eqtl_rsids', 'eqtl_source_variant']


def sha256(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def read_rows(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8', newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def key(row):
    return ':'.join([row['chr'], row['pos'], *sorted([row['effect_allele'], row['other_allele']])])


def identity(row, kind):
    return (row['locus'], row.get('gene_id', ''), row.get('tissue', ''), key(row))


def normalize(row, kind, regions, sdy):
    row = dict(row)
    row['chr'] = row['chr'].removeprefix('chr').upper()
    row['effect_allele'] = row['effect_allele'].upper()
    row['other_allele'] = row['other_allele'].upper()
    region = regions[row['locus']]
    if row['build'] != 'GRCh38' or region['build'] != 'GRCh38':
        raise ValueError('Genome build mismatch; no automatic liftover')
    pos = int(row['pos'])
    if row['chr'] != region['chr'] or not int(region['start']) <= pos <= int(region['end']):
        raise ValueError('Variant outside its declared region')
    row['pos'] = str(pos)
    if kind == 'eqtl':
        expected = f"chr{row['chr']}_{pos}_{row['other_allele']}_{row['effect_allele']}"
        if row['source_variant'] != expected:
            raise ValueError('GTEx variant ID contradicts GRCh38 ref/alt columns')
        if row['source_gene_id'].split('.')[0] != row['gene_id'].split('.')[0]:
            raise ValueError('Gene ID mapping mismatch')
        if float(row['sdY']) != sdy[(row['gene_id'], row['tissue'])]:
            raise ValueError('Expression SD differs from measured phenotype SD')
    return row


def reasons(row, kind):
    result = []
    for name in ['beta', 'se', 'p'] + (['sdY'] if kind == 'eqtl' else []):
        try:
            value = float(row[name])
            good = math.isfinite(value)
            if name in ('se', 'sdY'):
                good &= value > 0
            if name == 'se':
                good &= math.isfinite(value * value) and value * value > 0
            if name == 'p':
                good &= 0 < value <= 1
        except (ValueError, OverflowError):
            good = False
        if not good:
            result.append('invalid_' + name)
    ea, oa = row['effect_allele'], row['other_allele']
    if ea == oa:
        result.append('identical_alleles')
    if ea not in 'ACGT' or oa not in 'ACGT' or len(ea) != 1 or len(oa) != 1:
        result.append('non_snv')
    elif ''.join(sorted([ea, oa])) in ('AT', 'CG'):
        result.append('palindromic_snp')
    if kind == 'gwas' and row['hm_code'] not in {str(i) for i in range(1, 14) if i != 9}:
        result.append('unharmonized_gwas')
    return result


def deduplicate(rows, kind):
    """Collapse only exact duplicates except rsID. Quarantine ALL conflicting rows."""
    groups = defaultdict(list)
    for row in rows:
        groups[identity(row, kind)].append(row)
    kept, removed = [], []
    for values in groups.values():
        first = values[0]
        ignored = {'source_rsid', 'source_row'}
        signature = {k: v for k, v in first.items() if k not in ignored}
        if any({k: v for k, v in r.items() if k not in ignored} != signature for r in values[1:]):
            removed.extend((r, ['conflicting_duplicate'], '') for r in values)
            continue
        merged = dict(first)
        merged['source_rsid'] = ';'.join(sorted({r['source_rsid'] for r in values}))
        kept.append(merged)
        removed.extend((r, ['duplicate_rsid_alias'], first['source_row']) for r in values[1:])
    return kept, removed


def align(g, e):
    if (g['locus'], g['build'], key(g)) != (e['locus'], e['build'], key(e)):
        raise ValueError('Cannot align different variants/builds')
    same = (g['effect_allele'], g['other_allele']) == (e['effect_allele'], e['other_allele'])
    flip = (g['effect_allele'], g['other_allele']) == (e['other_allele'], e['effect_allele'])
    if not (same or flip):
        raise ValueError('Unresolved allele mismatch')
    return dict(locus=e['locus'], gene_id=e['gene_id'], tissue=e['tissue'], key=key(e),
                chr=e['chr'], pos=e['pos'], build=e['build'],
                effect_allele=g['effect_allele'], other_allele=g['other_allele'],
                beta_gwas=g['beta'], se_gwas=g['se'], p_gwas=g['p'],
                effect_allele_eqtl_original=e['effect_allele'], other_allele_eqtl_original=e['other_allele'],
                beta_eqtl_original=e['beta'], beta_eqtl_aligned=format(-float(e['beta']), '.17g') if flip else e['beta'],
                se_eqtl=e['se'], p_eqtl=e['p'], sdY=e['sdY'], allele_flipped=str(flip).upper(),
                gwas_source_row=g['source_row'], eqtl_source_row=e['source_row'],
                gwas_rsids=g['source_rsid'], eqtl_rsids=e['source_rsid'], eqtl_source_variant=e['source_variant'])


def writer(stack, path, columns):
    opener = gzip.open if path.suffix == '.gz' else open
    handle = stack.enter_context(opener(path, 'wt', encoding='utf-8', newline=''))
    result = csv.DictWriter(handle, columns, delimiter='\t', lineterminator='\n', extrasaction='ignore')
    result.writeheader()
    return result


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # A failed rerun must not leave a previous completion marker looking current.
    summary_path = OUT / 'harmonization_summary.json'
    if summary_path.exists():
        summary_path.replace(OUT / 'harmonization_summary.previous.json')
    acquisition = json.loads((ACQ / 'acquisition_summary.json').read_text(encoding='utf-8'))
    inputs = [ROOT / item['path'] for item in acquisition['files']]
    for item, path in zip(acquisition['files'], inputs):
        if sha256(path) != item['sha256']:
            raise ValueError(f'Acquisition checksum mismatch: {path}')
    regions = {r['locus']: r for r in read_rows(ACQ / 'regions.tsv')}
    targets = list(read_rows(ACQ / 'targets.tsv'))
    expected_targets = {(r['locus'], r['gene_id'], r['tissue']) for r in targets}
    sdy = {(r['gene_id'], r['tissue']): float(r['sdY']) for r in read_rows(ACQ / 'gtex_expression_sdY.tsv')}
    counts, reason_counts, raw_counts = Counter(), Counter(), Counter()
    gwas, gwas_positions, lead_status = {}, {}, defaultdict(set)
    coverage = []
    observed_targets = set()
    with ExitStack() as stack:
        audit = writer(stack, OUT / 'excluded_rows.tsv.gz',
                       ['dataset', 'source_row', 'reasons', 'kept_source_row'] + COMMON +
                       ['gene_id', 'tissue', 'source_rsid', 'source_variant'])
        shared = writer(stack, OUT / 'shared.aligned.tsv.gz', PAIR_COLS)
        unmatched = writer(stack, OUT / 'unmatched_eqtl.tsv.gz',
                           ['locus', 'gene_id', 'tissue', 'key', 'source_row', 'reason'])

        def reject(kind, row, why, kept=''):
            audit.writerow(dict(row, dataset=kind, reasons=';'.join(why), kept_source_row=kept))
            counts[kind + '_excluded_rows'] += 1
            reason_counts.update(kind + ':' + reason for reason in why)

        for kind in ('gwas', 'eqtl'):
            path = ACQ / (kind + '.regional.tsv.gz')
            columns = list(next(read_rows(path))) + ['key', 'source_row']
            cleaned = writer(stack, OUT / (kind + '.cleaned.tsv.gz'), columns)

            def indexed():
                for number, original in enumerate(read_rows(path), 1):
                    row = normalize(original, kind, regions, sdy)
                    if kind == 'eqtl' and (row['locus'], row['gene_id'], row['tissue']) not in expected_targets:
                        raise ValueError('Unexpected eQTL target')
                    row['source_row'] = number  # 1-based data row; excludes header
                    counts[kind + '_raw_rows'] += 1
                    raw_counts[(kind, row['locus'], row.get('gene_id', ''), row.get('tissue', ''))] += 1
                    yield row

            seen_blocks = set()
            for block, iterator in itertools.groupby(indexed(), lambda r: (r['locus'], r.get('tissue', ''))):
                if block in seen_blocks:
                    raise ValueError('Input blocks not contiguous; would miss duplicates')
                seen_blocks.add(block)
                locus, tissue = block
                rows = list(iterator)
                for row in rows:
                    if kind == 'gwas':
                        for rsid in row['source_rsid'].split(';'):
                            lead_status[(locus, rsid)].update(reasons(row, kind) or ['retained'])
                unique, removed = deduplicate(rows, kind)
                for row, why, kept in removed:
                    reject(kind, row, why, kept)
                good = []
                for row in unique:
                    why = reasons(row, kind)
                    if why:
                        reject(kind, row, why)
                    else:
                        row['key'] = key(row)
                        good.append(row)
                        cleaned.writerow(row)
                        counts[kind + '_clean_rows'] += 1
                if kind == 'gwas':
                    gwas[locus] = {r['key']: r for r in good}
                    gwas_positions[locus] = {r['pos'] for r in good}
                    continue
                by_gene = defaultdict(list)
                for row in good:
                    by_gene[row['gene_id']].append(row)
                for target in [t for t in targets if (t['locus'], t['tissue']) == block]:
                    gene = target['gene_id']
                    observed_targets.add((locus, gene, tissue))
                    a, b = gwas.get(locus, {}), by_gene[gene]
                    pairs, misses = [], Counter()
                    for row in b:
                        if row['key'] in a:
                            pair = align(a[row['key']], row)
                            pairs.append(pair)
                            shared.writerow(pair)
                        else:
                            reason = ('gwas_region_unavailable' if locus in acquisition['unavailable_gwas_regions']
                                      else 'allele_mismatch_at_position' if row['pos'] in gwas_positions.get(locus, set())
                                      else 'absent_from_clean_gwas')
                            misses[reason] += 1
                            unmatched.writerow(dict(row, reason=reason))
                    n = len(pairs)
                    counts['shared_rows'] += n
                    flips = sum(p['allele_flipped'] == 'TRUE' for p in pairs)
                    counts['flipped_rows'] += flips
                    region = regions[locus]
                    start, end = int(region['start']), int(region['end'])
                    bins = {min(9, 10 * (int(p['pos']) - start) // (end - start + 1)) for p in pairs}
                    shared_rsids = {rs for p in pairs for rs in p['gwas_rsids'].split(';')} | {rs for p in pairs for rs in p['eqtl_rsids'].split(';')}
                    # Plan joins multiple lead IDs with semicolons.
                    leads = set(region['lead_rsids'].split(';'))
                    status = ('gwas_region_unavailable' if locus in acquisition['unavailable_gwas_regions']
                              else 'meets_min_shared_snps' if n >= 50 else 'insufficient_shared_snps')
                    coverage.append(dict(target, n_gwas_raw=raw_counts[('gwas', locus, '', '')],
                        n_eqtl_raw=raw_counts[('eqtl', locus, gene, tissue)], n_gwas_clean=len(a), n_eqtl_clean=len(b),
                        n_shared=n, n_flipped=flips, fraction_gwas_shared=n / len(a) if a else 0,
                        fraction_eqtl_shared=n / len(b) if b else 0,
                        shared_start=min((int(p['pos']) for p in pairs), default=''),
                        shared_end=max((int(p['pos']) for p in pairs), default=''),
                        covered_bins_of_10=len(bins), n_shared_region_leads=len(leads & shared_rsids),
                        missing_region_leads=';'.join(sorted(leads - shared_rsids)),
                        unmatched_position_or_missing=misses['absent_from_clean_gwas'],
                        allele_mismatch_at_position=misses['allele_mismatch_at_position'], status=status))
                print(f'{kind}: {locus} {tissue} QC complete', flush=True)
        if observed_targets != expected_targets:
            raise ValueError('Missing or unexpected gene/tissue/locus groups')
        for kind in ('gwas', 'eqtl'):
            if counts[kind + '_raw_rows'] != counts[kind + '_clean_rows'] + counts[kind + '_excluded_rows']:
                raise ValueError('Row accounting failed')
        cov = writer(stack, OUT / 'group_qc.tsv', list(coverage[0]))
        cov.writerows(coverage)
        lead = writer(stack, OUT / 'lead_variant_qc.tsv', ['locus', 'rsid', 'gwas_qc'])
        for locus, region in regions.items():
            for rsid in region['lead_rsids'].split(';'):
                lead.writerow(dict(locus=locus, rsid=rsid,
                                   gwas_qc=';'.join(sorted(lead_status.get((locus, rsid), {'not_available'})))))
    eligible = [r for r in coverage if r['status'] == 'meets_min_shared_snps']
    paths = sorted(OUT.glob('*.tsv*'))
    inputs += [ACQ / 'regions.tsv', ACQ / 'targets.tsv', ACQ / 'acquisition_summary.json', Path(__file__)]
    summary = dict(completed_utc=datetime.now(timezone.utc).isoformat(), status='harmonized_not_coloc_run',
                   build='GRCh38', effect_orientation='GWAS effect allele in shared.aligned.tsv.gz',
                   cleaned_eqtl_orientation='Native GTEx ALT; unchanged, including unmatched rows',
                   matching_key='chromosome:1-based_position:alphabetically_sorted_allele_pair',
                   counts=dict(counts), exclusion_reason_counts=dict(reason_counts),
                   groups=len(coverage), groups_meeting_min_shared_snps=len(eligible),
                   genes_meeting_min_shared_snps=len({r['gene_id'] for r in eligible}), min_shared_snps=50,
                   n_groups_with_allele_mismatch=sum(r['allele_mismatch_at_position'] > 0 for r in coverage),
                   limitations=['Provisional merged +/-500kb windows; LD boundaries not established',
                                '50 shared SNPs is an operational minimum, not proof of adequate coverage',
                                'Palindromic SNPs and non-SNVs excluded; no strand complementation or indel normalization',
                                'Source GRCh38 harmonization used; no independent reference FASTA check',
                                'No ancestry-specific allele-frequency/LD comparison; mixed-ancestry GWAS',
                                'X GWAS missing for SHROOM2/GPR143; no colocalization performed'],
                   inputs=[dict(path=p.relative_to(ROOT).as_posix(), sha256=sha256(p)) for p in inputs],
                   files=[dict(path=p.relative_to(ROOT).as_posix(), bytes=p.stat().st_size, sha256=sha256(p)) for p in paths])
    summary_path.write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: summary[k] for k in ('status', 'counts', 'exclusion_reason_counts', 'groups_meeting_min_shared_snps', 'genes_meeting_min_shared_snps')}, indent=2))


if __name__ == '__main__':
    main()
