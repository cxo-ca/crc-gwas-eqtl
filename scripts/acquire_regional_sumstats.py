"""Acquire unthresholded regional GWAS/GTEx v8 data; do not run colocalization.

Run from the repository root. The fixed +/-500 kb download windows are provisional
and must be reviewed for LD and cis-window coverage before scientific analysis.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import json
import math
import time
from pathlib import Path
from datetime import datetime, timezone

from lib_remote_tabix import RemoteTabix

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/processed/acquisition'
RAW = ROOT / 'data/raw/summary_statistics'
GWAS_URL = ('https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/'
            'GCST90129001-GCST90130000/GCST90129505/harmonised/GCST90129505.h.tsv.gz')
EQTL_BASE = 'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/'
SOURCES = {'GCST90129505': GWAS_URL,
           **{t: EQTL_BASE + t + '.tsv.gz' for t in ('Colon_Sigmoid', 'Colon_Transverse')}}


def read_tsv(path):
    with open(path, encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def utc():
    return datetime.now(timezone.utc).isoformat()


def finish_file(partial, destination):
    # Windows scanners may briefly hold a newly closed compressed file.
    for attempt in range(10):
        try:
            partial.replace(destination)
            return
        except PermissionError:
            if attempt == 9:
                raise
            time.sleep(0.25 * (attempt + 1))


def make_plan():
    evidence = read_tsv(ROOT / 'results/public/variant_gene_evidence.tsv')
    leads = sorted({(r['chr'].removeprefix('chr'), int(r['pos']), r['rsid']) for r in evidence},
                   key=lambda r: (23 if r[0] == 'X' else int(r[0]), r[1]))
    regions = []
    for chrom, pos, rsid in leads:
        start, end = max(1, pos - 500_000), pos + 500_000
        if regions and regions[-1]['chr'] == chrom and start <= regions[-1]['end']:
            regions[-1]['end'] = max(end, regions[-1]['end'])
            regions[-1]['lead_rsids'] += ';' + rsid
        else:
            regions.append(dict(chr=chrom, start=start, end=end, lead_rsids=rsid))
    for region in regions:
        region['locus'] = f"chr{region['chr']}_{region['start']}_{region['end']}"
        region['build'] = 'GRCh38'
        region['definition'] = 'merged lead +/-500kb download window; provisional, not LD-defined'
    targets = []
    for r in evidence:
        region = next(w for w in regions if r['rsid'] in w['lead_rsids'].split(';'))
        row = dict(locus=region['locus'], gene_id=r['gene_id'],
                   source_gene_id=r['gene_id'].split('.')[0], gene_symbol=r['gene_symbol'],
                   tissue=r['tissue'])
        if row not in targets:
            targets.append(row)
    write_tsv(OUT / 'regions.tsv', regions,
              ['locus', 'chr', 'start', 'end', 'lead_rsids', 'build', 'definition'])
    write_tsv(OUT / 'targets.tsv', targets,
              ['locus', 'gene_id', 'source_gene_id', 'gene_symbol', 'tissue'])
    return regions, targets


def acquire_source(name, regions, targets, offline):
    selected = regions if name == 'GCST90129505' else [
        r for r in regions if any(t['tissue'] == name and t['locus'] == r['locus'] for t in targets)]
    source = RemoteTabix(SOURCES[name], RAW / 'range_cache' / name, offline=offline)
    header = source.get_header().lstrip('#').rstrip('\r\n')
    columns = header.split('\t')
    records = []
    for region in selected:
        locus = region['locus']
        genes = sorted({t['source_gene_id'] for t in targets if t['tissue'] == name and t['locus'] == locus})
        query = dict(url=SOURCES[name], chromosome=region['chr'], start=region['start'],
                     end=region['end'], gene_filter=genes, pvalue_filter=None)
        path = RAW / 'regions' / name / (locus + '.tsv.gz')
        meta_path = path.with_suffix(path.suffix + '.source.json')
        if region['chr'] not in source.references and 'chr' + region['chr'] not in source.references:
            meta = dict(query=query, rows=0, downloaded_utc=utc(),
                        path=path.relative_to(ROOT).as_posix(), status='reference_absent_from_source_index',
                        available_references=list(source.references), source=name)
            meta_path.parent.mkdir(parents=True, exist_ok=True)
            meta_path.write_text(json.dumps(meta, indent=2) + '\n', encoding='utf-8')
            records.append(meta)
            print(f'{name} {locus}: chromosome absent from source index; recorded as unavailable', flush=True)
            continue
        if path.exists() and meta_path.exists():
            meta = json.loads(meta_path.read_text(encoding='utf-8'))
            if meta['query'] != query or sha256(path) != meta['sha256']:
                raise ValueError(f'Existing region cache differs or is damaged: {path}')
            records.append(meta)
            print(f'{name} {locus}: cached {meta["rows"]:,} rows', flush=True)
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        count = 0
        gene_column = columns.index('gene_id') if genes else None
        partial = path.with_suffix(path.suffix + '.part')
        with gzip.open(partial, 'wt', encoding='utf-8', newline='') as handle:
            handle.write(header + '\n')
            for line in source.fetch(region['chr'], region['start'], region['end']):
                line = line.rstrip('\r\n')
                if gene_column is not None and line.split('\t')[gene_column].split('.')[0] not in genes:
                    continue
                handle.write(line + '\n')
                count += 1
        if not count:
            raise ValueError(f'No selected associations for {name} {locus}; investigate, not a success')
        finish_file(partial, path)
        meta = dict(query=query, rows=count, downloaded_utc=utc(),
                    path=path.relative_to(ROOT).as_posix(), sha256=sha256(path),
                    bytes=path.stat().st_size, header=columns)
        meta_path.write_text(json.dumps(meta, indent=2) + '\n', encoding='utf-8')
        records.append(meta)
        print(f'{name} {locus}: saved {count:,} rows', flush=True)
    return records


def number(text):
    try:
        return float(text)
    except (ValueError, TypeError):
        return math.nan


def validate_records(regions, targets):
    records = []
    for name, url in SOURCES.items():
        for region in regions:
            genes = sorted({t['source_gene_id'] for t in targets if t['tissue'] == name and t['locus'] == region['locus']})
            if name != 'GCST90129505' and not genes:
                continue
            path = RAW / 'regions' / name / (region['locus'] + '.tsv.gz')
            meta_path = path.with_suffix(path.suffix + '.source.json')
            if not meta_path.exists():
                raise FileNotFoundError(f'Expected acquisition is incomplete: {meta_path}')
            rec = json.loads(meta_path.read_text(encoding='utf-8'))
            expected = dict(url=url, chromosome=region['chr'], start=region['start'],
                            end=region['end'], gene_filter=genes, pvalue_filter=None)
            if rec['query'] != expected or rec['path'] != path.relative_to(ROOT).as_posix():
                raise ValueError(f'Acquisition query differs from current plan: {meta_path}')
            if rec.get('status') == 'reference_absent_from_source_index':
                source = RemoteTabix(url, RAW / 'range_cache' / name, offline=True)
                if region['chr'] in source.references or 'chr' + region['chr'] in source.references:
                    raise ValueError(f'Unavailable status contradicts cached index: {meta_path}')
            elif not path.exists() or sha256(path) != rec['sha256']:
                raise ValueError(f'Raw region checksum mismatch: {path}')
            records.append(rec)
    return records


def stage(records, targets):
    """Map column names only; retain invalid/duplicate records for subsequent QC."""
    sdy_path = OUT / 'gtex_expression_sdY.tsv'
    if not sdy_path.exists():
        raise FileNotFoundError('Acquire normalized expression with fetch_gtex_expression.py first')
    expr_manifest = json.loads((RAW / 'GTEx_v8/expression/acquisition_manifest.json').read_text(encoding='utf-8'))
    if sha256(sdy_path) != expr_manifest['output_sha256']:
        raise ValueError('Measured sdY file checksum differs from expression manifest')
    sdy_rows = read_tsv(sdy_path)
    sdy = {(r['gene_id'], r['tissue']): r for r in sdy_rows}
    if len(sdy) != len(sdy_rows):
        raise ValueError('Duplicated gene/tissue entries in measured sdY file')
    target_map = {(r['locus'], r['source_gene_id'], r['tissue']): r['gene_id'] for r in targets}
    common = ['locus', 'chr', 'pos', 'effect_allele', 'other_allele', 'beta', 'se', 'p', 'build']
    gcols = common + ['source_rsid', 'hm_code', 'hm_coordinate_conversion']
    ecols = common + ['gene_id', 'tissue', 'sdY', 'source_gene_id', 'source_variant',
                      'source_rsid', 'maf', 'an']
    audit = []
    target_counts = {key: 0 for key in target_map}
    for kind, columns in [('gwas', gcols), ('eqtl', ecols)]:
        dest = OUT / (kind + '.regional.tsv.gz')
        partial = dest.with_suffix(dest.suffix + '.part')
        with gzip.open(partial, 'wt', encoding='utf-8', newline='') as output:
            writer = csv.DictWriter(output, columns, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            for rec in records:
                if rec.get('status') == 'reference_absent_from_source_index':
                    continue
                name = Path(rec['path']).parent.name
                if (kind == 'gwas') != (name == 'GCST90129505'):
                    continue
                locus = Path(rec['path']).name.removesuffix('.tsv.gz')
                counts = dict(source=name, locus=locus, rows=0, nonsignificant_p_gt_0_05=0,
                              missing_or_invalid_beta_se=0, invalid_p=0, missing_sdY=0,
                              identical_alleles=0, duplicate_variant_gene_rows=0)
                seen = set()
                with gzip.open(ROOT / rec['path'], 'rt', encoding='utf-8', newline='') as handle:
                    for r in csv.DictReader(handle, delimiter='\t'):
                        if kind == 'gwas':
                            row = dict(locus=locus, chr=r['chromosome'], pos=r['base_pair_location'],
                                       effect_allele=r['effect_allele'], other_allele=r['other_allele'],
                                       beta=r['beta'], se=r['standard_error'], p=r['p_value'],
                                       build='GRCh38', source_rsid=r.get('rsid', r.get('rs_id', '')),
                                       hm_code=r.get('hm_code', ''),
                                       hm_coordinate_conversion=r.get('hm_coordinate_conversion', ''))
                            ident = (row['chr'], row['pos'], row['effect_allele'], row['other_allele'])
                        else:
                            source_gene = r['gene_id'].split('.')[0]
                            key = (locus, source_gene, name)
                            gene = target_map[key]
                            sd = sdy.get((gene, name), {}).get('sdY', 'NA')
                            row = dict(locus=locus, chr=r['chromosome'], pos=r['position'],
                                       effect_allele=r['alt'], other_allele=r['ref'], beta=r['beta'],
                                       se=r['se'], p=r['pvalue'], build='GRCh38', gene_id=gene,
                                       tissue=name, sdY=sd, source_gene_id=r['gene_id'],
                                       source_variant=r['variant'], source_rsid=r['rsid'],
                                       maf=r['maf'], an=r['an'])
                            ident = (gene, r['variant'])
                            target_counts[key] += 1
                            counts['missing_sdY'] += not (math.isfinite(number(sd)) and number(sd) > 0)
                        counts['rows'] += 1
                        counts['nonsignificant_p_gt_0_05'] += number(row['p']) > 0.05
                        counts['invalid_p'] += not (0 < number(row['p']) <= 1)
                        counts['missing_or_invalid_beta_se'] += not (
                            math.isfinite(number(row['beta'])) and math.isfinite(number(row['se'])) and number(row['se']) > 0)
                        counts['identical_alleles'] += row['effect_allele'] == row['other_allele']
                        counts['duplicate_variant_gene_rows'] += ident in seen
                        seen.add(ident)
                        writer.writerow(row)
                audit.append(counts)
        finish_file(partial, dest)
    write_tsv(OUT / 'acquisition_qc.tsv', audit, list(audit[0]))
    missing_gwas = {Path(r['path']).name.removesuffix('.tsv.gz') for r in records
                    if Path(r['path']).parent.name == 'GCST90129505' and r.get('status') == 'reference_absent_from_source_index'}
    coverage = [dict(**r, rows=target_counts[(r['locus'], r['source_gene_id'], r['tissue'])],
                     gwas_available=r['locus'] not in missing_gwas) for r in targets]
    write_tsv(OUT / 'target_coverage.tsv', coverage, list(coverage[0]))
    if any(r['rows'] == 0 for r in coverage):
        raise ValueError('Some candidate gene/tissue targets have no rows; inspect target_coverage.tsv')
    if any(r['missing_sdY'] for r in audit):
        raise ValueError('Measured phenotype SD is missing for acquired eQTL targets')
    availability = [dict(source=Path(r['path']).parent.name,
                         locus=Path(r['path']).name.removesuffix('.tsv.gz'),
                         status=r.get('status', 'acquired'), rows=r['rows']) for r in records]
    write_tsv(OUT / 'source_availability.tsv', availability, list(availability[0]))
    paired = [r for r in coverage if r['gwas_available']]
    summary = dict(completed_utc=utc(), status='acquired_available_regions_for_qc_with_documented_gaps' if missing_gwas else 'acquired_for_qc_not_coloc_ready',
                   studies=SOURCES, regions=len({r['locus'] for r in targets}),
                   genes=len({r['gene_id'] for r in targets}), gene_tissue_pairs=len({(r['gene_id'], r['tissue']) for r in targets}),
                   target_groups=len(targets), pvalue_filter=None,
                   unavailable_gwas_regions=sorted(missing_gwas),
                   paired_regions=len({r['locus'] for r in paired}),
                   paired_genes=len({r['gene_id'] for r in paired}),
                   paired_gene_tissue_pairs=len({(r['gene_id'], r['tissue']) for r in paired}),
                   note='Provisional merged +/-500kb windows. No allele harmonization, deduplication, QC exclusions or colocalization performed.',
                   rows={kind: sum(r['rows'] for r in audit if (r['source'] == 'GCST90129505') == (kind == 'gwas')) for kind in ('gwas', 'eqtl')},
                   files=[dict(path=p.relative_to(ROOT).as_posix(), bytes=p.stat().st_size, sha256=sha256(p))
                          for p in (OUT/'gwas.regional.tsv.gz', OUT/'eqtl.regional.tsv.gz', sdy_path)])
    (OUT / 'acquisition_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--plan-only', action='store_true')
    parser.add_argument('--workers', type=int, default=3, choices=range(1, 4))
    parser.add_argument('--source', choices=list(SOURCES))
    parser.add_argument('--stage-only', action='store_true')
    args = parser.parse_args()
    regions, targets = make_plan()
    print(f'Plan: {len(regions)} regions, {len(targets)} gene/tissue/locus targets', flush=True)
    if args.plan_only:
        return
    if args.stage_only:
        records = validate_records(regions, targets)
    else:
        names = [args.source] if args.source else list(SOURCES)
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            jobs = [pool.submit(acquire_source, name, regions, targets, args.offline) for name in names]
            records = [r for job in jobs for r in job.result()]
        if args.source:
            return
        records = validate_records(regions, targets)
    stage(records, targets)


if __name__ == '__main__':
    main()
