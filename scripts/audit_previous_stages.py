"""Trace prior-stage omissions to indexed source records, without changing results."""
import csv
import gzip
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from acquire_regional_sumstats import SOURCES, RAW
from lib_remote_tabix import RemoteTabix
from report_real_coloc import rows, write_tsv, sha, group

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/final'


def variant(r, raw_kind):
    if raw_kind == 'eqtl':
        return ':'.join([r['chromosome'].removeprefix('chr'), r['position'], *sorted([r['ref'], r['alt']])])
    return ':'.join([r['chromosome'], r['base_pair_location'], *sorted([r['effect_allele'], r['other_allele']])])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    run = ROOT / (ROOT / 'results/coloc/latest_real_run.txt').read_text().strip()
    checked = 0
    for name in ('input_sha256.tsv', 'output_sha256.tsv'):
        for r in rows(run / name):
            if sha(ROOT / r['path']) != r['sha256']:
                raise ValueError('Changed analysis artifact: ' + r['path'])
            checked += 1
    groups = list(rows(run / 'group_summary.tsv'))
    sources = {name: RemoteTabix(url, RAW / 'range_cache' / name, offline=True) for name, url in SOURCES.items()}
    headers = {name: s.get_header().lstrip('#').split('\t') for name, s in sources.items()}
    query_cache = {}

    def fetch(name, key):
        chrom, pos, *_ = key.split(':')
        ident = name, chrom, pos
        if ident not in query_cache:
            query_cache[ident] = [dict(zip(headers[name], line.rstrip('\r\n').split('\t')))
                                  for line in sources[name].fetch(chrom, int(pos), int(pos))]
        return query_cache[ident]

    sd = {(r['gene_id'], r['tissue']): r for r in rows(ROOT / 'data/processed/acquisition/gtex_expression_sdY.tsv')}
    audit = []
    for i, r in enumerate(groups):
        if r['strongest_clean_gwas_in_shared'] == 'True':
            continue
        key, gene, tissue = r['strongest_clean_gwas_key'], r['gene_id'], r['tissue']
        source_rows = fetch(tissue, key)
        exact = [a for a in source_rows if variant(a, 'eqtl') == key]
        target = [a for a in exact if a['gene_id'].split('.')[0] == gene.split('.')[0]]
        if target:
            raise ValueError(f'Potential pipeline omission: {key} {gene} {tissue}; resolve before final report')
        dist = int(key.split(':')[1]) - int(sd[(gene, tissue)]['tss'])
        reason = ('outside_documented_cis_window' if abs(dist) > 1_000_000
                  else 'variant_present_but_target_pair_absent' if exact
                  else 'allele_pair_absent_at_present_position' if source_rows
                  else 'position_absent_from_imported_source')
        raw_g = [a for a in fetch('GCST90129505', key) if variant(a, 'gwas') == key]
        if len(raw_g) != 1:
            raise ValueError('GWAS source peak not uniquely reproduced')
        g = raw_g[0]
        audit.append(dict(locus=r['locus'], gene_id=gene, gene_symbol=r['gene_symbol'], tissue=tissue,
            missing_trait='eqtl', key=key, source_rsid=g.get('rsid', ''), source_p=g['p_value'],
            source_beta=g['beta'], source_se=g['standard_error'],
            source_effect_allele_frequency=g.get('effect_allele_frequency', ''),
            tss_distance=dist, source_rows_at_position=len(source_rows),
            source_rows_matching_alleles=len(exact), source_rows_matching_target=len(target),
            reason=reason, evidence='Offline replay of original indexed bytes; all genes at exact position',
            recovered_rows=0))
        if (i + 1) % 25 == 0:
            print(f'Audited {i + 1}/{len(groups)} groups', flush=True)
    for r in groups:
        if r['strongest_clean_eqtl_in_shared'] == 'True':
            continue
        key = r['strongest_clean_eqtl_key']
        source_rows = fetch('GCST90129505', key)
        exact = [a for a in source_rows if variant(a, 'gwas') == key]
        if exact:
            raise ValueError('eQTL peak exists in GWAS source; investigate QC instead of reporting source absence')
        audit.append(dict(locus=r['locus'], gene_id=r['gene_id'], gene_symbol=r['gene_symbol'], tissue=r['tissue'],
            missing_trait='gwas', key=key, source_rsid='', source_p=r['strongest_clean_eqtl_p'],
            source_beta='', source_se='', source_effect_allele_frequency='', tss_distance='',
            source_rows_at_position=len(source_rows), source_rows_matching_alleles=0,
            source_rows_matching_target=0,
            reason='allele_pair_absent_at_present_position' if source_rows else 'position_absent_from_harmonized_gwas',
            evidence='Offline replay of original indexed GWAS bytes at exact position', recovered_rows=0))
    write_tsv(OUT / 'missing_peak_audit.tsv', audit)
    leads = list(rows(run / 'lead_variant_qc.tsv'))
    evidence = list(rows(ROOT / 'results/public/variant_gene_evidence.tsv'))
    lead_audit = []
    for r in leads:
        if r['gwas_qc'] != 'palindromic_snp':
            continue
        ev = [a for a in evidence if a['rsid'] == r['rsid']]
        keys = {':'.join([a['chr'].removeprefix('chr'), a['pos'], *sorted([a['ref'], a['alt']])]) for a in ev}
        assert len(keys) == 1
        key = keys.pop()
        original = [a for a in fetch('GCST90129505', key) if variant(a, 'gwas') == key]
        assert len(original) == 1
        a = original[0]
        lead_audit.append(dict(**r, key=key, effect_allele=a['effect_allele'], other_allele=a['other_allele'],
            source_p=a['p_value'], hm_code=a['hm_code'],
            affected_target_pairs=len({(x['gene_id'], x['tissue']) for x in ev}),
            disposition='intentional_QC_exclusion; not restored without strand/frequency evidence'))
    write_tsv(OUT / 'excluded_lead_audit.tsv', lead_audit)
    missing_x = 'X' not in sources['GCST90129505'].references
    assert missing_x and '23' not in sources['GCST90129505'].references
    qc = list(rows(ROOT / 'data/processed/harmonization/group_qc.tsv'))
    assert len(qc) == 147 and len(groups) == 144
    summary = dict(status='audit_complete_with_unresolved_source_and_method_limits',
        completed_utc=datetime.now(timezone.utc).isoformat(), verified_artifact_hashes=checked,
        missing_gwas_peak_groups=sum(r['missing_trait'] == 'eqtl' for r in audit),
        missing_eqtl_peak_groups=sum(r['missing_trait'] == 'gwas' for r in audit),
        reason_counts=dict(Counter(r['missing_trait'] + ':' + r['reason'] for r in audit)),
        pipeline_omissions_found=0, recovered_rows=0, excluded_leads=len(lead_audit),
        excluded_lead_target_pairs=sum(r['affected_target_pairs'] for r in lead_audit),
        x_reference_absent_from_gwas=True, x_skipped_groups=3,
        scope='Audit of reported missing peaks, excluded leads and X availability; not every absent SNP genome-wide',
        unresolved=['Original GWAS allele frequencies absent where NA; no cohort-level genotype/QC data',
                    'No LD/recombination boundary validation or conditional/SuSiE multi-signal analysis',
                    'No independent reference-FASTA check or palindromic/indel rescue',
                    'Candidate set selected via significant eQTL overlap; not genome-wide gene screening',
                    'No independent cohort or functional validation'],
        input_run=run.relative_to(ROOT).as_posix(), indexed_source_queries=len(query_cache))
    (OUT / 'audit_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    write_tsv(OUT / 'audit_input_sha256.tsv',
              [dict(path=p.relative_to(ROOT).as_posix(), sha256=sha(p)) for p in
               [run / 'group_summary.tsv', run / 'input_sha256.tsv', run / 'output_sha256.tsv',
                Path(__file__), ROOT / 'scripts/lib_remote_tabix.py',
                *[s.manifest_path for s in sources.values()]]])
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
