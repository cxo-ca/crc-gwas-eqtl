"""Validate final accounting, source/audit provenance, report links and exports."""
import json
from collections import Counter
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path

from report_real_coloc import rows, sha, write_tsv, group

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/final'


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        for key, value in attrs:
            if (tag == 'a' and key == 'href') or (tag == 'img' and key == 'src'):
                self.links.append(value)


def main():
    summary = json.loads((OUT / 'final_summary.json').read_text(encoding='utf-8'))
    run = ROOT / summary['input_run']
    checked = 0
    for path in (run / 'input_sha256.tsv', run / 'output_sha256.tsv', OUT / 'audit_input_sha256.tsv'):
        for r in rows(path):
            assert sha(ROOT / r['path']) == r['sha256'], r['path']
            checked += 1
    genes = list(rows(OUT / 'candidates.tsv'))
    groups = list(rows(OUT / 'candidate_groups.tsv'))
    original = list(rows(ROOT / 'results/public/candidate_genes.tsv'))
    targets = list(rows(ROOT / 'data/processed/acquisition/targets.tsv'))
    fits = list(rows(run / 'coloc_sensitivity.tsv'))
    assert len(genes) == len({r['gene_id'] for r in genes}) == len(original) == 102
    assert {r['gene_id'] for r in genes} == {r['gene_id'] for r in original}
    assert len(groups) == len({group(r) for r in groups}) == len(targets) == 147
    assert {group(r) for r in groups} == {group(r) for r in targets}
    assert dict(Counter(r['evidence_tier'] for r in genes)) == summary['gene_tiers']
    assert dict(Counter(r['evidence_tier'] for r in groups)) == summary['group_tiers']
    fit_map = {(group(r), float(r['p12_requested'])): r for r in fits}
    for r in groups:
        if r['evidence_tier'].startswith('E'):
            assert r['PP_H4_baseline'] == '' and r['representative_key'] == ''
            continue
        p = [float(fit_map[(group(r), prior)]['PP_H4']) for prior in (1e-6, 1e-5, 1e-4)]
        assert float(r['PP_H4_min']) == min(p)
        assert float(r['PP_H4_baseline']) == p[1]
        assert float(r['PP_H4_max']) == max(p)
        assert 0 < float(r['representative_snp_pp_given_H4']) <= 1
        product = float(r['beta_gwas']) * float(r['beta_eqtl_aligned'])
        assert r['representative_direction'] == ('concordant' if product > 0 else 'discordant' if product < 0 else 'zero_effect')
    report = ROOT / 'reports/final_report.html'
    parser = Links()
    parser.feed(report.read_text(encoding='utf-8'))
    local_links = [v for v in parser.links if not v.startswith(('https://', 'http://', '#'))]
    for link in local_links:
        assert (report.parent / link).is_file(), link
    assert (ROOT / 'reports/final_candidate_sensitivity.png').read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
    assert (ROOT / 'reports/final_candidate_sensitivity.pdf').read_bytes().startswith(b'%PDF-')
    validation = dict(status='passed', completed_utc=datetime.now(timezone.utc).isoformat(),
        source_checksum_records_checked=checked, candidate_genes=102, candidate_groups=147,
        existing_coloc_results_unchanged=True, local_report_links_checked=len(local_links),
        posterior_values_preserved=True, unanalysed_groups_have_no_posterior=True)
    (OUT / 'final_validation.json').write_text(json.dumps(validation, indent=2) + '\n')
    summary['status'] = 'final_exploratory_report_complete_with_explicit_limits'
    summary['final_validation'] = validation
    (OUT / 'final_summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    paths = [p for p in OUT.iterdir() if p.is_file() and p.name != 'final_manifest.tsv']
    paths += list((ROOT / 'reports').glob('final_*'))
    paths += [ROOT / 'scripts' / name for name in ('audit_previous_stages.py', 'build_final_report.py',
              'plot_final_candidates.R', 'verify_final_report.py')]
    write_tsv(OUT / 'final_manifest.tsv', [dict(path=p.relative_to(ROOT).as_posix(), bytes=p.stat().st_size,
                                             sha256=sha(p)) for p in sorted(paths)])
    print(json.dumps(validation, indent=2))


if __name__ == '__main__':
    main()
