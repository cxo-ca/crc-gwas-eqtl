"""Reclassify all discovery candidates after omission audit; preserve the coloc run."""
import csv
import html
import itertools
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from report_real_coloc import rows, write_tsv, sha, group, labf, logsum

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/final'
LABELS = {
    'A_conditional_prior_stable': 'A: 세 prior 지지·자료 제한 조건부',
    'B_baseline_prior_sensitive': 'B: 기본 prior 지지·민감',
    'C_permissive_prior_only': 'C: 큰 prior에서만 지지',
    'D_no_threshold_support': 'D: 시험 범위 기준 미충족',
    'E_data_unavailable': 'E: 자료 부재·보류',
}
DIRECTIONS = {
    'concordant': '발현 증가 방향과 CRC odds 증가 방향 일치',
    'discordant': '발현 증가 방향과 CRC odds 감소 방향 일치',
    'zero_effect': '효과값 0으로 방향 미정',
}


def tier(category):
    return {'robust_across_tested_priors': 'A_conditional_prior_stable',
            'baseline_only_sensitive': 'B_baseline_prior_sensitive',
            'permissive_prior_only': 'C_permissive_prior_only',
            'below_threshold': 'D_no_threshold_support'}.get(category, 'E_data_unavailable')


def main():
    audit = json.loads((OUT / 'audit_summary.json').read_text())
    assert audit['pipeline_omissions_found'] == 0
    run = ROOT / audit['input_run']
    for manifest in [OUT / 'audit_input_sha256.tsv', run / 'input_sha256.tsv', run / 'output_sha256.tsv']:
        for item in rows(manifest):
            if sha(ROOT / item['path']) != item['sha256']:
                raise ValueError('Changed audited input: ' + item['path'])
    initial = list(rows(ROOT / 'results/public/candidate_genes.tsv'))
    targets = list(rows(ROOT / 'data/processed/acquisition/targets.tsv'))
    groups = list(rows(run / 'group_summary.tsv'))
    skips = list(rows(run / 'skipped_groups.tsv'))
    posterior_check = json.loads((run / 'posterior_validation.json').read_text())
    assert posterior_check['status'] == 'passed'
    peak_audit = {(group(r), r['missing_trait']): r for r in rows(OUT / 'missing_peak_audit.tsv')}
    warnings = defaultdict(set)
    for r in rows(run / 'warnings.tsv'):
        if r['message'].startswith('minimum p'):
            warnings[group(r)].add('weak_signal_warning')
    representatives = {}
    for k, it in itertools.groupby(rows(ROOT / 'data/processed/harmonization/shared.aligned.tsv.gz'), group):
        block = list(it)
        logs = [labf(float(r['beta_gwas']), float(r['se_gwas']), .2) +
                labf(float(r['beta_eqtl_aligned']), float(r['se_eqtl']), .15 * float(r['sdY'])) for r in block]
        idx = max(range(len(block)), key=lambda i: logs[i])
        r = block[idx]
        product = float(r['beta_gwas']) * float(r['beta_eqtl_aligned'])
        representatives[k] = dict(representative_key=r['key'], representative_rsids=r['gwas_rsids'],
            effect_allele=r['effect_allele'], other_allele=r['other_allele'],
            beta_gwas=r['beta_gwas'], beta_eqtl_aligned=r['beta_eqtl_aligned'],
            representative_snp_pp_given_H4=math.exp(logs[idx] - logsum(logs)),
            representative_direction='concordant' if product > 0 else 'discordant' if product < 0 else 'zero_effect')
    final_groups = []
    for r in groups:
        k = group(r)
        flags = set(warnings[k])
        if r['strongest_clean_gwas_in_shared'] != 'True':
            flags.add('missing_strongest_gwas')
        if r['strongest_clean_eqtl_in_shared'] != 'True':
            flags.add('missing_strongest_eqtl')
        if r['missing_region_leads']:
            flags.add('region_lead_missing')
        if r['prior_adjusted'] == 'True':
            flags.add('package_prior_adjustment')
        if int(r['covered_bins_of_10']) < 10:
            flags.add('partial_window_coverage')
        final_groups.append(dict(r, evidence_tier=tier(r['category']),
            gwas_peak_absence_reason=peak_audit.get((k, 'eqtl'), {}).get('reason', ''),
            eqtl_peak_absence_reason=peak_audit.get((k, 'gwas'), {}).get('reason', ''),
            review_flags=';'.join(sorted(flags)), **representatives[k]))
    cols = list(final_groups[0])
    for r in skips:
        row = {c: '' for c in cols}
        row.update({c: r[c] for c in ('locus', 'gene_id', 'gene_symbol', 'tissue')})
        row.update(evidence_tier='E_data_unavailable', category='not_analyzed',
                   review_flags='gwas_X_region_unavailable', n_shared=0)
        final_groups.append(row)
    final_groups.sort(key=lambda r: (r['evidence_tier'], -float(r['PP_H4_min'] or -1), r['gene_id'], r['tissue']))
    assert {group(r) for r in final_groups} == {group(r) for r in targets}
    assert len(final_groups) == len(targets) == 147
    write_tsv(OUT / 'candidate_groups.tsv', final_groups)
    gene_groups = defaultdict(list)
    for r in final_groups:
        gene_groups[r['gene_id']].append(r)
    candidates = []
    for original in initial:
        gs = gene_groups[original['gene_id']]
        best = min(gs, key=lambda r: (r['evidence_tier'], -float(r['PP_H4_min'] or -1), -float(r['PP_H4_baseline'] or -1)))
        supported = [r for r in gs if r['evidence_tier'][0] in 'AB']
        directions = {r['representative_direction'] for r in supported}
        candidates.append(dict(gene_id=original['gene_id'], gene_symbol=original['gene_symbol'],
            evidence_tier=best['evidence_tier'], classification=LABELS[best['evidence_tier']],
            target_tissues=';'.join(sorted(r['tissue'] for r in gs)),
            analyzed_groups=sum(r['evidence_tier'][0] != 'E' for r in gs),
            prior_stable_tissues=';'.join(sorted(r['tissue'] for r in gs if r['evidence_tier'][0] == 'A')),
            baseline_supported_tissues=';'.join(sorted(r['tissue'] for r in supported)),
            representative_tissue=best['tissue'], representative_locus=best['locus'],
            representative_PP_H4_baseline=best['PP_H4_baseline'], representative_PP_H4_min=best['PP_H4_min'],
            representative_PP_H4_max=best['PP_H4_max'], representative_PP_H3=best['PP_H3_baseline'],
            representative_key=best['representative_key'], representative_rsids=best['representative_rsids'],
            representative_effect_allele=best['effect_allele'], representative_other_allele=best['other_allele'],
            representative_beta_gwas=best['beta_gwas'], representative_beta_eqtl=best['beta_eqtl_aligned'],
            representative_snp_pp_given_H4=best['representative_snp_pp_given_H4'],
            representative_direction=best['representative_direction'],
            supported_groups_direction_discordant=len(directions - {'zero_effect'}) > 1,
            representative_n_shared=best['n_shared'],
            representative_fraction_eqtl_shared=best['fraction_eqtl_shared'],
            representative_peak_absence_reason=best['gwas_peak_absence_reason'],
            review_flags=';'.join(sorted({flag for r in gs for flag in r['review_flags'].split(';') if flag})),
            initial_lead_variants=original['n_lead_variants'],
            interpretation='conditional_exploratory_evidence_not_causal_gene_confirmation'))
    candidates.sort(key=lambda r: (r['evidence_tier'], -float(r['representative_PP_H4_min'] or -1), r['gene_id']))
    assert len(candidates) == len({r['gene_id'] for r in candidates}) == 102
    write_tsv(OUT / 'candidates.tsv', candidates)
    issues = [
        ('processing_integrity', 'checked', '181 input/output checksum records verified; 96 indexed coordinate queries replayed offline', 'No omission found within audited scope; no data rescued'),
        ('missing_GWAS_peaks', 'unresolved_source_coverage', '137 groups: 131 positions absent from imported eQTL source; 6 outside cis window', 'Interpret PP.H4 conditional on shared variants; do not impute missing effects'),
        ('missing_eQTL_peaks', 'unresolved_source_coverage', '4 groups: peak position absent from harmonised GWAS source', 'Retain group-level coverage limitation'),
        ('palindromic_leads', 'intentional_QC_exclusion', '4 lead SNPs; 36 direct discovery gene/tissue pairs affected', 'Keep excluded until strand/frequency evidence supports rescue'),
        ('X_chromosome', 'unresolved_source_absence', '2 genes / 3 groups; source GWAS index contains autosomes only', 'Keep as unanalysed, not negative results'),
        ('GWAS_extreme_effects', 'source_quality_unresolved', 'Some source beta values are extreme; per-cohort estimates/QC not available', 'No arbitrary removal; verify original study quality before causal claims'),
        ('LD_and_multiple_signals', 'not_performed', 'Distance windows; no ancestry-matched LD or conditional/SuSiE analysis', 'Required for stronger multi-signal interpretation'),
        ('reference_and_indels', 'not_performed', 'Provider GRCh38 alignment used; no independent FASTA/indel normalization', 'No claim that all variants or orientations were rescued'),
        ('phenotype_scope', 'documented', 'Combined colorectal cancer study; no separate colon/rectal subtype analysis', 'Do not infer subtype-specific or tumor-specific causality'),
        ('independent_validation', 'not_performed', 'Tissue donors overlap; no external cohort or functional experiments', 'Two tissues are not independent replication'),
        ('discovery_selection', 'documented', '102 candidates selected from significant lead/eQTL overlap', 'Not genome-wide gene discovery; no gene-level FDR or enrichment claims'),
    ]
    write_tsv(OUT / 'issue_ledger.tsv', [dict(item=a, status=b, evidence=c, disposition=d) for a,b,c,d in issues])
    gene_counts = Counter(r['evidence_tier'] for r in candidates)
    group_counts = Counter(r['evidence_tier'] for r in final_groups)
    a_rows = [r for r in candidates if r['evidence_tier'][0] == 'A']
    summary = dict(status='report_built_pending_final_validation',
        completed_utc=datetime.now(timezone.utc).isoformat(), input_run=audit['input_run'],
        candidate_genes=102, candidate_groups=147, analyzed_genes=100, analyzed_groups=144,
        gene_tiers=dict(gene_counts), group_tiers=dict(group_counts),
        prior_stable_genes=[r['gene_symbol'] for r in a_rows],
        audited_pipeline_omissions=0, recovered_rows=0, external_or_functional_validation=False,
        audit_scope=audit['scope'], all_original_candidates_accounted_for=True,
        direction_scope='Association signs at the top conditional-H4 SNP of the representative group; not mediation or gene-wide direction',
        posterior_validation=posterior_check)
    (OUT / 'final_summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')

    def number(value):
        return f'{float(value):.4f}' if value != '' else '—'

    md = ['# 대장암 GWAS·대장 eQTL 후보 재정리 최종 보고서', '',
          '2026-09-23 | 현재 확보 자료에 대한 탐색적 분석의 최종 정리', '',
          '## 먼저 확인한 전 단계 누락', '',
          '체크섬 181건과 원본 색인의 좌표 조회 96건을 대조했습니다. 확인 대상에서 추출·계산 과정의 누락은 발견되지 않았으며 복구된 행은 0개입니다. 모든 미공유 SNP를 전장 유전체 수준에서 점검했다는 뜻은 아닙니다.', '',
          '- 최강 GWAS 변이가 빠진 137개 조합: eQTL 원본에 위치 없음 131개, 유전자 TSS ±1 Mb cis 범위 밖 6개.',
          '- 최강 eQTL 변이가 빠진 4개 조합: 사용한 harmonised GWAS 원본에 위치 없음.',
          '- 모호한 lead 4개(rs2581817, rs209142, rs9271363, rs7398375)는 의도적 QC 제외. 직접 연결된 36개 발견 조합에 영향.',
          '- SHROOM2 두 조직·GPR143 Transverse의 3개 조합은 GWAS X 자료 부재로 보류.',
          '- 점검한 137개 GWAS peak 행의 효과 대립유전자 빈도(EAF)는 모두 NA였습니다. 빈도 대조로 누락 원인이나 모호한 SNP 방향을 해결하지 못했습니다.',
          '- LD·복수 신호·독립 reference FASTA·indel 복구·외부 및 기능 검증은 미수행 항목으로 남겼습니다.', '',
          '## 자료·방법·검증', '',
          'GCST90129505의 GRCh38 대장암 GWAS와 GTEx v8 Colon_Sigmoid/Colon_Transverse를 사용했습니다. 1차 중첩 탐색에서 선정한 102개 후보 중 100개·144개 조합에 실제 coloc을 수행했습니다. 다른 암 아형이나 종양 특이 효과를 따로 분석하지 않았습니다.', '',
          'lead ±500 kb 병합 거리 창, beta/SE² 및 실측 sdY, coloc 5.2.3의 형질별 단일 원인 변이 가정을 사용했습니다. p1=p2=10⁻⁴, 요청 p12=10⁻⁶·10⁻⁵·10⁻⁴로 432회 실행했습니다. 29개 조합의 패키지 prior 자동 조정은 실제 적용값과 함께 원래 결과에 남아 있습니다. p값으로 변이를 선별하지 않았습니다.', '',
          '432개 posterior 모두 별도 Python 계산과 일치했습니다(최대 절대 차이 약 5×10⁻¹⁶). 수치 재현성 검증은 원본 품질이나 생물학적 타당성을 보증하지 않습니다.', '',
          '## 후보 분류', '',
          '원래 후보 102개와 대상 147개 조합을 모두 유지했습니다. 유전자별 분류는 가장 강한 지지를 보이는 조합을 요약하며 조직별 posterior를 합산하지 않습니다. 같은 등급에서는 대표 조합의 최소 PP.H4 순서입니다. 확률이나 별도 종합점수를 새로 추정한 순위가 아닙니다.', '',
          '|분류|유전자|조합|', '|---|---:|---:|']
    md += [f'|{LABELS[k]}|{gene_counts[k]}|{group_counts[k]}|' for k in LABELS]
    md += ['', '기본 prior에서 기준을 충족한 18개 유전자 중 9개는 작은 prior에서 지지가 낮아집니다. A 등급 9개 유전자도 공통 변이 집합과 구간 가정에 조건부인 검토 후보입니다.', '',
           '## 세 prior에서 기준을 유지한 후보', '',
           '|유전자|지지 조직|대표 기본 PP.H4|대표 최소 PP.H4|대표 SNP|연관 방향|', '|---|---|---:|---:|---|---|']
    for r in a_rows:
        md.append(f"|{r['gene_symbol']}|{r['prior_stable_tissues']}|{number(r['representative_PP_H4_baseline'])}|{number(r['representative_PP_H4_min'])}|{r['representative_rsids']}|{DIRECTIONS[r['representative_direction']]}|")
    md += ['', '방향은 대표 조합에서 H4 조건부 SNP 확률이 가장 큰 변이의 같은 효과 대립유전자 기준 beta 부호를 비교한 것입니다. 대표 변이의 SNP.PP.H4는 H4가 참이라는 조건부 값이며, 해당 변이의 인과 확률이나 유전자 전체의 매개 효과가 아닙니다. 조직별 방향과 원래 beta는 전체 표에 보존했습니다.', '',
           '## 해석과 다음 검증', '',
           '자료 누락을 해결했다고 간주해 후보를 확정하지 않았습니다. 특히 137개 조합의 최강 GWAS 변이 부재와 원본의 극단 beta 추정치, 제외된 lead, 거리 구간/복수 신호가 후속 검토 대상입니다. A 등급부터 이 항목을 확인한 뒤 ancestry에 맞는 LD·조건부 분석, 독립 대장 eQTL 및 기능 검증으로 이어가야 합니다. 두 조직은 기증자를 공유하므로 독립 재현이 아닙니다.', '',
           'D 등급은 시험한 모델에서 기준 미충족이며 생물학적 관련성이 없다는 뜻이 아닙니다. E 등급은 미분석이며 음성 결과로 취급하지 않습니다. 유의 eQTL 중첩으로 후보를 선정했으므로 전체 유전자 FDR·enrichment·약물 표적 확정으로 확대 해석하지 않습니다.', '',
           '## 파일과 출처', '',
           '- [전체 102개 후보 TSV](../results/final/candidates.tsv)',
           '- [전체 147개 조합·방향·제한 TSV](../results/final/candidate_groups.tsv)',
           '- [누락 원본 대조 TSV](../results/final/missing_peak_audit.tsv)',
           '- [제외 lead 대조 TSV](../results/final/excluded_lead_audit.tsv)',
           '- [미해결 항목과 처리 방침](../results/final/issue_ledger.tsv)',
           '- [coloc 실행 결과와 그림](coloc_report.html)',
           '- [세부 coloc 방법](../docs/COLOC_RESULTS.md)',
           '- [민감도 그림 PNG](final_candidate_sensitivity.png) · [PDF](final_candidate_sensitivity.pdf)',
           '- [GWAS 원본](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90129001-GCST90130000/GCST90129505/harmonised/)',
           '- [GTEx imported 원본](https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/)',
           '- [GTEx cis 1 Mb 설정](https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/README.md)',
           '- [coloc prior 민감도](https://chr1swallace.github.io/coloc/articles/a04_sensitivity.html)', '']
    report_md = ROOT / 'reports/final_report.md'
    report_md.write_text('\n'.join(md), encoding='utf-8')
    # Render the report prose/tables without external assets or a Markdown dependency.
    import re
    def inline(text):
        escaped = html.escape(text)
        return re.sub(r'\[([^\]]+)\]\(([^)]+)\)', r'<a href="\2">\1</a>', escaped)
    body, in_table, in_list = [], False, False
    for line in md:
        is_table, is_list = line.startswith('|'), line.startswith('- ')
        if in_table and not is_table:
            body.append('</tbody></table>'); in_table = False
        if in_list and not is_list:
            body.append('</ul>'); in_list = False
        if is_table:
            if line.startswith('|---'):
                continue
            if not in_table:
                body.append('<table><tbody>'); in_table = True
            body.append('<tr>' + ''.join('<td>' + inline(c.strip()) + '</td>' for c in line.strip('|').split('|')) + '</tr>')
        elif is_list:
            if not in_list:
                body.append('<ul>'); in_list = True
            body.append('<li>' + inline(line[2:]) + '</li>')
        elif line.startswith('## '):
            body.append('<h2>' + inline(line[3:]) + '</h2>')
        elif line.startswith('# '):
            body.append('<h1>' + inline(line[2:]) + '</h1>')
        elif line:
            body.append('<p>' + inline(line) + '</p>')
    if in_list: body.append('</ul>')
    if in_table: body.append('</tbody></table>')
    body.append('<img src="final_candidate_sensitivity.png" alt="11개 조합의 prior별 PP.H4 범위" style="max-width:100%">')
    body.append('<h2>전체 후보 검색</h2><input id="search" aria-label="후보 검색" placeholder="유전자·조직·분류 검색"><table id="candidates"><thead><tr><th>유전자</th><th>분류</th><th>대표 조직</th><th>기본 PP.H4</th><th>최소 PP.H4</th><th>검토 표시</th></tr></thead><tbody>')
    for r in candidates:
        body.append('<tr>' + ''.join('<td>' + html.escape(str(v)) + '</td>' for v in
            [r['gene_symbol'], r['classification'], r['representative_tissue'],
             number(r['representative_PP_H4_baseline']), number(r['representative_PP_H4_min']), r['review_flags']]) + '</tr>')
    body.append('</tbody></table>')
    style = 'body{font:16px/1.7 system-ui,sans-serif;color:#243342;max-width:1350px;margin:32px auto;padding:0 24px}h1{color:#125867}h2{margin-top:36px}table{border-collapse:collapse;width:100%;font-size:13px}td,th{padding:8px;border-bottom:1px solid #d8e1e5;text-align:left}th{background:#eaf1f4}a{color:#087e8b}input{padding:10px;width:320px}@media print{input{display:none}body{font-size:11px;margin:0}h2{break-after:avoid}tr{break-inside:avoid}}'
    (ROOT / 'reports/final_report.html').write_text('<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>CRC 후보 최종 보고서</title><style>' + style + '</style>' + ''.join(body) +
        "<script>document.getElementById('search').addEventListener('input',function(){let q=this.value.toLowerCase();document.querySelectorAll('#candidates tbody tr').forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(q));});</script></html>", encoding='utf-8')
    print(json.dumps(summary, indent=2, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
