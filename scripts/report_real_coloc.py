"""Validate all real coloc posteriors independently, then write a reviewable report."""
import csv
import gzip
import hashlib
import html
import itertools
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def rows(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8', newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def write_tsv(path, data):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, list(data[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(data)


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def group(r):
    return r['locus'], r['gene_id'], r['tissue']


def logsum(values):
    peak = max(values)
    return peak + math.log(math.fsum(math.exp(x - peak) for x in values))


def labf(beta, se, prior_sd):
    variance = se * se
    shrink = prior_sd * prior_sd / (prior_sd * prior_sd + variance)
    return 0.5 * (math.log1p(-shrink) + shrink * (beta / se) ** 2)


def independent_posteriors(stats, p1, p2, p12):
    s1, s2, shared = stats
    total = s1 + s2
    # The off-diagonal configurations in H3 exclude the shared causal SNP.
    different = total + math.log(-math.expm1(shared - total)) if shared < total else -math.inf
    logs = [0, math.log(p1) + s1, math.log(p2) + s2,
            math.log(p1) + math.log(p2) + different, math.log(p12) + shared]
    norm = logsum(logs)
    return [math.exp(x - norm) for x in logs]


def main():
    run = ROOT / (ROOT / 'results/coloc/latest_real_run.txt').read_text().strip()
    if not run.resolve().is_relative_to((ROOT / 'results/coloc').resolve()):
        raise ValueError('Invalid run path')
    run_info = json.loads((run / 'run_summary.json').read_text())
    assert run_info['status'] == 'coloc_completed' and run_info['real_data']
    for item in rows(run / 'input_sha256.tsv'):
        assert sha(ROOT / item['path']) == item['sha256'], item['path']
    fits = list(rows(run / 'coloc_sensitivity.tsv'))
    by_group = defaultdict(list)
    for r in fits:
        by_group[group(r)].append(r)
    stats, shared_counts, shared_keys = {}, {}, {}
    source = ROOT / 'data/processed/harmonization/shared.aligned.tsv.gz'
    for key, block in itertools.groupby(rows(source), group):
        if key in stats:
            raise ValueError('Non-contiguous shared group')
        abf1, abf2, keys = [], [], set()
        for r in block:
            keys.add(r['key'])
            abf1.append(labf(float(r['beta_gwas']), float(r['se_gwas']), 0.2))
            abf2.append(labf(float(r['beta_eqtl_aligned']), float(r['se_eqtl']), 0.15 * float(r['sdY'])))
        stats[key] = logsum(abf1), logsum(abf2), logsum([a + b for a, b in zip(abf1, abf2)])
        shared_counts[key] = len(abf1)
        shared_keys[key] = keys
    assert set(by_group) == set(stats)
    gwas_peaks, eqtl_peaks = {}, {}
    for kind, peaks in [('gwas', gwas_peaks), ('eqtl', eqtl_peaks)]:
        for r in rows(ROOT / f'data/processed/harmonization/{kind}.cleaned.tsv.gz'):
            k = r['locus'] if kind == 'gwas' else group(r)
            if k not in peaks or float(r['p']) < float(peaks[k]['p']):
                peaks[k] = r
    max_error = 0
    adjusted = set()
    for r in fits:
        key = group(r)
        n = int(r['nsnps'])
        assert n == shared_counts[key]
        requested = [float(r[k + '_requested']) for k in ('p1', 'p2', 'p12')]
        actual = [float(r[k]) for k in ('p1', 'p2', 'p12')]
        for req, act in zip(requested, actual):
            expected = 1 / (n + 1) if n * req >= 1 else req
            assert math.isclose(act, expected, rel_tol=1e-12)
        if any(not math.isclose(a, b, rel_tol=1e-12) for a, b in zip(requested, actual)):
            adjusted.add(key)
        expected = independent_posteriors(stats[key], *actual)
        observed = [float(r[f'PP_H{i}']) for i in range(5)]
        assert all(math.isfinite(v) and 0 <= v <= 1 for v in observed)
        assert abs(sum(observed) - 1) < 1e-10
        err = max(abs(a - b) for a, b in zip(expected, observed))
        max_error = max(max_error, err)
        assert err < 1e-7, (key, err)
    qc = list(rows(run / 'qc.tsv'))
    pre = {group(r): r for r in rows(run / 'preanalysis_qc.tsv')}
    assert len(fits) == run_info['sensitivity_rows'] == 432
    assert len(by_group) == run_info['groups_analyzed'] == 144
    assert len(qc) == len(pre) == 147
    summary_rows, skipped = [], []
    labels = {
        'robust_across_tested_priors': '세 사전확률 모두 PP.H4 ≥ 0.8',
        'baseline_only_sensitive': '기본 기준 충족·사전확률에 민감',
        'permissive_prior_only': '큰 p12에서만 기준 충족',
        'below_threshold': '시험한 범위에서 기준 미충족',
    }
    for i, q in enumerate(qc, 1):
        key = group(q)
        p = pre[key]
        assert int(q['n_shared']) == int(p['n_shared'])
        if q['status'] != 'analyzed':
            assert key not in by_group and int(q['n_shared']) == 0
            skipped.append(dict(locus=key[0], gene_id=key[1], gene_symbol=p['gene_symbol'],
                                tissue=key[2], reason=p['status']))
            continue
        fs = by_group[key]
        assert sorted(float(r['p12_requested']) for r in fs) == [1e-6, 1e-5, 1e-4]
        baseline = next(r for r in fs if float(r['p12_requested']) == 1e-5)
        low, high = min(float(r['PP_H4']) for r in fs), max(float(r['PP_H4']) for r in fs)
        base = float(baseline['PP_H4'])
        category = ('robust_across_tested_priors' if low >= .8 else 'baseline_only_sensitive'
                    if base >= .8 else 'permissive_prior_only' if high >= .8 else 'below_threshold')
        plot = f'locus_{i:03d}.png'
        assert (run / plot).exists()
        summary_rows.append(dict(locus=key[0], gene_id=key[1], gene_symbol=p['gene_symbol'], tissue=key[2],
            n_shared=int(q['n_shared']), PP_H4_baseline=base, PP_H3_baseline=float(baseline['PP_H3']),
            PP_H0_baseline=float(baseline['PP_H0']), PP_H1_baseline=float(baseline['PP_H1']),
            PP_H2_baseline=float(baseline['PP_H2']), PP_H4_min=low, PP_H4_max=high, category=category,
            fraction_gwas_shared=p['fraction_gwas_shared'], fraction_eqtl_shared=p['fraction_eqtl_shared'],
            covered_bins_of_10=p['covered_bins_of_10'], missing_region_leads=p['missing_region_leads'],
            strongest_clean_gwas_key=gwas_peaks[key[0]]['key'], strongest_clean_gwas_p=gwas_peaks[key[0]]['p'],
            strongest_clean_gwas_in_shared=gwas_peaks[key[0]]['key'] in shared_keys[key],
            strongest_clean_eqtl_key=eqtl_peaks[key]['key'], strongest_clean_eqtl_p=eqtl_peaks[key]['p'],
            strongest_clean_eqtl_in_shared=eqtl_peaks[key]['key'] in shared_keys[key],
            prior_adjusted=key in adjusted, plot=plot))
    summary_rows.sort(key=lambda r: (-r['PP_H4_min'], -r['PP_H4_baseline'], r['gene_id'], r['tissue']))
    write_tsv(run / 'group_summary.tsv', summary_rows)
    write_tsv(run / 'skipped_groups.tsv', skipped)
    genes = defaultdict(list)
    for r in summary_rows:
        genes[r['gene_id']].append(r)
    gene_rows = []
    for gene, gs in genes.items():
        gene_rows.append(dict(gene_id=gene, gene_symbol=gs[0]['gene_symbol'],
            analyzed_tissues=';'.join(sorted(r['tissue'] for r in gs)),
            robust_tissues=';'.join(sorted(r['tissue'] for r in gs if r['PP_H4_min'] >= .8)),
            baseline_supported_tissues=';'.join(sorted(r['tissue'] for r in gs if r['PP_H4_baseline'] >= .8)),
            max_PP_H4_baseline=max(r['PP_H4_baseline'] for r in gs),
            best_min_PP_H4=max(r['PP_H4_min'] for r in gs)))
    gene_rows.sort(key=lambda r: -r['best_min_PP_H4'])
    write_tsv(run / 'gene_summary.tsv', gene_rows)
    categories = Counter(r['category'] for r in summary_rows)
    validation = dict(status='passed', method='Independent Python Wakefield ABF and H0-H4 enumeration from aligned inputs',
                      fits_checked=len(fits), shared_rows_checked=sum(shared_counts.values()),
                      max_absolute_posterior_error=max_error, tolerance=1e-7)
    (run / 'posterior_validation.json').write_text(json.dumps(validation, indent=2) + '\n')
    result = dict(status='completed_and_validated', real_data=True, run=run.relative_to(ROOT).as_posix(),
        completed_utc=datetime.now(timezone.utc).isoformat(), coloc_version=run_info['coloc_version'],
        analyzed_genes=len(genes), analyzed_groups=len(summary_rows), skipped_groups=len(skipped),
        fits=len(fits), group_categories=dict(categories),
        baseline_supported_genes=sum(bool(r['baseline_supported_tissues']) for r in gene_rows),
        robust_genes=sum(bool(r['robust_tissues']) for r in gene_rows),
        robust_both_tissues_genes=sum(';' in r['robust_tissues'] for r in gene_rows),
        groups_with_adjusted_priors=len(adjusted),
        groups_missing_strongest_clean_gwas=sum(not r['strongest_clean_gwas_in_shared'] for r in summary_rows),
        groups_missing_strongest_clean_eqtl=sum(not r['strongest_clean_eqtl_in_shared'] for r in summary_rows),
        validation=validation)
    (run / 'analysis_summary.json').write_text(json.dumps(result, indent=2) + '\n')
    prefix = '../' + run.relative_to(ROOT).as_posix() + '/'
    esc = html.escape
    table_rows = []
    for r in summary_rows:
        values = [r['gene_symbol'], r['tissue'], f"{r['PP_H4_baseline']:.4f}",
                  f"{r['PP_H3_baseline']:.4f}", f"{r['PP_H4_min']:.4f}–{r['PP_H4_max']:.4f}",
                  labels[r['category']], r['n_shared'],
                  f"{float(r['fraction_eqtl_shared']):.1%}", r['missing_region_leads'] or '없음',
                  '포함' if r['strongest_clean_gwas_in_shared'] else '누락']
        table_rows.append('<tr>' + ''.join('<td>' + esc(str(v)) + '</td>' for v in values)
                          + f'<td><a href="{esc(prefix + r["plot"])}">구간 그림</a></td></tr>')
    html_text = f'''<!doctype html><html lang="ko"><meta charset="utf-8">
<title>실제 CRC GWAS–eQTL 공위치 분석</title><style>
body{{font-family:system-ui,sans-serif;margin:36px;color:#243342;line-height:1.6}}
table{{border-collapse:collapse;font-size:13px;width:100%}}th,td{{border-bottom:1px solid #ddd;padding:8px;text-align:left}}
th{{position:sticky;top:0;background:#eaf1f4}}a{{color:#087e8b}}.note{{background:#fff5dd;padding:18px}}input{{padding:10px;width:340px}}
</style><h1>실제 CRC GWAS–eQTL 공위치 분석</h1>
<p>100개 유전자 · 144개 유전자–조직 조합 · 432회 분석. coloc {esc(run_info['coloc_version'])}.</p>
<p>기본 p12=10⁻⁵에서 PP.H4 ≥ 0.8인 조합: <b>{categories['robust_across_tested_priors'] + categories['baseline_only_sensitive']}</b>.
시험한 p12=10⁻⁶, 10⁻⁵, 10⁻⁴ 모두에서 기준을 충족: <b>{categories['robust_across_tested_priors']}개 조합 · {result['robust_genes']}개 유전자</b>.</p>
<p>전체 432개 posterior를 별도 Python 계산으로 검증했습니다. 최대 절대 오차 {max_error:.2g}.</p>
<div class="note">이 결과는 거리로 정한 구간과 단일 원인 변이 가정 아래의 탐색적 결과입니다.
PP.H4는 공유 원인 변이 가설의 사후확률이며, 발현의 인과적 매개나 원인 유전자 확정을 뜻하지 않습니다.
조직 간 표본 공유 때문에 두 조직 관찰은 독립 재현이 아닙니다.
모호한 lead SNP 4개와 indel이 제외됐고, SHROOM2·GPR143의 X염색체 3개 조합은 GWAS 부재로 분석하지 못했습니다.
구간별 공유 비율, 제외된 lead, 복수 신호와 ancestry/LD 차이를 후속 검토해야 합니다.</div>
<p>p1=p2=10⁻⁴ 요청. {len(adjusted)}개 조합에서 SNP 수에 따른 패키지 자동 prior 조정이 있었으며 실제 적용값을 원본 결과에 기록했습니다.</p>
<p>각 정제 구간에서 p값이 가장 작은 GWAS 변이가 공통 표에 없는 조합은 {result['groups_missing_strongest_clean_gwas']}개,
eQTL 변이가 없는 조합은 {result['groups_missing_strongest_clean_eqtl']}개입니다. 아래 최강 GWAS 열과 구간 그림을 함께 확인하세요.
원본 GWAS에는 매우 큰 beta 추정치도 있어, 후속 후보 해석 시 원 연구의 빈도·추정 품질을 검토해야 합니다.</p>
<p><a href="{prefix}coloc_sensitivity.tsv">전체 H0–H4·사전확률 결과</a> ·
<a href="{prefix}group_summary.tsv">조합별 요약</a> · <a href="{prefix}gene_summary.tsv">유전자별 요약</a> ·
<a href="{prefix}skipped_groups.tsv">분석 보류</a> · <a href="{prefix}warnings.tsv">패키지 경고</a> ·
<a href="{prefix}posterior_validation.json">독립 검증</a></p>
<p>기본=요청 p12 10⁻⁵. PP.H4 범위는 시험한 세 prior의 최솟값–최댓값. 최소 PP.H4 내림차순.</p>
<input id="filter" placeholder="유전자·조직 검색" aria-label="결과 검색"><table><thead><tr>
<th>유전자</th><th>조직</th><th>기본 PP.H4</th><th>기본 PP.H3</th><th>PP.H4 범위</th><th>판정</th>
<th>공통 SNP</th><th>eQTL 공유 비율</th><th>공통 표에 없는 구간 lead</th><th>최강 GWAS 변이</th><th>그림</th></tr></thead>
<tbody>{''.join(table_rows)}</tbody></table>
<p>방법: <a href="https://chr1swallace.github.io/coloc/articles/a02_data.html">coloc 입력 정의</a> ·
<a href="https://chr1swallace.github.io/coloc/articles/a04_sensitivity.html">사전확률 민감도</a>.
실행 경로: {esc(run.relative_to(ROOT).as_posix())}</p>
<script>document.getElementById('filter').addEventListener('input',function(){{
const q=this.value.toLowerCase();document.querySelectorAll('tbody tr').forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(q));}});</script></html>'''
    report = ROOT / 'reports/coloc_report.html'
    report.write_text(html_text, encoding='utf-8')
    paths = sorted(p for p in run.iterdir() if p.is_file() and p.name != 'output_sha256.tsv')
    paths += [report, Path(__file__)]
    write_tsv(run / 'output_sha256.tsv', [dict(path=p.relative_to(ROOT).as_posix(), sha256=sha(p)) for p in paths])
    print(json.dumps(result, indent=2))
    print('Candidates meeting all tested priors:', ', '.join(r['gene_symbol'] for r in gene_rows if r['robust_tissues']))


if __name__ == '__main__':
    main()
