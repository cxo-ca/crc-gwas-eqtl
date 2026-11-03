"""Render a Korean one-page research overview from preserved result tables.

Run from any directory: py scripts/plot_research_overview.py
Requires matplotlib; uses Malgun Gothic on Windows or an installed Korean font.
No analysis inputs or existing result tables are changed.
"""

import csv
import json
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Rectangle


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports" / "overview"
INK = "#183047"
MUTED = "#53677A"
TEAL = "#087F75"
GRID = "#E4EBEF"
STABLE = "A_conditional_prior_stable"
TIERS = [
    (STABLE, "세 설정 모두에서 지지", TEAL),
    ("B_baseline_prior_sensitive", "기본 설정 지지 · 설정에 민감", "#4D80AB"),
    ("C_permissive_prior_only", "큰 사전확률 설정에서만 지지", "#B98838"),
    ("D_no_threshold_support", "시험한 설정에서 기준 미충족", "#8B9BA9"),
    ("E_data_unavailable", "자료 부재로 보류", "#CAD4DB"),
]


def read_tsv(relative):
    with (ROOT / relative).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def configure_font():
    malgun = Path("C:/Windows/Fonts/malgun.ttf")
    if malgun.is_file():
        font_manager.fontManager.addfont(str(malgun))
        family = font_manager.FontProperties(fname=str(malgun)).get_name()
    else:
        available = {font.name for font in font_manager.fontManager.ttflist}
        family = next(
            (name for name in ("Malgun Gothic", "Apple SD Gothic Neo", "Noto Sans CJK KR")
             if name in available), None
        )
        if family is None:
            raise RuntimeError("Install a Korean font before rendering the overview")
    plt.rcParams.update({
        "font.family": family,
        "font.size": 10,
        "text.color": INK,
        "axes.labelcolor": MUTED,
        "xtick.color": MUTED,
        "ytick.color": INK,
        "axes.unicode_minus": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def rounded(fig, x, y, width, height, color):
    fig.add_artist(FancyBboxPatch(
        (x, y), width, height, transform=fig.transFigure,
        boxstyle="round,pad=0.007,rounding_size=0.012",
        facecolor=color, edgecolor="none", zorder=0,
    ))


def main():
    configure_font()
    summary = json.loads((ROOT / "results/final/final_summary.json").read_text(encoding="utf-8"))
    candidates = read_tsv("results/final/candidates.tsv")
    groups = read_tsv("results/final/candidate_groups.tsv")
    gene_counts = Counter(row["evidence_tier"] for row in candidates)
    group_counts = Counter(row["evidence_tier"] for row in groups)
    assert len(candidates) == summary["candidate_genes"]
    assert len(groups) == summary["candidate_groups"]
    assert gene_counts == summary["gene_tiers"]
    assert group_counts == summary["group_tiers"]
    observed = [row for row in groups if row["PP_H4_baseline"] not in ("", "NA")]
    assert len(observed) == summary["analyzed_groups"]
    stable = [row for row in candidates if row["evidence_tier"] == STABLE]
    assert [row["gene_symbol"] for row in stable] == summary["prior_stable_genes"]
    peak_reasons = Counter(row["gwas_peak_absence_reason"] for row in observed)
    absent = sum(row["strongest_clean_gwas_in_shared"].lower() == "false" for row in observed)
    assert absent == (
        peak_reasons["position_absent_from_imported_source"]
        + peak_reasons["outside_documented_cis_window"]
    )
    fits = summary["posterior_validation"]["fits_checked"]
    assert fits == len(observed) * 3

    fig = plt.figure(figsize=(15.2, 10.8), facecolor="white")
    fig.text(.045, .964, "CRC  /  GWAS × eQTL  /  RESEARCH OVERVIEW", fontsize=10,
             color=TEAL, weight="bold")
    fig.text(.045, .919, "대장암 후보 유전자, 어떤 근거로 좁혔을까?",
             fontsize=25, weight="bold")
    fig.text(.045, .882, "공개자료를 연결하고, 같은 유전 신호인지 비교하고, 가정을 바꾸어 후보를 검토한 연구",
             fontsize=11.5, color=MUTED)
    fig.text(.955, .964, "분석 결과 기준  2026.09.23", fontsize=9, ha="right", color=MUTED)

    cards = [
        ("01  공개자료 연결", "GWAS × 대장 eQTL", "GCST90129505 · GTEx v8 두 조직"),
        ("02  후보 탐색 · 자료 정렬", f"{len(candidates)}개 후보 유전자", "위치·대립유전자 정렬 · 2개 자료 부재"),
        ("03  실제 공위치 분석", f"{summary['analyzed_genes']}개 유전자 · {len(observed)}개 조합", "공통 변이에서 공유 원인 가설 평가"),
        ("04  사전확률 민감도 평가", f"{len(stable)}개 유전자 · {group_counts[STABLE]}개 조합", "세 설정 모두에서 PP.H4 ≥ 0.8 유지"),
    ]
    for index, (label, value, note) in enumerate(cards):
        x = .045 + index * .233
        rounded(fig, x, .742, .209, .101, "#EAF5F2" if index == 3 else "#F0F4F7")
        fig.text(x + .012, .819, label, fontsize=9, color=TEAL, weight="bold")
        fig.text(x + .012, .786, value, fontsize=13.3, weight="bold")
        fig.text(x + .012, .759, note, fontsize=8.1, color=MUTED)
        if index < 3:
            fig.text(x + .218, .783, "→", fontsize=16, ha="center", color=TEAL)

    fig.text(.045, .694, "A   전체 후보는 어떻게 분류되었나", fontsize=13, weight="bold")
    fig.text(.045, .668, "단위: 유전자 · 조직별 결과 중 가장 강한 지지를 기준으로 분류", fontsize=8.9, color=MUTED)
    ax = fig.add_axes([.245, .439, .225, .188])
    values = [gene_counts[key] for key, _, _ in TIERS]
    ax.barh(range(5), values, color=[color for _, _, color in TIERS], height=.48)
    ax.set_yticks(range(5), [label for _, label, _ in TIERS], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, max(values) * 1.12)
    ax.set_xticks([0, 25, 50, 75])
    ax.tick_params(axis="both", length=0, labelsize=8.5, pad=6)
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, color=GRID)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for index, value in enumerate(values):
        ax.text(value + 1.2, index, str(value), va="center", fontsize=10, weight="bold")
    fig.text(.045, .411, "기준 미충족은 생물학적 무관함을 뜻하지 않으며, 자료 부재는 미분석입니다.",
             fontsize=8.4, color=MUTED)

    fig.text(.545, .694, "B   우선 검토 후보 9개의 조직별 지지", fontsize=13, weight="bold")
    fig.text(.545, .668, "각 칸은 세 사전확률 설정 중 최소 PP.H4 · 0.8 이상이면 모두 기준 충족", fontsize=8.9, color=MUTED)
    tissues = ["Colon_Sigmoid", "Colon_Transverse"]
    lookup = {(row["gene_id"], row["tissue"]): row for row in observed}
    matrix = np.full((len(stable), len(tissues)), np.nan)
    for i, gene in enumerate(stable):
        for j, tissue in enumerate(tissues):
            row = lookup.get((gene["gene_id"], tissue))
            if row is not None:
                matrix[i, j] = float(row["PP_H4_min"])
    cmap = LinearSegmentedColormap.from_list("evidence", ["#F0F5F6", "#96CCC0", TEAL])
    cmap.set_bad("#F3F4F5")
    heat = fig.add_axes([.710, .293, .225, .327])
    picture = heat.imshow(np.ma.masked_invalid(matrix), cmap=cmap, vmin=0, vmax=1, aspect="auto")
    heat.set_xticks([0, 1], ["S상결장", "횡행결장"], fontsize=10)
    heat.xaxis.tick_top()
    heat.set_yticks(range(len(stable)), [row["gene_symbol"] for row in stable], fontsize=9)
    heat.tick_params(axis="both", length=0, pad=8)
    heat.set_xticks(np.arange(-.5, 2, 1), minor=True)
    heat.set_yticks(np.arange(-.5, len(stable), 1), minor=True)
    heat.grid(which="minor", color="white", linewidth=3)
    heat.tick_params(which="minor", length=0)
    for spine in heat.spines.values():
        spine.set_visible(False)
    for i in range(len(stable)):
        for j in range(2):
            value = matrix[i, j]
            if np.isnan(value):
                heat.add_patch(Rectangle((j-.5, i-.5), 1, 1, hatch="///", fill=False,
                                         edgecolor="#D5DCE1", linewidth=0))
                label, color = "—", MUTED
            else:
                label, color = f"{value:.3f}", "white" if value >= .7 else INK
            heat.text(j, i, label, ha="center", va="center", fontsize=9.5, color=color, weight="bold")
    cax = fig.add_axes([.710, .257, .225, .010])
    bar = fig.colorbar(picture, cax=cax, orientation="horizontal", ticks=[0, .5, 1])
    bar.outline.set_visible(False)
    bar.ax.tick_params(labelsize=7, length=0, pad=3)
    fig.text(.545, .216, "빗금 / — : 초기 후보에 해당 조합이 없어 미분석 · 0과 구분", fontsize=8.4, color=MUTED)
    fig.text(.545, .194, "두 조직은 기증자를 공유하므로 독립적인 반복 검증으로 세지 않습니다.", fontsize=8.4, color=MUTED)

    fig.text(.045, .366, "C   공통 입력 자료의 범위도 함께 확인", fontsize=13, weight="bold")
    fig.text(.045, .337, f"분석한 {len(observed)}개 조합 중 {absent}개에서 구간 최강 GWAS 변이가 공통 입력에 없음",
             fontsize=9.3, color=MUTED)
    coverage = fig.add_axes([.045, .278, .425, .035])
    present = len(observed) - absent
    coverage_values = [peak_reasons["position_absent_from_imported_source"],
                       peak_reasons["outside_documented_cis_window"], present]
    left = 0
    for value, color in zip(coverage_values, ["#CFD8E0", "#D6B375", TEAL]):
        coverage.barh(0, value, left=left, color=color, height=.8)
        left += value
    coverage.set_xlim(0, len(observed))
    coverage.set_ylim(-.5, .5)
    coverage.set_axis_off()
    fig.text(.045, .252, f"eQTL 원본에 위치 없음 {coverage_values[0]}  ·  cis 검사 범위 밖 {coverage_values[1]}  ·  공통 입력에 있음 {present}",
             fontsize=8.6, color=MUTED)
    fig.text(.045, .229, "점검 범위에서 처리 과정 누락은 발견되지 않았습니다.", fontsize=8.6, color=MUTED)
    fig.text(.045, .197, f"계산 검증  |  {summary['posterior_validation']['shared_rows_checked']:,}행 입력 · posterior {fits}건 독립 재계산 일치",
             fontsize=8.8, color=TEAL, weight="bold")

    rounded(fig, .045, .060, .910, .098, "#F0F4F7")
    fig.text(.058, .133, "현재 성과", fontsize=10, weight="bold", color=TEAL)
    fig.text(.129, .133, f"후속 검증 후보 {len(stable)}개 선정 · {len(observed)}개 조합 × 사전확률 3개 설정 = {fits}회 통계 분석",
             fontsize=10, weight="bold")
    fig.text(.058, .105, "해석 범위", fontsize=9, weight="bold")
    fig.text(.129, .105, "현재 공통 변이와 단일 원인 변이 모델에 조건부인 결과 · 원인 유전자나 발현 매개 효과의 확정은 아님",
             fontsize=9, color=MUTED)
    fig.text(.058, .078, "다음 검증", fontsize=9, weight="bold")
    fig.text(.129, .078, "원본 품질·변이 범위 보완 → LD·복수 신호 분석 → 독립 자료 재현 → 기능 실험", fontsize=9, color=MUTED)
    fig.text(.045, .029, "자료: results/final/{candidates.tsv, candidate_groups.tsv, final_summary.json} · 공개자료 기반 분석 · 기능 실험 미수행",
             fontsize=7.6, color=MUTED)

    OUT.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        fig.savefig(OUT / f"research_overview.{extension}", dpi=200, facecolor="white")
    plt.close(fig)
    print(json.dumps({
        "files": [str((OUT / f"research_overview.{ext}").relative_to(ROOT)) for ext in ("png", "pdf")],
        "candidate_genes": len(candidates), "analyzed_groups": len(observed),
        "prior_stable_genes": len(stable), "prior_stable_groups": group_counts[STABLE],
        "missing_gwas_peak_groups": absent,
        "heatmap_observed_cells": int(np.isfinite(matrix).sum()),
        "heatmap_not_targeted_cells": int(np.isnan(matrix).sum()),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
