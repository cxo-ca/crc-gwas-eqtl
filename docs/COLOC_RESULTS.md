# 실제 자료 공위치 분석

GCST90129505 대장암 GWAS와 GTEx v8 대장 두 조직의 정제 구간 요약통계로 `coloc::coloc.abf`를 실행합니다. 1차 유의 변이 중첩 후보표와 별도의 실제 자료 분석입니다.

이 실행을 보존한 채 후속 [원본 누락 점검](PREVIOUS_STAGE_AUDIT.md)과 [후보 재정리·최종 보고서](../reports/final_report.html)를 완료했습니다. 137개 peak 누락 중 131개는 사용한 eQTL 원본에 위치가 없었고, 6개는 해당 유전자의 cis 범위 밖이었습니다. 점검 범위에서 추출·계산 누락은 발견되지 않았으며 누락 통계를 임의로 복구하거나 posterior를 변경하지 않았습니다.

## 2026-09-23 완료 결과

실행 폴더: [`20260923_094957_real`](../results/coloc/20260923_094957_real/analysis_summary.json).

|항목|결과|
|---|---:|
|실제 분석|100개 유전자·144개 조합|
|사전확률 3개 설정의 전체 실행|432회|
|기본 p12에서 PP.H4 ≥ 0.8|18개 유전자·23개 조합|
|시험한 세 p12 모두 PP.H4 ≥ 0.8|9개 유전자·11개 조합|
|기본 기준 충족하지만 작은 p12에서 미충족|12개 조합|
|큰 p12에서만 기준 충족|8개 조합|
|시험한 범위에서 기준 미충족|113개 조합|
|자료 부재로 보류|2개 유전자·3개 조합|

세 설정 모두 기준을 충족한 유전자는 **ACTR1B, LRP1, COLCA2, LMOD1, ATF1, COLCA1, LINC02044, MZF1-AS1, RP11-351I24.1**입니다. ACTR1B와 ATF1은 두 조직에서 충족했습니다. 이는 시험한 prior 범위에서의 지지이며 독립 재현이나 원인 유전자 확정은 아닙니다.

**공유 범위에 중요한 제한이 있습니다.** 144개 조합 중 137개에서 정제 구간의 최저 p값 GWAS 변이가 공통 표에 없었고, 4개에서는 최저 p값 eQTL 변이가 없었습니다. 구간별 해당 변이 키·p값·포함 여부를 `group_summary.tsv`에 기록했습니다. 일부 원본 GWAS에는 매우 큰 beta도 존재합니다(예: rs149255017, beta −30.5238, SE 0.827156). 임의의 효과크기 기준으로 삭제하지 않았으며, 후속 후보 평가에서는 원 연구의 빈도·추정 품질과 누락 변이를 검토해야 합니다. PP.H4 지지는 현재 공통 SNP 집합에 조건부입니다.

기존 R 테스트를 통과했고, 실제 **432개 H0–H4 posterior 모두를 Python으로 독립 재계산**해 확인했습니다. 최대 절대 차이는 약 **5.0×10⁻¹⁶**입니다. 입력·출력 SHA-256과 144개 구간 그림을 저장했습니다. 계산 재현성 검증이 생물학적 타당성 검증을 대신하지는 않습니다.

29개 조합에서 패키지의 SNP 수에 따른 prior 조정이 있었습니다. `warnings.tsv`에는 prior 조정 경고 203건과 약한 연관 경고 15건(5개 조합×3개 설정)이 기록됐습니다. 실행 로그의 반복적인 행 이름 경고는 경고표의 그룹 행을 재활용할 때 발생하는 R 표시 경고이며, 그룹 열과 posterior 대조는 모두 통과했습니다.

## 방법과 실행 범위

- GRCh38의 36개 상염색체 구간, 100개 유전자·144개 유전자–조직 조합을 사용합니다. 공통 SNP 통계는 630,083행이며 유전자·조직별 반복을 포함합니다.
- GWAS는 case-control의 log odds beta와 SE², eQTL은 원래 slope·SE²와 공식 정규화 발현행렬에서 측정한 sdY를 사용합니다. N/MAF로 효과크기나 SE를 추정하지 않습니다.
- 단일 구간 내 형질별 최대 하나의 원인 변이를 가정하는 ABF 분석입니다. LD를 입력하지 않는 현재 방법으로 다중 신호를 분리하지 않습니다.
- p1=p2=10⁻⁴를 요청하고, p12=10⁻⁶, 10⁻⁵, 10⁻⁴의 3개 설정을 시험합니다. 기본 결과는 요청 p12=10⁻⁵입니다. 표의 `p1/p2/p12`는 **실제 적용된 값**, `*_requested`는 요청한 값입니다. 설치된 coloc 5.2.3은 `SNP 수 × prior ≥ 1`이면 해당 prior를 `1/(SNP 수 + 1)`로 조정합니다.
- PP.H4 ≥ 0.8을 운영 기준으로 사용합니다. 세 설정 모두 충족하는 경우와 사전확률에 민감한 경우를 구분하며, 유전자 수준 다중검정 FDR로 해석하지 않습니다.

[공식 입력 설명](https://chr1swallace.github.io/coloc/articles/a02_data.html), [사전확률 민감도 설명](https://chr1swallace.github.io/coloc/articles/a04_sensitivity.html)을 참고했습니다. 설치된 5.2.3의 실제 함수와 기본 효과크기 사전분포(GWAS SD 0.2, eQTL SD 0.15×sdY)를 확인하고 버전 변경 없이 실행했습니다.

## 결과 파일

- [검색 가능한 HTML 결과 보고서](../reports/coloc_report.html)
- 실제 실행 폴더는 `results/coloc/latest_real_run.txt`에 기록합니다. 실행마다 별도의 시각 폴더를 생성합니다.
- `coloc_sensitivity.tsv`: 432개 설정별 PP.H0–H4, 요청/실제 사전확률.
- `group_summary.tsv`, `gene_summary.tsv`: 조합별·유전자별 요약. 서로 다른 조합의 posterior를 합치지 않습니다.
- `qc.tsv`, `preanalysis_qc.tsv`, `lead_variant_qc.tsv`, `skipped_groups.tsv`: 공통 변이 수·공유 비율·lead 소실·보류 사유.
- `locus_*.png`: 각 조합의 GWAS/eQTL 구간 그림. 회색은 QC 통과 후 한쪽 자료에만 있는 변이, 색 점은 공통 변이입니다.
- `warnings.tsv`: 약한 연관 및 사전확률 자동 조정 등 패키지 경고를 조합/요청 p12별로 기록합니다.
- `harmonized.tsv`: 실제 계산에 전달한 공통 통계.
- `analysis_settings.rds`, `metadata.R`, `sessionInfo.txt`, `input_sha256.tsv`, `output_sha256.tsv`: 출처·설정·버전·재현성 기록.
- `posterior_validation.json`: 정렬 입력에서 별도 Python ABF 계산으로 전체 432개 posterior를 대조한 기록.
- `analysis_summary.json`: 결과 정리와 독립 검증까지 성공했을 때의 완료 요약.

## 재실행

추가 다운로드가 필요 없습니다.

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla scripts/05_run_real_coloc.R
py scripts/report_real_coloc.py
```

새 실행기는 이미 검증한 정제 압축 파일을 직접 읽고 모든 확보/정렬 체크섬을 확인합니다. `data/input/`에 대용량 복사본을 만들지 않습니다. `scripts/04_run_coloc.R`은 기존 사용자 제공 입력용 경로로 유지합니다.

## 해석 제한

현재 구간은 lead ±500 kb 창을 합친 탐색적 거리 구간입니다. LD·재조합 경계를 확정하지 않았고, 복수 신호와 GWAS/eQTL ancestry 차이를 모델링하지 않았습니다. QC 통과 GWAS의 25.1–60.4%, eQTL의 74.8–98.0%가 조합별 공통 SNP에 해당합니다. 제공된 cis 검사 범위 밖의 변이를 추가하거나 누락 효과를 추정하지 않았습니다.

A/T·C/G lead SNP **rs2581817, rs209142, rs9271363, rs7398375** 및 indel이 제외됐습니다. SHROOM2 두 조직과 GPR143 Transverse는 X염색체 GWAS 부재로 보류합니다. 정렬 과정의 모든 제외 이유는 [QC 기록](DATA_HARMONIZATION.md)에 남아 있습니다.

높은 PP.H4는 현재 자료·구간·가정에서 공유 원인 변이 가설을 지지한다는 뜻입니다. 발현이 질병 위험을 매개한다거나 특정 유전자가 원인이라는 확정은 아닙니다. 낮은 PP.H4도 공통 변이 누락이나 약한 신호의 영향을 받을 수 있습니다. 두 조직은 표본을 공유하므로 독립 재현으로 세지 않습니다. LD/복수 신호 검토와 독립·기능 검증은 후속 단계입니다.
