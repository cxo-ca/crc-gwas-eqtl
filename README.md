# 대장암 GWAS·eQTL 기반 유전자 후보 탐색

R로 공개 대장암 GWAS 변이와 대장 조직 eQTL을 연결하고, 구간 전체 요약통계가 준비되면 공위치 분석(colocalization)으로 후보 근거를 평가하는 프로젝트입니다.

## 현재 완료된 분석

2026-09-23에 **앞 단계 누락 점검과 후보 재정리·최종 보고서 작성**까지 완료했습니다. 점검한 누락은 사용한 원본의 부재·cis 범위·의도적 QC 제외로 확인됐으며, 확인 범위에서 처리 과정 누락이나 복구 행은 없었습니다. 원래 후보 102개를 모두 유지해 prior 지지가 안정적인 조건부 후보 9개, 기본 prior 지지·민감 후보 9개, 큰 prior에서만 지지 7개, 기준 미충족 75개, 자료 부재 보류 2개로 정리했습니다.

- **[최종 보고서](reports/final_report.html)** · [Markdown](reports/final_report.md)
- **[처음 보는 사람을 위한 6쪽 소개 자료](reports/introduction/research_intro.pdf)** · [웹 문서](reports/introduction/research_intro.html) · [3분 발표 설명문](reports/introduction/research_intro_script.md)
- [전 단계 누락 점검](docs/PREVIOUS_STAGE_AUDIT.md)
- [전체 후보 102개](results/final/candidates.tsv) · [조합별 147개 결과·효과 방향·제한](results/final/candidate_groups.tsv)

현재 확보한 자료에 대한 탐색적 분석의 최종 정리입니다. LD·복수 신호·원본 추정 품질·독립/기능 검증은 남아 있으며 원인 유전자 확정을 뜻하지 않습니다.

2026-09-22에 **실제 공개 자료의 1차 후보 탐색**을 실행했습니다.

|항목|결과|
|---|---:|
|GWAS Catalog의 유의한 단일 rsID (p ≤ 5×10⁻⁸)|174|
|GTEx v8에 매핑된 rsID|169|
|매핑되지 않은 rsID|5|
|변이–유전자–조직 연결|153|
|후보 유전자|102|
|두 대장 조직에서 관찰된 후보|45|

- [HTML 분석 보고서](reports/public_report.html)
- [전체 후보 유전자 표](results/public/candidate_genes.tsv)
- [변이별 근거 표](results/public/variant_gene_evidence.tsv)
- [품질 점검 요약](results/public/qc_summary.tsv)
- [원본 자료 체크섬](data/processed/download_manifest.tsv)

이 결과는 **GWAS 변이와 유의한 eQTL이 겹치는 후보**입니다. 인과 유전자 확정이나 실제 자료의 공위치 분석 결과가 아닙니다. 논문의 전체 연관 신호 수와 Catalog에서 이번에 수집한 레코드 수는 다를 수 있습니다. 본 분석은 지정한 연구의 현재 Catalog 등록 자료를 사용합니다.

## 후속 분석용 자료 확보

2026-09-23에 후속 분석용 **비유의 변이를 포함한 구간 요약통계**도 확보했습니다. GWAS 298,636행, eQTL 876,402행이며, **100개 유전자·144개 유전자–조직 조합**에 양쪽 자료가 있습니다. X염색체의 두 후보는 해당 GWAS 파일에 통계가 없어 별도 표시했습니다. 기존 153개 연관 대조와 오프라인 재현 검증을 통과했습니다. [확보 자료·검증 결과](docs/DATA_ACQUISITION.md)를 확인하세요.

이어 **좌표·변이 ID·대립유전자 정렬과 중복/결측 QC**를 수행해 공통 통계 630,083행을 만들었습니다. 100개 유전자·144개 조합 모두 최소 공통 SNP 수를 충족했습니다. 제외한 lead SNP 4개와 X 구간 부재, 구간별 공유 비율은 [기준 정렬 결과](docs/DATA_HARMONIZATION.md)에 기록했습니다.

## 실제 자료 공위치 분석 완료

2026-09-23에 **100개 유전자·144개 조합에 대해 사전확률 3개 설정, 총 432회 coloc**을 실행하고 별도 계산으로 전체 posterior를 검증했습니다. 기본 p12에서 18개 유전자·23개 조합이 PP.H4 ≥ 0.8을 충족했고, 시험한 세 설정 모두에서 기준을 유지한 것은 **9개 유전자·11개 조합**입니다.

[실제 coloc 결과 보고서](reports/coloc_report.html) · [방법·후보·검증·해석 제한](docs/COLOC_RESULTS.md)

현재 공통 SNP 집합과 거리 구간에 조건부인 탐색적 결과입니다. 137개 조합에서 정제 GWAS의 최강 변이가 공통 표에 없으므로 coverage·원본 추정 품질·복수 신호를 후속 검토해야 합니다. 원인 유전자 확정 결과는 아닙니다.

## RStudio에서 시작

`crc-gwas-eqtl.Rproj`를 열고 프로젝트 루트에서 실행합니다.

```r
source("scripts/setup.R")
# 저장된 자료를 재분석: 네트워크가 필요하지 않습니다.
source("scripts/02_analyze_public.R")
```

Windows PowerShell에서 전체 흐름을 실행할 수도 있습니다.

공개 자료의 수집·분석을 한 번에 재실행하려면 `Rscript run.R --offline`을 사용합니다.

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla scripts/setup.R --coloc
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla scripts/01_fetch_public.R --offline
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla scripts/02_analyze_public.R
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla scripts/03_run_demo.R
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla tests/run_tests.R
```

다른 컴퓨터에서는 설치된 `Rscript` 경로로 바꾸세요. Rscript가 PATH에 등록되어 있으면 `Rscript scripts/02_analyze_public.R`처럼 실행합니다. `--offline` 없이 수집 스크립트를 실행하면 **없는 캐시만** 내려받습니다. 기존 캐시는 덮어쓰지 않습니다. API 실패를 빈 결과로 처리하지 않습니다. 재분석 결과표는 같은 이름으로 갱신됩니다.

기본 설정은 [config.R](config.R)에 있습니다. 패키지는 `.Rlib/`에 설치하며 이미 사용할 수 있는 패키지는 기존 R 라이브러리에서 이용합니다. 실제 버전은 `logs/package_versions.csv`와 각 단계의 `sessionInfo`에 기록합니다. 이는 실행 환경 기록이며 자동 버전 고정용 lockfile은 아닙니다. 한글 사용자 경로를 처리하기 위해 Windows에서는 bootstrap이 UTF-8 locale을 설정합니다.

## 분석 구성

1. **공개 자료 수집**: `GCST90129505`의 등록 연관 변이를 수집하고 p ≤ 5×10⁻⁸로 선별합니다. GTEx 공식 API로 rsID를 GRCh38 ref/alt 변이에 매핑합니다.
2. **조직별 후보 탐색**: GTEx v8 `Colon_Transverse`, `Colon_Sigmoid`의 유의한 cis-eQTL과 직접 연결합니다. 버전을 명시해 API 기본 버전이 바뀌어도 혼합하지 않습니다. 모든 페이지를 수집합니다.
3. **후보 정리**: 두 조직 관찰 여부 → 최소 eQTL p값 → 최소 GWAS p값 순으로 정렬합니다. 순서는 검토 편의를 위한 것이며 통계적 인과확률이나 검증된 종합점수가 아닙니다.
4. **공위치 분석 준비**: `coloc::coloc.abf`를 사용하며 H0–H4와 p12 사전확률 민감도를 저장합니다. 합성 자료에서 공통 피크와 서로 다른 피크를 구별하는지 검증합니다.

기존 후보 탐색 결과인 `results/public/`은 유의한 lead 변이와 유의 eQTL만 포함하므로 **coloc에 바로 넣을 수 없습니다**. `data/processed/acquisition/`에는 확보 원본을, `data/processed/harmonization/`에는 QC·정렬 결과를 저장합니다. 두 단계 모두 p값에 따른 유의성 선별을 하지 않았습니다. 공위치 분석 전에는 임시 구간과 coverage, 제외된 lead를 검토해야 합니다.

## 실제 자료의 공위치 분석

이미 확보·정렬한 실제 자료는 아래 명령으로 분석합니다. 결과와 해석 범위는 [실제 분석 문서](docs/COLOC_RESULTS.md)에 있습니다.

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla scripts/05_run_real_coloc.R
py scripts/report_real_coloc.py
```

별도의 사용자 제공 자료를 분석하려면 아래 입력 규격과 기존 실행기를 사용합니다.

[입력 자료 설명](docs/DATA_DICTIONARY.md)에 맞춰 아래 파일을 준비합니다.

```text
data/input/gwas.tsv
data/input/eqtl.tsv
data/input/metadata.R
```

`metadata.example.R`를 참고하여 자료 출처, 전체 구간 자료 여부, 대립유전자 방향 등을 기록한 뒤 실행합니다.

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla scripts/04_run_coloc.R
```

결과는 `results/coloc/<실행시각>/`에 저장됩니다. `results/demo/`는 모두 **합성 자료 결과**이며 실제 후보표에 포함하지 않습니다.

## 폴더

|경로|내용|
|---|---|
|`R/`|API·입출력·QC·대립유전자 정렬·coloc 함수|
|`scripts/`|설치, 수집, 실제 공개 자료 분석, 예제, 사용자 자료 분석|
|`data/raw/`|공개 API 원본 JSON, URL·수집시각|
|`data/processed/`|정리한 공개 자료, 매핑, 다운로드 체크섬|
|`data/input/`|사용자 구간 요약통계 입력|
|`results/public/`|실제 공개 자료 후보표·그림|
|`results/demo/`|합성 자료 공위치 분석·민감도·QC|
|`reports/`|한국어 HTML 보고서|
|`docs/`|연구 설계·자료 사전·출처|
|`tests/`|분석 오류 방지 점검|

자세한 해석과 후속 연구는 [연구 설계](docs/ANALYSIS_PLAN.md), 원문 및 데이터 링크는 [자료 출처](docs/SOURCES.md)에 정리했습니다.

## GitHub 공개 범위

저장소에는 분석 코드, 설정, 테스트, 문서, 보고서, 발표 자료, 최종 후보표와 재현에 필요한 기록을 포함합니다. 공개 원자료를 내려받아 만든 `data/raw/summary_statistics/` 캐시와 GitHub의 단일 파일 제한을 넘는 `results/coloc/*/harmonized.tsv`는 커밋하지 않습니다. 원자료 URL·체크섬·수집 방법은 문서와 manifest에 남아 있고, 같은 정렬 입력의 압축본은 `data/processed/harmonization/shared.aligned.tsv.gz`에 보존합니다. 자세한 목록은 [GitHub 게시 기록](docs/GITHUB_PUBLICATION.md)에 있습니다.
