# 전체 구간 요약통계 확보

## 2026-09-23 확보 결과

공개 경로에서 제공되는 대상 자료 수집을 완료했습니다. 이 문서와 통합 자료는 확보 당시 상태를 기록합니다. 이후 완료한 QC·대립유전자 정렬은 [기준 맞추기 결과](DATA_HARMONIZATION.md)에, 실제 공위치 분석은 [coloc 결과](COLOC_RESULTS.md)에 별도로 기록했습니다.

|항목|확보 결과|
|---|---:|
|GWAS 구간|36개 상염색체 구간|
|GWAS 행|298,636|
|eQTL 구간|61개 조직별 구간(37개 유전체 창)|
|eQTL 행|876,402|
|GWAS와 eQTL 모두 확보한 유전자|100개|
|양쪽 자료가 있는 유전자–조직 조합|144개|
|eQTL만 포함한 전체 대상|102개 유전자·147개 조합|

- [GWAS 통합 자료](../data/processed/acquisition/gwas.regional.tsv.gz)
- [eQTL 통합 자료](../data/processed/acquisition/eqtl.regional.tsv.gz)
- [최종 확보 요약](../data/processed/acquisition/acquisition_summary.json)
- [QC 수량](../data/processed/acquisition/acquisition_qc.tsv)
- [기존 API 대조 결과](../data/processed/acquisition/original_api_validation.json): 기존 153개 연관의 p값과 효과크기 모두 일치(상대 허용오차 10⁻⁴).
- [오프라인 재현 검증](../data/processed/acquisition/offline_validation.json): 두 통합 TSV와 발현 SD 표의 내용 SHA-256이 재실행 전후 동일. gzip 시간 정보는 비교에서 제외했습니다.

비유의 행(p > 0.05)도 GWAS 244,133행, eQTL 589,686행 포함됩니다. eQTL에서 동일 변이–유전자 키 중복 799행과 beta/SE 검사 문제 1행을 발견했으며, 원본 그대로 보존하고 QC 표에 기록했습니다. 후속 QC에서 검토해야 합니다. X염색체 GWAS가 없는 SHROOM2와 GPR143는 아래 제한 사항대로 별도 표시했습니다.

## 범위와 실행

기존 공개 후보 102개에 연결된 41개 lead 변이를 기준으로 GRCh38 좌표의 ±500 kb 창을 만들고 겹치는 창을 합쳤습니다. 총 37개 구간, 147개 유전자–조직 조합입니다. 이 범위는 **다운로드를 위한 임시 창**이며 LD나 재조합 경계를 검증한 최종 분석 구간이 아닙니다.

유의성 필터 없이 해당 구간의 GWAS 변이 전체와 대상 유전자의 eQTL 변이 전체를 추출합니다. 다른 유전자의 eQTL은 최종 추출 표에 포함하지 않습니다. 원래 GTEx cis 검사 범위 밖의 변이를 새로 만들어 채우지 않습니다.

```powershell
py scripts/fetch_gtex_expression.py
py scripts/acquire_regional_sumstats.py
```

Python 3 표준 라이브러리와 `curl.exe`를 사용합니다. 이미 검증된 구간과 바이트 캐시는 재사용합니다. 실패한 전송을 빈 결과로 대체하지 않습니다. 각 원본의 HTTP 범위 응답, 압축 블록 CRC, SHA-256, 파일 크기와 응답 버전 정보를 검사합니다. 동시 네트워크 요청은 최대 4개입니다.

긴 전송을 위한 작업 실행기는 `py scripts/run_acquisition_job.py`입니다. 현재 Windows에서 숨김 백그라운드 프로세스로 실행하도록 구성했으며, 진행 상태는 `data/processed/acquisition/download_job_status.json`, 시도별 로그는 `logs/download_job_*.log`에 기록합니다. 네트워크 오류만 최대 12회, 전체 최대 8시간 동안 재시도합니다. 다른 데이터 검증 오류는 중단합니다. `running`과 `retrying_network`는 진행 중이고, `completed_available_data`는 알려진 자료 미제공 구간을 제외한 확보 완료입니다. 컴퓨터가 꺼지거나 절전 상태이면 전송은 진행되지 않습니다. 작업 중복 실행은 파일 잠금으로 차단합니다.

자료 확보가 끝난 후 네트워크 없이 재생성할 수 있습니다.

```powershell
py scripts/fetch_gtex_expression.py --offline
py scripts/acquire_regional_sumstats.py --offline
```

## 원본 출처

### GWAS

- 연구: Fernandez-Rozadilla et al., Nature Genetics (2023), `GCST90129505`.
- 표본: 100,204 cases / 154,587 controls; 유럽·동아시아 메타분석.
- [GRCh38 harmonised 파일과 색인](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90129001-GCST90130000/GCST90129505/harmonised/).
- [연구 README와 academic license](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90129001-GCST90130000/GCST90129505/README.txt).
- 제공된 MD5와 비교한 메타데이터: `716307e778ac79d018db01f4982001b4`.
- 제공된 MD5와 비교한 Tabix 색인: `fc94254cca9f1c834a4513de254d9840`.

전체 압축 파일은 640,903,823 bytes입니다. 최초 전체 전송은 연결 종료로 중단되어 `.part`만 남았으며 **완료된 원본으로 취급하지 않습니다**. 이후 검증된 색인으로 필요한 구간만 받습니다. 따라서 전체 파일 MD5 검증을 완료했다고 주장하지 않으며, 구간별 바이트 체크섬과 BGZF CRC를 기록합니다.

이 파일에는 academic/noncommercial 이용 조건이 있으므로 원문 README와 YAML을 원본 폴더에 보존했습니다.

### GTEx v8 eQTL

- [eQTL Catalogue의 imported studies 설명](https://www.ebi.ac.uk/eqtl/Studies/).
- [GTEx v8 원본 imported 요약통계](https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/): `Colon_Sigmoid.tsv.gz`, `Colon_Transverse.tsv.gz` 및 각 `.tbi`.
- [열 정의](https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/Columns.md): GRCh38, ALT가 효과 대립유전자, `beta`, `se`, `pvalue` 사용.
- [공식 imported GTEx 분석 예제](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tutorials/tabix_use_case.md).

`GTEx_V8` imported 자료는 GTEx v8 원본 결과의 재배포입니다. Catalogue가 별도로 재분석한 `GTEx` 연구로 바꾸지 않았습니다. 표준오차는 원본 `se`를 사용하며 p값에서 추정하지 않습니다. 원본의 결측 SE와 다중 rsID에 의한 중복을 보존하고 수량을 보고합니다. 이를 제거하거나 allele을 정렬하는 것은 후속 QC 단계입니다.

### 발현형질 표준편차

[공식 GTEx 정규화 발현행렬 TAR](https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL_expression_matrices.tar)에서 다음 두 멤버만 받습니다.

- `Colon_Sigmoid.v8.normalized_expression.bed.gz`: 318 samples.
- `Colon_Transverse.v8.normalized_expression.bed.gz`: 368 samples.

TAR 전체 1,134,960,640 bytes 중 두 발현 멤버의 합계는 50,440,989 bytes입니다. 공개된 정규화 발현값으로 유전자별 sample SD(N−1 분모)를 직접 계산합니다. `sdY=1`을 임의로 넣지 않습니다. [Broad 정규화 설명](https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/README.md).

Sigmoid에서는 전체 후보 102개 중 99개의 발현 행이 있고 Transverse에서는 102개 모두 있습니다. 기존 탐색에서 실제 관찰된 **147개 유전자–조직 조합에는 모두 발현 행이 있습니다**. 발현 파일과 후보표의 버전 있는 유전자 ID가 다르면 중단합니다.

## 알려진 범위 제한

GWAS 표준화 파일의 검증된 색인에는 1–22 염색체만 있습니다. 다음 X 구간은 이 경로에서 제공되지 않습니다.

- 구간: `chrX_9295858_10295858`, lead `rs2732875`.
- 영향: `SHROOM2` 두 조직과 `GPR143` Transverse, 총 2개 유전자·3개 조합.

eQTL은 X 구간도 수집하지만 GWAS를 추정하거나 다른 연구로 몰래 대체하지 않습니다. 현재 공개 경로로 양쪽 자료를 확보할 수 있는 목표는 **36개 상염색체 구간, 100개 유전자, 144개 유전자–조직 조합**입니다. X 구간의 GWAS는 별도 원자료 확보가 필요합니다.

## 산출물과 완료 판정

|경로|내용|
|---|---|
|`data/processed/acquisition/regions.tsv`|구간 정의와 연결된 lead 목록|
|`data/processed/acquisition/targets.tsv`|대상 유전자·조직·구간|
|`data/raw/summary_statistics/range_cache/`|원본 압축 바이트, 색인, HTTP 헤더, 체크섬|
|`data/raw/summary_statistics/regions/`|구간별 원본 열 추출본과 출처 JSON|
|`data/raw/summary_statistics/GTEx_v8/expression/`|두 조직 발현 원본과 TAR 범위 manifest|
|`data/processed/acquisition/gtex_expression_sdY.tsv`|실측 발현 표준편차와 출처|
|`data/processed/acquisition/gwas.regional.tsv.gz`|후속 QC용 GWAS 표|
|`data/processed/acquisition/eqtl.regional.tsv.gz`|후속 QC용 eQTL 표|
|`data/processed/acquisition/acquisition_qc.tsv`|행 수, 비유의 변이, 결측·중복 수|
|`data/processed/acquisition/source_availability.tsv`|출처별 구간 확보 또는 명시적 미제공 상태|
|`data/processed/acquisition/target_coverage.tsv`|대상별 eQTL 행 수와 GWAS 제공 여부|
|`data/processed/acquisition/acquisition_summary.json`|모든 예정 파일 검증 후 생성되는 완료 요약|

스크립트와 다운로드 계획만 있다고 완료된 것은 아닙니다. 완료 요약과 대상 coverage를 확인해야 합니다. 출력은 **후속 QC용 자료**이며 `data/input/metadata.R`의 전체 구간·allele 검증 완료 플래그를 자동으로 켜지 않습니다. 공위치 분석은 이 단계에서 실행하지 않습니다.
