# 입력 자료와 출력 해석

## 공위치 분석 입력

탭으로 구분한 UTF-8 TSV를 사용합니다. **동일한 GRCh38 구간 전체**의 요약통계가 필요합니다. 유의한 변이만 남긴 자료, lead SNP 목록, GTEx significant pairs만으로 분석하지 않습니다. 두 자료에 공통으로 관측된 변이가 구간을 충분히 덮는지 별도로 확인해야 합니다.

### 공통 필수 열

|열|형식|의미|
|---|---|---|
|`locus`|문자|같은 분석 구간 식별자. 유전자와 별개이며 한 염색체의 지역 구간이어야 함|
|`chr`|문자/정수|1–22, X, Y. `chr` 접두사 허용|
|`pos`|양의 정수|1-based GRCh38 좌표|
|`effect_allele`|문자|beta의 효과 대립유전자; 정방향 기준|
|`other_allele`|문자|대조 대립유전자; 정방향 기준|
|`beta`|실수|GWAS: log(OR), eQTL: 제공된 정규화 발현 척도의 slope|
|`se`|양의 실수|beta 표준오차. 내부에서 제곱해 varbeta로 전달|
|`p`|(0,1]|연관 p값. 0으로 반올림하지 않음|
|`build`|문자|`GRCh38`|

### eQTL 추가 필수 열

|열|형식|의미|
|---|---|---|
|`gene_id`|문자|유전자 ID. 자료 버전의 Ensembl/GENCODE 버전 유지 권장|
|`tissue`|문자|조직 ID|
|`sdY`|양의 실수|분석에 쓰인 발현형질의 표준편차; 유전자·조직별 상수|

`sdY=1`은 발현형질을 실제로 표준화한 경우 또는 자료 제공자의 분석 정의로 타당성을 확인한 경우만 사용합니다. GTEx의 NES가 있다는 이유로 무조건 1을 넣지 않습니다. 실제 발현 척도와 표준편차가 없으면 현재 경로로는 진행하지 않습니다. 원본 eQTL의 slope SE를 사용하며 유의 p값에서 임의로 추정하지 않습니다.

GWAS는 case-control 자료이며 beta가 log odds 척도라고 가정합니다. OR 값을 그대로 beta에 넣으면 안 됩니다. 정량형 GWAS나 hazard ratio 자료는 이 스크립트의 범위가 아닙니다. 현재 beta/SE와 명시적 eQTL sdY를 사용하는 coloc ABF 경로에서는 N/MAF를 추정하지 않습니다.

```text
locus	chr	pos	effect_allele	other_allele	beta	se	p	build
L1	1	10000001	A	C	0.10	0.02	5.733e-7	GRCh38
```

위 한 줄은 형식 예시이며 분석에 필요한 변이 범위·수에 미치지 못합니다. 실행 가능한 합성 입력은 `scripts/03_run_demo.R` 실행 후 `data/processed/demo_*.tsv`에 생성됩니다.

## QC 정책

확보 자료에 적용한 사전 QC와 정렬 결과는 [DATA_HARMONIZATION.md](DATA_HARMONIZATION.md)에 있습니다. 사전 처리기는 rsID만 다른 동일 통계 중복을 합치고 제외 사유를 기록합니다. 아래 R 입력 검증기는 처리 후에도 남은 중복·결측을 오류로 중단합니다. `*.cleaned.tsv.gz`는 각각 원래 효과 대립유전자를 유지하며, `shared.aligned.tsv.gz`의 `beta_eqtl_aligned`가 GWAS 기준 효과값입니다.

- genome build 불일치, 필수값 결측, SE ≤ 0, p 범위 오류, 동일 분석 그룹의 중복 변이는 오류로 중단합니다.
- 양쪽 자료의 forward strand를 전제로 같은 위치와 대립유전자 쌍을 연결합니다. rsID만으로 coloc 변이를 합치지 않습니다.
- 대립유전자 순서가 뒤바뀐 경우 eQTL beta 부호를 반전합니다. strand complement 및 liftover는 자동 수행하지 않습니다.
- A/T, C/G의 모호한 SNP와 indel은 현재 구현에서 제외하고 제외 수를 QC에 저장합니다. 대립유전자 빈도로 모호성을 해결하거나 indel 정규화를 하는 기능은 아직 없습니다. 원인 변이 제외 가능성을 반드시 검토합니다.
- 그룹마다 공통 SNP 50개 미만이면 건너뜁니다. 50은 실행을 위한 운영 기준이며 충분한 coverage나 통계적 타당성을 보장하지 않습니다. 공유 비율·그림·원자료 coverage를 확인합니다.
- locus는 연구자가 사전에 정의한 단일 지역 구간이어야 합니다. 한 염색체 전체를 한 locus로 넣지 않습니다. 거리 창의 예로 lead ±500 kb를 검토할 수 있지만 LD/재조합 구조에 맞게 정해야 합니다. 거리로 자르는 작업과 LD clumping은 다른 절차입니다.

## coloc 결과

|열|해석|
|---|---|
|`PP_H0`|두 형질 모두 구간 내 연관 없음|
|`PP_H1`|GWAS만 연관|
|`PP_H2`|eQTL만 연관|
|`PP_H3`|양쪽이 서로 다른 원인 변이와 연관|
|`PP_H4`|양쪽이 하나의 원인 변이를 공유|
|`p1`, `p2`, `p12`|계산에 사용한 사전확률|
|`passes_pp4`|설정한 PP.H4 운영 기준 충족 여부; 기본 0.8|

각 locus·gene·tissue에 p12 = 10⁻⁶, 10⁻⁵, 10⁻⁴ 결과를 모두 저장합니다. PP.H4 기준은 보편적인 확정 기준이 아닙니다. PP.H3와 H0–H2, 사전확률 민감도, coverage를 함께 판단합니다. 현재 `coloc.abf`는 형질별 구간 내 최대 하나의 인과변이를 가정합니다. 다중 신호에는 ancestry에 맞는 LD와 조건부 분석 또는 SuSiE 기반 확장이 필요합니다.

## 공개 후보표

`min_gwas_p`와 `min_eqtl_p`는 해당 후보에 연결된 레코드의 최솟값이며 서로 다른 변이·조직에서 나올 수 있습니다. `best_eqtl_variant`는 최소 eQTL p값의 rsID입니다. `n_lead_variants`는 rsID 수이며 독립적인 LD 신호 수가 아닙니다. GTEx의 `nes`는 대체 대립유전자 기준의 정규화 효과입니다. 이 프로젝트는 Catalog 위험 대립유전자와 NES를 결합한 질병 방향을 보고하지 않습니다.

유의 eQTL API의 필터를 통과한 결과를 이용하므로 이 선택된 표에 BH 보정을 다시 적용하여 전체 유전자 수준 FDR로 해석하지 않습니다. `n_tissues=2`도 공유 기증자 때문에 독립 재현을 뜻하지 않습니다.
