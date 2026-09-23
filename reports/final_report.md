# 대장암 GWAS·대장 eQTL 후보 재정리 최종 보고서

2026-09-23 | 현재 확보 자료에 대한 탐색적 분석의 최종 정리

## 먼저 확인한 전 단계 누락

체크섬 181건과 원본 색인의 좌표 조회 96건을 대조했습니다. 확인 대상에서 추출·계산 과정의 누락은 발견되지 않았으며 복구된 행은 0개입니다. 모든 미공유 SNP를 전장 유전체 수준에서 점검했다는 뜻은 아닙니다.

- 최강 GWAS 변이가 빠진 137개 조합: eQTL 원본에 위치 없음 131개, 유전자 TSS ±1 Mb cis 범위 밖 6개.
- 최강 eQTL 변이가 빠진 4개 조합: 사용한 harmonised GWAS 원본에 위치 없음.
- 모호한 lead 4개(rs2581817, rs209142, rs9271363, rs7398375)는 의도적 QC 제외. 직접 연결된 36개 발견 조합에 영향.
- SHROOM2 두 조직·GPR143 Transverse의 3개 조합은 GWAS X 자료 부재로 보류.
- 점검한 137개 GWAS peak 행의 효과 대립유전자 빈도(EAF)는 모두 NA였습니다. 빈도 대조로 누락 원인이나 모호한 SNP 방향을 해결하지 못했습니다.
- LD·복수 신호·독립 reference FASTA·indel 복구·외부 및 기능 검증은 미수행 항목으로 남겼습니다.

## 자료·방법·검증

GCST90129505의 GRCh38 대장암 GWAS와 GTEx v8 Colon_Sigmoid/Colon_Transverse를 사용했습니다. 1차 중첩 탐색에서 선정한 102개 후보 중 100개·144개 조합에 실제 coloc을 수행했습니다. 다른 암 아형이나 종양 특이 효과를 따로 분석하지 않았습니다.

lead ±500 kb 병합 거리 창, beta/SE² 및 실측 sdY, coloc 5.2.3의 형질별 단일 원인 변이 가정을 사용했습니다. p1=p2=10⁻⁴, 요청 p12=10⁻⁶·10⁻⁵·10⁻⁴로 432회 실행했습니다. 29개 조합의 패키지 prior 자동 조정은 실제 적용값과 함께 원래 결과에 남아 있습니다. p값으로 변이를 선별하지 않았습니다.

432개 posterior 모두 별도 Python 계산과 일치했습니다(최대 절대 차이 약 5×10⁻¹⁶). 수치 재현성 검증은 원본 품질이나 생물학적 타당성을 보증하지 않습니다.

## 후보 분류

원래 후보 102개와 대상 147개 조합을 모두 유지했습니다. 유전자별 분류는 가장 강한 지지를 보이는 조합을 요약하며 조직별 posterior를 합산하지 않습니다. 같은 등급에서는 대표 조합의 최소 PP.H4 순서입니다. 확률이나 별도 종합점수를 새로 추정한 순위가 아닙니다.

|분류|유전자|조합|
|---|---:|---:|
|A: 세 prior 지지·자료 제한 조건부|9|11|
|B: 기본 prior 지지·민감|9|12|
|C: 큰 prior에서만 지지|7|8|
|D: 시험 범위 기준 미충족|75|113|
|E: 자료 부재·보류|2|3|

기본 prior에서 기준을 충족한 18개 유전자 중 9개는 작은 prior에서 지지가 낮아집니다. A 등급 9개 유전자도 공통 변이 집합과 구간 가정에 조건부인 검토 후보입니다.

## 세 prior에서 기준을 유지한 후보

|유전자|지지 조직|대표 기본 PP.H4|대표 최소 PP.H4|대표 SNP|연관 방향|
|---|---|---:|---:|---|---|
|ACTR1B|Colon_Sigmoid;Colon_Transverse|0.9996|0.9963|rs11692435|발현 증가 방향과 CRC odds 감소 방향 일치|
|LRP1|Colon_Transverse|0.9977|0.9771|rs715948|발현 증가 방향과 CRC odds 증가 방향 일치|
|COLCA2|Colon_Transverse|0.9971|0.9714|rs3087967|발현 증가 방향과 CRC odds 감소 방향 일치|
|LMOD1|Colon_Sigmoid|0.9920|0.9254|rs4494199|발현 증가 방향과 CRC odds 감소 방향 일치|
|ATF1|Colon_Sigmoid;Colon_Transverse|0.9846|0.8644|rs11169571|발현 증가 방향과 CRC odds 증가 방향 일치|
|COLCA1|Colon_Transverse|0.9834|0.8555|rs3087967|발현 증가 방향과 CRC odds 감소 방향 일치|
|LINC02044|Colon_Transverse|0.9794|0.8263|rs13086367|발현 증가 방향과 CRC odds 증가 방향 일치|
|MZF1-AS1|Colon_Sigmoid|0.9763|0.8049|rs73068329|발현 증가 방향과 CRC odds 감소 방향 일치|
|RP11-351I24.1|Colon_Transverse|0.9759|0.8020|rs57153895|발현 증가 방향과 CRC odds 증가 방향 일치|

방향은 대표 조합에서 H4 조건부 SNP 확률이 가장 큰 변이의 같은 효과 대립유전자 기준 beta 부호를 비교한 것입니다. 대표 변이의 SNP.PP.H4는 H4가 참이라는 조건부 값이며, 해당 변이의 인과 확률이나 유전자 전체의 매개 효과가 아닙니다. 조직별 방향과 원래 beta는 전체 표에 보존했습니다.

## 해석과 다음 검증

자료 누락을 해결했다고 간주해 후보를 확정하지 않았습니다. 특히 137개 조합의 최강 GWAS 변이 부재와 원본의 극단 beta 추정치, 제외된 lead, 거리 구간/복수 신호가 후속 검토 대상입니다. A 등급부터 이 항목을 확인한 뒤 ancestry에 맞는 LD·조건부 분석, 독립 대장 eQTL 및 기능 검증으로 이어가야 합니다. 두 조직은 기증자를 공유하므로 독립 재현이 아닙니다.

D 등급은 시험한 모델에서 기준 미충족이며 생물학적 관련성이 없다는 뜻이 아닙니다. E 등급은 미분석이며 음성 결과로 취급하지 않습니다. 유의 eQTL 중첩으로 후보를 선정했으므로 전체 유전자 FDR·enrichment·약물 표적 확정으로 확대 해석하지 않습니다.

## 파일과 출처

- [전체 102개 후보 TSV](../results/final/candidates.tsv)
- [전체 147개 조합·방향·제한 TSV](../results/final/candidate_groups.tsv)
- [누락 원본 대조 TSV](../results/final/missing_peak_audit.tsv)
- [제외 lead 대조 TSV](../results/final/excluded_lead_audit.tsv)
- [미해결 항목과 처리 방침](../results/final/issue_ledger.tsv)
- [coloc 실행 결과와 그림](coloc_report.html)
- [세부 coloc 방법](../docs/COLOC_RESULTS.md)
- [민감도 그림 PNG](final_candidate_sensitivity.png) · [PDF](final_candidate_sensitivity.pdf)
- [GWAS 원본](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90129001-GCST90130000/GCST90129505/harmonised/)
- [GTEx imported 원본](https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/)
- [GTEx cis 1 Mb 설정](https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/README.md)
- [coloc prior 민감도](https://chr1swallace.github.io/coloc/articles/a04_sensitivity.html)
