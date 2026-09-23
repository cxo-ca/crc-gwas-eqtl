# 자료 출처

확인 및 최초 수집: 2026-09-22. 원본 응답은 `data/raw/GCST90129505/gtex_v8/`, 요청 URL과 실제 UTC 수집시각은 각 `.source.txt`에 기록합니다. 다운로드 원본의 MD5는 `data/processed/download_manifest.tsv`에 저장합니다.

2026-09-23에 추가 확보한 전체 구간 GWAS·GTEx v8 imported eQTL 및 정규화 발현행렬의 출처와 이용 조건은 [전체 구간 자료 확보](DATA_ACQUISITION.md)에 기록했습니다. 원본과 HTTP 범위별 증거는 `data/raw/summary_statistics/`, 통합 TSV·검증 결과는 `data/processed/acquisition/`에 있습니다.

|자료|원 출처|이 프로젝트에서의 사용|
|---|---|---|
|대장암 GWAS|[Fernandez-Rozadilla et al., Nature Genetics 55, 89–99 (2023)](https://doi.org/10.1038/s41588-022-01222-9)|기준 연구. 100,204 cases / 154,587 controls의 유럽·동아시아 연구|
|GWAS Catalog 연구|[GCST90129505](https://www.ebi.ac.uk/gwas/studies/GCST90129505)|등록된 단일 변이 연관성의 탐색; full summary statistics와 구분|
|Catalog API|[연구 메타데이터](https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90129505), [연관 변이](https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90129505/associations?size=500)|연구 정보, p값, rsID, 원본 레코드 링크|
|GTEx 공개 자료|[GTEx 다운로드](https://gtexportal.org/home/downloads/adult-gtex/overview)|v8의 정상 대장 조직 eQTL. full all-pairs 자료 접근·파일 크기·비용 조건은 해당 출처에서 확인|
|GTEx 공식 API|[OpenAPI 규격](https://gtexportal.org/api/v2/openapi.json)|`dataset/variant`로 매핑, `association/singleTissueEqtl`로 유의 eQTL 조회. 두 endpoint 모두 `datasetId=gtex_v8` 지정|
|coloc 입력 요건|[Data structures](https://chr1swallace.github.io/coloc/articles/a02_data.html)|같은 지역의 조밀한 전체 변이 요약통계, beta/varbeta 및 발현형질 척도|
|coloc 모형|[Single causal variant analysis](https://chr1swallace.github.io/coloc/articles/a03_enumeration.html)|형질별 단일 인과변이 가정의 H0–H4 비교|
|coloc 민감도|[Sensitivity to prior assumptions](https://chr1swallace.github.io/coloc/articles/a04_sensitivity.html)|p12 변경에 대한 민감도 검토|

API endpoint 설명에 rsID 언급이 있어도 실제 파라미터·반환값은 규격으로 확인해야 합니다. 이 프로젝트는 rsID를 `dataset/variant?snpId=...`에서 GTEx variant ID로 변환한 뒤 해당 ID를 `singleTissueEqtl?variantId=...`에 전달합니다. API에 인식되지 않는 필터를 넣어 빈 응답을 잘못 해석하는 것을 방지합니다.
