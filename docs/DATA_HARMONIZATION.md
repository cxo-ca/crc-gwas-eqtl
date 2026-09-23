# GWAS와 eQTL 기준 맞추기

확보한 구간 요약통계를 같은 GRCh38 위치·대립유전자 기준으로 연결합니다. 원본 확보 자료는 보존하고 결과는 `data/processed/harmonization/`에 따로 저장합니다. 이 정렬 단계 이후 완료한 실제 분석은 [coloc 결과](COLOC_RESULTS.md)에 별도로 기록했습니다.

## 2026-09-23 결과

|항목|결과|
|---|---:|
|정제 GWAS|255,354행|
|정제 eQTL|693,951행|
|양쪽 자료를 정렬한 공통 통계|630,083행(유전자·조직별 반복 포함), 고유 SNP 111,807개|
|최소 공통 SNP 50개를 충족한 대상|100개 유전자·144개 조합|
|조합별 공통 SNP 수|1,024–11,421개|
|정제 eQTL 중 GWAS와 공통인 비율|74.8–98.0%|
|정제 GWAS 중 각 eQTL과 공통인 비율|25.1–60.4%|
|추가 beta 부호 반전|0행: 공통 변이의 효과 대립유전자가 이미 일치|

GWAS에서 A/T·C/G SNP 43,282행을 제외했습니다. eQTL에서는 중복 rsID 별칭 799행을 합치고, A/T·C/G SNP 117,448행과 비단일염기 변이 64,204행을 제외했습니다. SE가 없는 1행은 A/T·C/G 제외와 겹치며, 총 eQTL 제외 행은 182,451행입니다. 통계값이 충돌하는 중복은 없었습니다. 정제 eQTL 158행은 정제 GWAS의 같은 위치에 다른 대립유전자 쌍만 있어 연결하지 않았습니다.

다운로드 기준 lead 41개 중 36개가 GWAS QC를 통과했습니다. **rs2581817, rs209142, rs9271363, rs7398375**는 A/T·C/G SNP여서 제외됐으며, rs2732875는 X염색체 GWAS 부재로 보류됐습니다. 공통 SNP 수를 충족한 조합에서도 lead 소실과 관측 범위 제한을 다음 분석에서 검토해야 합니다.

부호 반전·잘못된 build/좌표·중복 충돌 등 **7개 단위 테스트를 통과**했습니다. 기존 R 코드로 **147개 대상과 공통 통계 630,083행 전체를 독립 재계산해 일치**를 확인했습니다. 공통 표에도 비유의 통계(p > 0.05)가 GWAS 406,433행, eQTL 419,660행 남아 있습니다(유전자·조직별 반복 포함).

## 처리 기준

1. 확보 manifest의 SHA-256을 검사한 뒤, 두 자료의 GRCh38·1-based 좌표와 다운로드 구간 경계를 확인합니다. GTEx의 `chr_position_ref_alt` ID가 각 행의 좌표·REF·ALT와 일치하는지도 검사합니다. 유전자 버전 ID와 실측 발현 SD를 대조합니다.
2. 변이 키는 `염색체:위치:알파벳순_대립유전자1:대립유전자2`입니다. 마지막 두 값은 REF/ALT 순서를 뜻하지 않습니다. rsID만으로 합치지 않으며, 같은 위치라도 대립유전자 쌍이 다르면 연결하지 않습니다.
3. 같은 구간·유전자·조직·변이의 중복에서 rsID 이외의 모든 값이 같을 때만 한 행으로 합치고 rsID 별칭을 보존합니다. 통계값 등이 다른 중복은 어느 쪽도 임의 선택하지 않고 해당 키의 모든 행을 제외합니다.
4. 유한한 beta, 양의 유한한 SE와 분산, 0 < p ≤ 1, 양의 sdY를 확인합니다. 기존 정책대로 A/T·C/G SNP, 비단일염기 변이, 동일 대립유전자, harmonisation 실패 코드를 제외합니다. **p값에 따른 유의성 필터는 적용하지 않습니다.**
5. 공통 변이는 GWAS의 효과 대립유전자에 맞춥니다. eQTL의 효과/비효과 대립유전자가 반대면 beta에 −1을 곱하고, 원래 효과값과 대립유전자 및 반전 여부를 함께 저장합니다. SE·p·sdY는 바꾸지 않습니다.

GWAS는 제공자가 이미 GRCh38로 표준화한 파일을 사용합니다. `hm_code`는 제공자가 수행한 처리를 나타내므로 코드 11–13을 보고 다시 뒤집지 않습니다. [GWAS Catalog 처리 설명](https://ebispot.github.io/gwas-sumstats-harmoniser-documentation/Introduction/Harmonising-the-variants/)과 [코드 정의](https://ebispot.github.io/gwas-sumstats-harmoniser-documentation/Reference-guide/Hm_code/)를 따릅니다. GTEx ALT가 효과 대립유전자라는 정의 및 여러 rsID가 중복 행으로 저장된다는 설명은 [원본 열 정의](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/Columns.md)에 있습니다.

## 산출물

|파일|내용|
|---|---|
|[gwas.cleaned.tsv.gz](../data/processed/harmonization/gwas.cleaned.tsv.gz)|QC를 통과한 GWAS 전체 행. eQTL에 없는 변이도 보존|
|[eqtl.cleaned.tsv.gz](../data/processed/harmonization/eqtl.cleaned.tsv.gz)|QC를 통과한 eQTL 전체 행. **효과값은 원래 GTEx ALT 기준**이며 GWAS 미제공 X 구간도 보존|
|[shared.aligned.tsv.gz](../data/processed/harmonization/shared.aligned.tsv.gz)|공통 변이의 양쪽 통계, GWAS 기준 `beta_eqtl_aligned`, 원래 beta와 대립유전자, 반전 여부|
|[group_qc.tsv](../data/processed/harmonization/group_qc.tsv)|147개 대상별 원래/정제/공통 행 수, 공유 비율, 좌표 범위, 10등분 창 중 관측된 구간 수, lead 보존 여부|
|[excluded_rows.tsv.gz](../data/processed/harmonization/excluded_rows.tsv.gz)|제외한 행의 원본 행 번호와 사유. `source_row`는 헤더를 제외한 원본 통합 TSV의 1-based 행 번호|
|[unmatched_eqtl.tsv.gz](../data/processed/harmonization/unmatched_eqtl.tsv.gz)|정제 eQTL 중 공통 변이에 포함되지 않은 행: GWAS 구간 미제공, 같은 위치의 allele 불일치, 정제 GWAS에 변이 없음|
|[lead_variant_qc.tsv](../data/processed/harmonization/lead_variant_qc.tsv)|다운로드 창에 연결된 lead 변이의 GWAS QC 상태|
|[harmonization_summary.json](../data/processed/harmonization/harmonization_summary.json)|최종 수량·정책·입출력 SHA-256. 전체 작업이 성공한 경우에만 생성|
|[r_validation.json](../data/processed/harmonization/r_validation.json)|기존 R의 입력 검증·대립유전자 정렬 함수와 전체 공통 결과를 독립 대조한 기록|

한 행에 제외 사유가 둘 이상일 수 있어 사유별 합계와 제외 행 수는 다를 수 있습니다. 중복 제거를 먼저 수행하므로 중복 사유가 붙은 행에 allele 사유를 중복 집계하지 않습니다. `absent_from_clean_gwas`는 원본에 없거나 GWAS QC에서 제외된 경우를 포함하며, 원본 부재만 뜻하지 않습니다.

## 실행과 검증

추가 다운로드 없이 실행합니다.

```powershell
py scripts/harmonize_regional_sumstats.py
py -m unittest discover -s tests -p test_harmonize_regional_sumstats.py -v
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla scripts/verify_harmonization.R
```

정제 GWAS/eQTL 표는 기존 R `validate_stats` 입력 규격을 따릅니다. 공통 표는 이미 맞춘 효과값을 검토하기 위한 산출물입니다. R 검증에서는 원래 ALT 기준의 정제 eQTL부터 `harmonize_pair`를 다시 실행하여 공유 키·행 수·부호 반전·beta·SE·p를 모두 대조합니다.

## 해석 범위와 남은 작업

- 각 자료의 GRCh38 표기·제공자 표준화·상호 allele 일치를 검증했습니다. 별도 reference FASTA 대조, strand complement, liftover, indel 정규화는 수행하지 않았습니다.
- 최소 공통 SNP 50개는 기존 실행 기준이며 충분한 coverage의 보증이 아닙니다. `fraction_eqtl_shared`, `fraction_gwas_shared`, 좌표 범위와 lead 소실을 함께 검토해야 합니다. 10등분 창 지표도 관측 범위를 요약할 뿐 LD를 측정하지 않습니다.
- 다운로드 창은 lead ±500 kb를 합친 임시 범위입니다. LD·재조합 경계와 복수 신호는 별도 검토가 필요합니다. 두 자료의 ancestry별 allele frequency/LD는 이번 단계에서 맞추지 않았습니다.
- SHROOM2·GPR143의 X염색체 GWAS 부재는 유지됩니다. 모호한 SNP와 indel을 제외하므로 원인 후보가 빠질 수 있으며, lead QC 표에 기록합니다.
- 이 단계의 완료 범위는 **자료 기준 정렬과 QC**입니다. 후속 실제 coloc 실행기는 정제 파일을 직접 읽고 메타데이터를 각 실행 폴더에 저장합니다. [실제 분석 결과](COLOC_RESULTS.md)를 참고하세요.
