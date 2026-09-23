# 앞 단계 누락 점검

2026-09-23에 최종 후보 정리 전에 기존 분석의 입력·출력 SHA-256 181건을 대조하고, 저장된 원본 압축 바이트와 Tabix 색인으로 96개 좌표 조회를 오프라인 재실행했습니다. 대상은 보고된 최강 변이 누락, 제외된 lead, X염색체 자료 부재입니다. 전장 유전체의 모든 미공유 변이를 점검한 것은 아닙니다.

|점검 대상|확인 결과|처리|
|---|---|---|
|최강 GWAS 변이 누락 137개 조합|131개는 GTEx imported 원본에 위치 없음; 6개는 해당 유전자의 TSS ±1 Mb 밖|추출 코드의 누락으로 확인되지 않음. 가상 통계로 채우지 않음|
|최강 eQTL 변이 누락 4개 조합|사용한 harmonised GWAS 원본에 위치 없음|자료 범위 제한으로 유지|
|제외된 lead 4개|원본에서도 A/T·C/G, 의도적 QC 제외; 36개 직접 발견 조합에 영향|방향·빈도 근거 없이 복구하지 않음|
|SHROOM2·GPR143|검증된 GWAS 색인에 X/23 없음; 3개 조합 미분석|음성 결과로 분류하지 않음|
|효과 대립유전자 빈도|점검한 GWAS peak 137행의 EAF는 모두 NA|빈도 대조 미수행 상태를 명시|
|입력·출력 무결성|기록된 SHA-256 모두 일치|기존 coloc 결과 보존|

확인 범위에서 **처리 과정 누락 0건, 복구한 행 0개**입니다. 원본에 왜 특정 변이가 없는지(검정 대상·유전형 QC·빈도 필터 등)는 원자료 없이 단정하지 않습니다. 6개 cis 범위 판정은 [GTEx 공식 ±1 Mb 설정](https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/README.md)과 확보한 유전자 TSS를 대조했습니다.

일부 GWAS 원본 beta가 매우 큰 점도 남은 품질 검토 항목입니다. LD/재조합 경계·복수 신호·ancestry 차이·reference FASTA 대조·indel/모호한 SNP 복구·외부/기능 검증은 완료된 것으로 처리하지 않습니다. 102개 후보 자체가 유의 eQTL 중첩으로 선택됐다는 제한도 최종 보고서에 포함했습니다.

근거: [변이별 원본 대조](../results/final/missing_peak_audit.tsv), [제외 lead 대조](../results/final/excluded_lead_audit.tsv), [점검 요약](../results/final/audit_summary.json), [남은 검토 항목](../results/final/issue_ledger.tsv).

재실행: `py scripts/audit_previous_stages.py`. 기존 coloc 입력·출력 체크섬이 달라지거나 원본에 존재하는 누락 행을 찾으면 중단합니다.
