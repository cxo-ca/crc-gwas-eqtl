# 후보 재정리와 최종 보고서

[최종 HTML 보고서](../reports/final_report.html)와 [Markdown 보고서](../reports/final_report.md)에 누락 점검·방법·후보 분류·효과 방향·해석 제한을 정리했습니다.

원래 102개 후보를 모두 보존했습니다. 세 prior에서 지지되는 조건부 후보 9개, 기본 prior 지지·민감 후보 9개, 큰 prior에서만 지지 7개, 기준 미충족 75개, 자료 부재 보류 2개입니다. 대상 조합 147개도 모두 남겨 미분석 3개를 음성 결과와 구분합니다.

실행 순서는 다음과 같습니다. 추가 네트워크 전송이나 coloc 재계산을 하지 않습니다.

```powershell
py scripts/audit_previous_stages.py
py scripts/build_final_report.py
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla scripts/plot_final_candidates.R
py scripts/verify_final_report.py
```

`results/final/final_validation.json`은 후보 보존·기존 posterior 일치·효과 방향 부호·원본 해시·보고서 링크·PNG/PDF 형식을 확인합니다. 검증 후 `final_summary.json`에 완료 상태를 기록하고 `final_manifest.tsv`에 산출물과 생성 코드의 SHA-256을 저장합니다. 그림의 범위는 prior별 PP.H4 범위이며 신뢰구간이 아닙니다.

최종 보고서는 현재 자료의 탐색적 결과를 정리한 산출물입니다. 자료 미제공, LD/복수 신호, 원본 추정 품질, reference/indel 검토, 독립·기능 검증은 완료된 것으로 표시하지 않습니다. 세부 사유는 [점검 기록](PREVIOUS_STAGE_AUDIT.md)과 `results/final/issue_ledger.tsv`에 있습니다.
