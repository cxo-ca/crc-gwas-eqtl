source('R/bootstrap.R'); source('R/utils.R'); source('config.R')
collected <- readRDS('data/processed/public_collection_config.rds')
for (field in c('study', 'gtex_dataset', 'tissues', 'gwas_p')) {
  if (!identical(config[[field]], collected[[field]])) {
    stop('Configuration differs from collected data; rerun 01_fetch_public.R: ', field)
  }
}
g <- read_tsv('data/processed/gwas_lead_variants.tsv')
m <- read_tsv('data/processed/variant_mapping.tsv')
e <- read_tsv('data/processed/gtex_significant_eqtl.tsv')
joined <- merge(merge(g, m, by = 'rsid'), e, by = 'variant_id')
joined <- joined[order(joined$eqtl_p, joined$gwas_p), ]
joined$evidence <- rep('GWAS lead variant overlaps significant eQTL; not colocalization', nrow(joined))
write_tsv(joined, 'results/public/variant_gene_evidence.tsv')
candidates <- bind_rows(lapply(split(joined, joined$gene_id), function(d) {
  best <- d[which.min(d$eqtl_p), ]
  data.frame(gene_id = best$gene_id, gene_symbol = best$gene_symbol,
    n_tissues = length(unique(d$tissue)), n_lead_variants = length(unique(d$rsid)),
    tissues = paste(sort(unique(d$tissue)), collapse = ';'),
    min_gwas_p = min(d$gwas_p), min_eqtl_p = min(d$eqtl_p),
    best_eqtl_variant = best$rsid,
    evidence = 'overlap_only')
}))
if (!nrow(candidates)) candidates <- data.frame(gene_id = character(),
  gene_symbol = character(), n_tissues = integer(), n_lead_variants = integer(),
  tissues = character(), min_gwas_p = numeric(), min_eqtl_p = numeric(),
  best_eqtl_variant = character(), evidence = character())
# Sorting is a transparent review order, not a calibrated probability/score.
candidates <- candidates[order(-candidates$n_tissues, candidates$min_eqtl_p,
                               candidates$min_gwas_p, candidates$gene_id), ]
write_tsv(candidates, 'results/public/candidate_genes.tsv')
qc <- data.frame(metric = c('significant_catalog_rsids', 'mapped_rsids',
  'unmapped_rsids', 'variant_gene_tissue_links', 'candidate_genes', 'genes_in_both_tissues'),
  value = c(nrow(g), length(unique(m$rsid)), length(setdiff(g$rsid, m$rsid)),
    nrow(joined), nrow(candidates), sum(candidates$n_tissues == length(config$tissues))))
write_tsv(qc, 'results/public/qc_summary.tsv')
if (nrow(candidates)) {
  top <- head(candidates, 20)
  png('results/public/candidate_evidence.png', width = 1600, height = 1100, res = 150)
  par(mar = c(5, 9, 4, 2))
  barplot(rev(-log10(pmax(top$min_eqtl_p, .Machine$double.xmin))),
    names.arg = rev(top$gene_symbol), horiz = TRUE, las = 1,
    col = rev(ifelse(top$n_tissues == 2, '#087E8B', '#A6CDD1')),
    border = NA, xlab = '-log10(minimum nominal eQTL p-value)',
    main = 'CRC GWAS / colon eQTL overlap candidates')
  legend('bottomright', c('Both colon tissues', 'One colon tissue'),
    fill = c('#087E8B', '#A6CDD1'), bty = 'n')
  dev.off()
}
report <- c('<!doctype html><html lang="ko"><meta charset="utf-8">',
 '<title>대장암 GWAS·eQTL 후보 탐색</title>',
 '<style>body{font-family:Segoe UI,Malgun Gothic,sans-serif;max-width:1120px;margin:48px auto;padding:0 24px;color:#16313b;line-height:1.7;background:#f7fafb}h1{font-size:32px}h2{margin-top:36px}.note{padding:20px;background:#e6f2f3;border-left:5px solid #087e8b}table{border-collapse:collapse;width:100%;font-size:13px;background:white}td,th{text-align:left;padding:10px;border-bottom:1px solid #d8e4e8}th{background:#eaf1f4}img{max-width:100%}a{color:#087e8b}</style>',
 '<h1>대장암 GWAS·eQTL 기반 유전자 후보 탐색</h1>',
 paste0('<p>연구 ', config$study, ' · ', config$gtex_dataset,
   ' · Colon Transverse / Sigmoid · 생성: ', Sys.Date(), '</p>'),
 '<div class="note"><strong>실제 공개 자료를 이용한 1차 후보 탐색</strong><br>GWAS 유의 변이와 GTEx 유의 eQTL의 직접 중첩 결과입니다. 인과 유전자 확정, 발현의 질병 인과효과, 공위치 분석 결과를 뜻하지 않습니다.</div>',
 '<h2>분석 범위와 품질 점검</h2>', html_table(qc),
 '<p>GWAS p ≤ 5×10⁻⁸인 단일 rsID를 GTEx GRCh38 변이로 매핑했습니다. LD 대리 변이는 추가하지 않았습니다. 다중 대립유전자 rsID는 복수 매핑이 가능하므로 variant_gene_evidence.tsv의 ref/alt를 후속 확인해야 합니다. eQTL이 조회되지 않았다는 사실은 발현 조절 효과가 없다는 증거가 아닙니다.</p>',
 '<h2>후보 유전자</h2>',
 '<p>두 조직에서 관찰된 유전자, 최소 eQTL p값, 최소 GWAS p값 순으로 정렬했습니다. 이 순서는 검토 편의를 위한 것이며 독립적인 재현이나 인과확률 순위가 아닙니다. 두 조직의 기증자는 겹칠 수 있습니다. 아래에는 상위 30개를 표시합니다.</p>',
 html_table(candidates), if (nrow(candidates)) '<img src="../results/public/candidate_evidence.png" alt="후보별 eQTL 근거">',
 '<h2>후속 분석</h2><p>같은 genome build의 구간 전체 GWAS 및 유전자별 eQTL beta/SE 자료를 준비해 공위치 분석을 수행합니다. 유의한 SNP만으로 PP.H4를 계산하지 않습니다. 다중 인과변이, 조상집단 간 LD 차이, 건강한 대장 조직과 종양의 차이를 고려해야 합니다. NES는 GTEx 대체 대립유전자 기준이며 이 보고서에서는 대장암 위험 방향과 결합하지 않았습니다.</p>',
 '<h2>자료와 재현</h2><p><a href="../results/public/candidate_genes.tsv">전체 후보표</a> · <a href="../results/public/variant_gene_evidence.tsv">변이별 근거</a> · <a href="../data/processed/download_manifest.tsv">원본 체크섬</a></p>',
 '<p><a href="https://www.ebi.ac.uk/gwas/studies/GCST90129505">GWAS Catalog</a> · <a href="https://doi.org/10.1038/s41588-022-01222-9">원 연구</a> · <a href="https://gtexportal.org/api/v2/openapi.json">GTEx 공식 API</a> · <a href="https://chr1swallace.github.io/coloc/articles/a02_data.html">coloc 자료 요건</a></p></html>')
writeLines(enc2utf8(report), 'reports/public_report.html', useBytes = TRUE)
writeLines(capture.output(sessionInfo()), 'logs/public_analysis_sessionInfo.txt')
saveRDS(config, 'results/public/analysis_settings.rds')
manifest(c('data/processed/gwas_lead_variants.tsv', 'data/processed/variant_mapping.tsv',
  'data/processed/gtex_significant_eqtl.tsv', 'config.R'), 'results/public/input_manifest.tsv')
print(qc)
