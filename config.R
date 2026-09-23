config <- list(
  study = 'GCST90129505',
  gtex_dataset = 'gtex_v8',
  tissues = c('Colon_Transverse', 'Colon_Sigmoid'),
  gwas_p = 5e-8,
  genome_build = 'GRCh38',
  min_shared_snps = 50L,
  pp4_threshold = 0.8,
  p1 = 1e-4, p2 = 1e-4, p12 = c(1e-6, 1e-5, 1e-4),
  seed = 20260922L
)
