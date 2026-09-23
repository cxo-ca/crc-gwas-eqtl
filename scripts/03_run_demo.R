source('R/bootstrap.R'); source('R/utils.R'); source('R/coloc_pipeline.R'); source('config.R')
set.seed(config$seed)
n <- 240L
positions <- seq(10000001, by = 1000, length.out = n)
# Artificial regional patterns: one shared peak and one distinct peak.
# Synthetic test data only; these are not human gene or disease findings.
zg <- 0.6 + 9 * exp(-((seq_len(n) - 80)/15)^2) + rnorm(n, 0, 0.15)
g <- data.frame(locus = 'SIMULATED_LOCUS', chr = '1', pos = positions,
  effect_allele = 'A', other_allele = 'C', beta = zg * 0.025, se = 0.025,
  p = 2 * pnorm(-abs(zg)), build = 'GRCh38')
make_eqtl <- function(center, gene) {
  z <- 0.5 + 10 * exp(-((seq_len(n) - center)/15)^2) + rnorm(n, 0, 0.15)
  d <- g; d$beta <- z * 0.06; d$se <- 0.06; d$p <- 2 * pnorm(-abs(z))
  d$gene_id <- gene; d$tissue <- 'SIMULATED_COLON'; d$sdY <- 1
  flip <- seq(2, n, by = 2)
  d$effect_allele[flip] <- 'C'; d$other_allele[flip] <- 'A'; d$beta[flip] <- -d$beta[flip]
  d
}
e <- rbind(make_eqtl(80, 'SIM_SHARED'), make_eqtl(180, 'SIM_DISTINCT'))
write_tsv(g, 'data/processed/demo_gwas.tsv'); write_tsv(e, 'data/processed/demo_eqtl.tsv')
meta <- list(full_region_gwas = TRUE, full_region_eqtl = TRUE,
  allele_strand = 'forward', gwas_source = 'SYNTHETIC', eqtl_source = 'SYNTHETIC')
r <- run_coloc(g, e, meta, config, 'results/demo')
writeLines(c('SYNTHETIC DEMONSTRATION ONLY. Not colorectal cancer findings.',
  'SIM_SHARED and SIM_DISTINCT are artificial identifiers.',
  'The profiles test shared versus distinct peaks and allele swaps, not biological calibration.'),
  'results/demo/README.txt')
stopifnot(all(r$PP_H4[r$gene_id == 'SIM_SHARED'] > 0.8),
          all(r$PP_H3[r$gene_id == 'SIM_DISTINCT'] > 0.8))
print(r[c('gene_id', 'p12', 'PP_H3', 'PP_H4')])
