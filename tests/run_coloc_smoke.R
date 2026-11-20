# Actual R integration on preserved synthetic data, with output confined to tempdir().
source('tests/coloc/environment.R')
check_coloc_test_environment()
source('R/utils.R'); source('R/coloc_pipeline.R'); source('config.R')
input <- read_tsv('results/demo/harmonized.tsv')
trait_input <- function(kind) {
  common <- c('locus', 'chr', 'pos', 'build')
  trait <- c('effect_allele', 'other_allele', 'beta', 'se', 'p')
  d <- input[c(common, paste0(trait, '_', kind))]
  names(d) <- c(common, trait)
  if (kind == 'eqtl') d <- cbind(d, input[c('gene_id', 'tissue', 'sdY')])
  unique(d)
}
meta <- list(full_region_gwas = TRUE, full_region_eqtl = TRUE, allele_strand = 'forward',
             gwas_source = 'PRESERVED SYNTHETIC DEMO', eqtl_source = 'PRESERVED SYNTHETIC DEMO')
outdir <- tempfile('coloc-smoke-')
actual <- run_coloc(trait_input('gwas'), trait_input('eqtl'), meta, config, outdir)
expected <- read_tsv('results/demo/coloc_sensitivity.tsv')
stopifnot(nrow(actual) == 6L, nrow(expected) == 6L)
keys <- c('locus', 'gene_id', 'tissue', 'p1', 'p2', 'p12')
stopifnot(!anyDuplicated(actual[keys]), !anyDuplicated(expected[keys]))
compared <- merge(actual, expected, by = keys, suffixes = c('_actual', '_expected'))
stopifnot(nrow(compared) == 6L, all(compared$nsnps_actual == 240L))
for (hypothesis in paste0('PP_H', 0:4)) {
  observed <- compared[[paste0(hypothesis, '_actual')]]
  reference <- compared[[paste0(hypothesis, '_expected')]]
  stopifnot(all(is.finite(observed)), all(observed >= 0 & observed <= 1),
            all(abs(observed - reference) <= 1e-10 * reference))
}
stopifnot(all(abs(rowSums(actual[paste0('PP_H', 0:4)]) - 1) < 1e-12))
qc <- read_tsv(file.path(outdir, 'qc.tsv'))
group_keys <- c('locus', 'gene_id', 'tissue')
stopifnot(nrow(qc) == 2L, all(qc$n_shared == 240L), all(qc$n_flipped == 120L),
          all(qc$status == 'analyzed'), !anyDuplicated(qc[group_keys]),
          nrow(merge(qc[group_keys], unique(expected[group_keys]), by = group_keys)) == 2L)
aligned <- read_tsv(file.path(outdir, 'harmonized.tsv'))
row_keys <- c(group_keys, 'key')
stopifnot(nrow(aligned) == 480L, !anyDuplicated(aligned[row_keys]),
          !anyDuplicated(input[row_keys]))
aligned <- merge(aligned, input, by = row_keys, suffixes = c('_actual', '_expected'))
stopifnot(nrow(aligned) == 480L,
          isTRUE(all.equal(aligned$beta_eqtl_aligned_actual,
                           aligned$beta_eqtl_aligned_expected, tolerance = 1e-12)),
          identical(aligned$allele_flipped_actual, aligned$allele_flipped_expected))
cat('PASS: coloc 5.2.3 integration; six fits, 30 posteriors, allele flips and QC\n')
