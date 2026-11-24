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

# All-ineligible runs retain diagnostic QC in a fresh directory, including no-GWAS groups.
g <- trait_input('gwas')
e <- trait_input('eqtl')
ineligible <- e[e$pos %in% head(g$pos, config$min_shared_snps - 1L), ]
no_gwas <- head(e, 10L); no_gwas$locus <- 'SYNTHETIC_NO_GWAS'
ineligible <- rbind(ineligible, no_gwas)
must_fail <- function(expr, message) {
  err <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(err, 'error'), grepl(message, conditionMessage(err), fixed = TRUE))
}
fresh <- tempfile('coloc-ineligible-')
must_fail(run_coloc(g, ineligible, meta, config, fresh), 'No eligible coloc groups.')
skipped <- read_tsv(file.path(fresh, 'qc.tsv'))
stopifnot(nrow(skipped) == 3L, !anyDuplicated(skipped[group_keys]),
          nrow(merge(skipped[group_keys], unique(ineligible[group_keys]), by = group_keys)) == 3L,
          all(skipped$status == 'insufficient_shared_snps'),
          all(skipped$n_shared[skipped$locus == 'SIMULATED_LOCUS'] == 49L),
          sum(skipped$locus == 'SYNTHETIC_NO_GWAS') == 1L)
missing <- skipped[skipped$locus == 'SYNTHETIC_NO_GWAS', ]
stopifnot(missing$n_gwas == 0L, missing$n_shared == 0L, missing$n_eqtl == 10L,
          missing$fraction_gwas_shared == 0, missing$fraction_eqtl_shared == 0,
          nrow(read_tsv(file.path(fresh, 'warnings.tsv'))) == 0L,
          setequal(list.files(fresh), c('qc.tsv', 'warnings.tsv')))

# Reusing a successful directory must not pair old posteriors with newly skipped QC.
paths <- list.files(outdir, full.names = TRUE)
before <- tools::md5sum(paths)
must_fail(run_coloc(g, ineligible, meta, config, outdir), 'Existing output preserved')
stopifnot(identical(before, tools::md5sum(paths)),
          identical(paths, list.files(outdir, full.names = TRUE)))
cat('PASS: all-ineligible QC, missing GWAS and preservation of existing output\n')

# Invalid priors must be rejected before changed inputs can overwrite prior plots.
changed_g <- g; changed_g$p <- sqrt(changed_g$p)
paths <- list.files(outdir, full.names = TRUE)
before <- tools::md5sum(paths)
for (priors in list(NULL, numeric(), -1e-5, Inf, '1e-5')) {
  bad_config <- config; bad_config$p12 <- priors
  must_fail(run_coloc(changed_g, e, meta, bad_config, outdir), 'Invalid config$p12')
  stopifnot(identical(before, tools::md5sum(paths)),
            identical(paths, list.files(outdir, full.names = TRUE)))
}

# Preserve coloc's supported zero prior and SNP-count-based automatic adjustment.
one_gene <- e[e$gene_id == 'SIM_SHARED', ]
for (prior in c(0, 2)) {
  compatible <- config; compatible$p12 <- prior
  fit <- run_coloc(g, one_gene, meta, compatible, tempfile('coloc-prior-'))
  posterior <- unlist(fit[paste0('PP_H', 0:4)], use.names = FALSE)
  stopifnot(nrow(fit) == 1L, all(is.finite(posterior)),
            abs(sum(posterior) - 1) < 1e-12, fit$p12_requested == prior)
  if (prior == 0) stopifnot(fit$p12 == 0, fit$PP_H4 == 0)
  else stopifnot(fit$p12 == 1 / (nrow(one_gene) + 1))
}
cat('PASS: prior validation, output preservation, zero prior and automatic adjustment\n')
