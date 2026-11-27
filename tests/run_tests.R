source('R/bootstrap.R'); source('R/utils.R'); source('R/coloc_pipeline.R')
must_fail <- function(expr, expected) {
  err <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(err, 'error'),
            grepl(expected, conditionMessage(err), fixed = TRUE))
}
g <- data.frame(locus = 'L1', chr = 1, pos = 1:4, effect_allele = c('A','A','A','C'),
  other_allele = c('C','G','T','G'), beta = c(.2,.3,.1,.4), se = .05,
  p = c(.001,.002,.03,.01), build = 'GRCh38')
e <- g; e$gene_id <- 'TEST'; e$tissue <- 'TEST'; e$sdY <- 1
e$effect_allele[1] <- 'C'; e$other_allele[1] <- 'A'; e$beta[1] <- -.2
a <- validate_stats(g, 'gwas'); b <- validate_stats(e, 'eqtl')
h <- harmonize_pair(a, b)
stopifnot(nrow(h) == 2, sum(h$allele_flipped) == 1,
          isTRUE(all.equal(h$beta_eqtl_aligned[1], .2)))
bad <- g; bad$se[1] <- 0; must_fail(validate_stats(bad, 'gwas'), 'invalid SE')
bad <- g; bad$p[1] <- -1; must_fail(validate_stats(bad, 'gwas'), 'p must be in (0,1]')
bad <- g; bad$beta[1] <- NA
must_fail(validate_stats(bad, 'gwas'), 'missing values in required fields')
must_fail(validate_stats(rbind(g, g[1, ]), 'gwas'), 'duplicate variant within analysis group')
bad <- g; bad$chr[1] <- 2
must_fail(validate_stats(bad, 'gwas'), 'each locus must lie on a single chromosome')
stopifnot(nrow(harmonize_pair(a[0, ], b)) == 0)
# sdY describes the gene/tissue phenotype, including variants later excluded.
bad <- e; bad$sdY[1] <- 2
invalid_sdY <- list(shared = bad)
bad <- e; bad$sdY[3] <- 2
invalid_sdY$palindromic <- bad
bad <- e; bad$pos[1] <- 1000; bad$sdY[1] <- 2
invalid_sdY$unshared <- bad
other_locus <- e; other_locus$locus <- 'L2'
stopifnot(nrow(validate_stats(rbind(e, other_locus), 'eqtl')) == 2 * nrow(e))
other_locus$sdY <- 2
invalid_sdY$across_loci <- rbind(e, other_locus)
for (bad in invalid_sdY) {
  must_fail(validate_stats(bad, 'eqtl'), 'sdY must be constant within gene/tissue')
}
other_gene <- e; other_gene$gene_id <- 'OTHER_GENE'; other_gene$sdY <- 2
other_tissue <- e; other_tissue$tissue <- 'OTHER_TISSUE'; other_tissue$sdY <- 3
stopifnot(nrow(validate_stats(rbind(e, other_gene, other_tissue), 'eqtl')) == 3 * nrow(e))
source('config.R')
meta <- list(full_region_gwas = TRUE, full_region_eqtl = TRUE,
  allele_strand = 'forward', gwas_source = 'test', eqtl_source = 'test')
outdir <- tempfile('rejected-coloc-')
invalid_priors <- list(NULL, numeric(), NA_real_, NaN, Inf, -Inf, -1e-5,
                      c(1e-5, NA_real_), c(1e-5, -1e-5), '1e-5', TRUE, 1+0i, list(1e-5))
for (priors in invalid_priors) {
  bad_config <- config; bad_config$p12 <- priors
  must_fail(run_coloc(g, e, meta, bad_config, outdir), 'Invalid config$p12')
}
for (prior_name in c('p1', 'p2')) {
  for (priors in c(invalid_priors, list(c(1e-5, 1e-4)))) {
    bad_config <- config; bad_config[[prior_name]] <- priors
    must_fail(run_coloc(g, e, meta, bad_config, outdir), paste0('Invalid config$', prior_name))
  }
}
invalid_qc <- list(
  min_shared_snps = list(NULL, numeric(), NA_real_, NaN, Inf, -Inf, 0, -1, 1.5,
                         c(1, 2), '50', TRUE, 1+0i, list(1)),
  pp4_threshold = list(NULL, numeric(), NA_real_, NaN, Inf, -Inf, -.01, 1.01,
                      c(.5, .8), '0.8', TRUE, 1+0i, list(.8)))
for (field in names(invalid_qc)) {
  for (value in invalid_qc[[field]]) {
    bad_config <- config; bad_config[[field]] <- value
    must_fail(run_coloc(g, e, meta, bad_config, outdir), paste0('Invalid config$', field))
  }
}
for (bad in invalid_sdY) {
  must_fail(run_coloc(g, bad, meta, config, outdir), 'sdY must be constant within gene/tissue')
}
for (field in c('full_region_gwas', 'full_region_eqtl')) {
  bad_meta <- meta; bad_meta[[field]] <- FALSE
  must_fail(run_coloc(g, e, bad_meta, config, outdir),
            'Coloc requires unthresholded regional statistics for both traits')
}
bad_meta <- meta; bad_meta$allele_strand <- 'reverse'
must_fail(run_coloc(g, e, bad_meta, config, outdir), 'Confirm forward-strand alleles')
for (field in c('gwas_source', 'eqtl_source')) {
  bad_meta <- meta; bad_meta[[field]] <- ''
  must_fail(run_coloc(g, e, bad_meta, config, outdir), 'Record both data sources in metadata')
}
bad <- g; bad$build <- 'GRCh37'
must_fail(run_coloc(bad, e, meta, config, outdir), 'Genome build mismatch')
bad <- e; bad$build <- 'GRCh37'
must_fail(run_coloc(g, bad, meta, config, outdir), 'Genome build mismatch')
stopifnot(!dir.exists(outdir))
checks <- c('allele flips', 'ambiguous SNP removal', 'invalid input', 'duplicates',
  'single chromosome per locus', 'empty overlap',
  'sdY per gene/tissue', 'prior settings', 'QC settings', 'full-region guard', 'strand guard',
  'source metadata guard', 'build guard')
writeLines(c(paste('PASS:', paste(checks, collapse = ', ')),
  format(Sys.time(), tz = 'UTC', usetz = TRUE)), 'logs/test_results.txt')
cat(readLines('logs/test_results.txt'), sep = '\n')
