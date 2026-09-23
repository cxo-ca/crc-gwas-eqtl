source('R/bootstrap.R'); source('R/utils.R'); source('R/coloc_pipeline.R')
must_fail <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), 'error'))
g <- data.frame(locus = 'L1', chr = 1, pos = 1:4, effect_allele = c('A','A','A','C'),
  other_allele = c('C','G','T','G'), beta = c(.2,.3,.1,.4), se = .05,
  p = c(.001,.002,.03,.01), build = 'GRCh38')
e <- g; e$gene_id <- 'TEST'; e$tissue <- 'TEST'; e$sdY <- 1
e$effect_allele[1] <- 'C'; e$other_allele[1] <- 'A'; e$beta[1] <- -.2
a <- validate_stats(g, 'gwas'); b <- validate_stats(e, 'eqtl')
h <- harmonize_pair(a, b)
stopifnot(nrow(h) == 2, sum(h$allele_flipped) == 1,
          isTRUE(all.equal(h$beta_eqtl_aligned[1], .2)))
bad <- g; bad$se[1] <- 0; must_fail(validate_stats(bad, 'gwas'))
bad <- g; bad$p[1] <- -1; must_fail(validate_stats(bad, 'gwas'))
bad <- g; bad$beta[1] <- NA; must_fail(validate_stats(bad, 'gwas'))
must_fail(validate_stats(rbind(g, g[1, ]), 'gwas'))
bad <- g; bad$chr[1] <- 2; must_fail(validate_stats(bad, 'gwas'))
stopifnot(nrow(harmonize_pair(a[0, ], b)) == 0)
source('config.R')
if (requireNamespace('coloc', quietly = TRUE)) {
  meta <- list(full_region_gwas = FALSE, full_region_eqtl = TRUE,
    allele_strand = 'forward', gwas_source = 'test', eqtl_source = 'test')
  must_fail(run_coloc(g, e, meta, config, tempdir()))
  meta$full_region_gwas <- TRUE
  bad <- g; bad$build <- 'GRCh37'
  must_fail(run_coloc(bad, e, meta, config, tempdir()))
}
checks <- c('allele flips', 'ambiguous SNP removal', 'invalid input', 'duplicates',
  'single chromosome per locus', 'empty overlap',
  if (requireNamespace('coloc', quietly = TRUE)) c('full-region guard', 'build guard'))
writeLines(c(paste('PASS:', paste(checks, collapse = ', ')),
  format(Sys.time(), tz = 'UTC', usetz = TRUE)), 'logs/test_results.txt')
cat(readLines('logs/test_results.txt'), sep = '\n')
