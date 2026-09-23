# Independently validate the full real-data output using the existing R pipeline.
# This script does not invoke coloc or compute posterior probabilities.
source('R/bootstrap.R'); source('R/utils.R'); source('R/coloc_pipeline.R')
source('config.R')
base <- 'data/processed/harmonization'
read_gz <- function(name) {
  con <- gzfile(file.path(base, name), 'rt')
  on.exit(close(con))
  read.delim(con, stringsAsFactors = FALSE, check.names = FALSE, quote = '', comment.char = '')
}
g <- validate_stats(read_gz('gwas.cleaned.tsv.gz'), 'gwas')
e <- validate_stats(read_gz('eqtl.cleaned.tsv.gz'), 'eqtl')
z <- read_gz('shared.aligned.tsv.gz')
qc <- read_tsv(file.path(base, 'group_qc.tsv'))
stopifnot(all(g$build == 'GRCh38'), all(e$build == 'GRCh38'),
          all(usable_snvs(g)), all(usable_snvs(e)),
          !anyDuplicated(paste(z$locus, z$gene_id, z$tissue, z$key)),
          sum(qc$n_shared) == nrow(z))
group_key <- function(d) paste(d$locus, d$gene_id, d$tissue)
e_groups <- split(seq_len(nrow(e)), group_key(e))
z_groups <- split(seq_len(nrow(z)), group_key(z))
g_groups <- split(seq_len(nrow(g)), g$locus)
close_enough <- function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-12))
for (i in seq_len(nrow(qc))) {
  q <- qc[i, ]
  id <- group_key(q)
  a <- g[g_groups[[q$locus]] %||% integer(), ]
  b <- e[e_groups[[id]] %||% integer(), ]
  actual <- harmonize_pair(a, b)
  expected <- z[z_groups[[id]] %||% integer(), ]
  expected <- expected[order(expected$pos, expected$key), ]
  stopifnot(nrow(actual) == q$n_shared, nrow(actual) == nrow(expected),
            nrow(a) == q$n_gwas_clean, nrow(b) == q$n_eqtl_clean,
            sum(actual$allele_flipped) == q$n_flipped,
            all(actual$key == expected$key),
            all(actual$effect_allele_gwas == expected$effect_allele),
            all(actual$other_allele_gwas == expected$other_allele),
            close_enough(actual$beta_eqtl_aligned, expected$beta_eqtl_aligned),
            close_enough(actual$beta_gwas, expected$beta_gwas),
            close_enough(actual$se_eqtl, expected$se_eqtl),
            close_enough(actual$se_gwas, expected$se_gwas),
            close_enough(actual$p_eqtl, expected$p_eqtl),
            close_enough(actual$p_gwas, expected$p_gwas))
}
stopifnot(all(qc$status[qc$n_shared >= config$min_shared_snps] == 'meets_min_shared_snps'))
result <- list(status = 'passed', checked_groups = nrow(qc), checked_shared_rows = nrow(z),
               checked_flips = sum(z$allele_flipped),
               validation = 'Full output independently reproduced by R validate_stats and harmonize_pair',
               completed_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE))
jsonlite::write_json(result, file.path(base, 'r_validation.json'), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(base, 'validation_sessionInfo.txt'))
cat(jsonlite::toJSON(result, pretty = TRUE, auto_unbox = TRUE), '\n')
