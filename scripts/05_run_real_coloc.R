# Real-data exploratory coloc on the documented, unthresholded distance windows.
source('R/bootstrap.R'); source('R/utils.R'); source('R/coloc_pipeline.R'); source('config.R')
stopifnot(requireNamespace('digest', quietly = TRUE), requireNamespace('jsonlite', quietly = TRUE))
base <- 'data/processed/harmonization'
hs <- jsonlite::fromJSON(file.path(base, 'harmonization_summary.json'))
rv <- jsonlite::fromJSON(file.path(base, 'r_validation.json'))
stopifnot(hs$status == 'harmonized_not_coloc_run', rv$status == 'passed',
          rv$checked_shared_rows == hs$counts$shared_rows)
hash_records <- rbind(hs$inputs[c('path', 'sha256')], hs$files[c('path', 'sha256')])
for (i in seq_len(nrow(hash_records))) {
  stopifnot(digest::digest(file = hash_records$path[i], algo = 'sha256') == hash_records$sha256[i])
}
metadata <- list(
  full_region_gwas = TRUE, full_region_eqtl = TRUE, allele_strand = 'forward',
  full_region_definition = 'All available regional rows after documented QC; no significance threshold; does not imply all variants were measured',
  gwas_source = 'GCST90129505 GRCh38 harmonised summary statistics; see docs/DATA_ACQUISITION.md',
  eqtl_source = 'GTEx v8 original imported Colon_Sigmoid and Colon_Transverse statistics from eQTL Catalogue',
  gwas_ancestry = 'European and East Asian meta-analysis',
  eqtl_ancestry = 'GTEx v8 tissue cohorts; not ancestry-matched to the GWAS meta-analysis',
  region_definition = 'Exploratory merged lead +/-500kb GRCh38 windows in data/processed/acquisition/regions.tsv; not LD-defined',
  expression_scale_and_sdY_source = 'Official GTEx v8 normalized expression matrices; measured sample SD per versioned gene/tissue',
  limitations = c('Single causal variant per trait per region assumption; no SuSiE or conditional analysis',
    'Palindromic SNPs and indels excluded; four download lead SNPs removed',
    'SHROOM2 and GPR143 X-region GWAS absent; three gene/tissue groups skipped',
    'Shared coverage and ancestry/LD differences limit interpretation; not proof of causal mediation'))
read_clean <- function(kind) {
  con <- gzfile(file.path(base, paste0(kind, '.cleaned.tsv.gz')), 'rt')
  on.exit(close(con))
  d <- read.delim(con, stringsAsFactors = FALSE, check.names = FALSE, quote = '', comment.char = '')
  cols <- c('locus', 'chr', 'pos', 'effect_allele', 'other_allele', 'beta', 'se', 'p', 'build',
            if (kind == 'eqtl') c('gene_id', 'tissue', 'sdY'))
  d[cols]
}
outdir <- file.path('results/coloc', paste0(format(Sys.time(), '%Y%m%d_%H%M%S'), '_real'))
if (dir.exists(outdir)) stop('Output directory already exists')
dir.create(outdir, recursive = TRUE)
dput(metadata, file = file.path(outdir, 'metadata.R'))
write_tsv(read_tsv(file.path(base, 'group_qc.tsv')), file.path(outdir, 'preanalysis_qc.tsv'))
write_tsv(read_tsv(file.path(base, 'lead_variant_qc.tsv')), file.path(outdir, 'lead_variant_qc.tsv'))
cat('Output:', outdir, '\nReading real regional inputs...\n')
g <- read_clean('gwas'); e <- read_clean('eqtl')
result <- run_coloc(g, e, metadata, config, outdir)
qc <- read_tsv(file.path(outdir, 'qc.tsv'))
stopifnot(nrow(result) == hs$groups_meeting_min_shared_snps * length(config$p12),
          sum(qc$n_shared) == hs$counts$shared_rows,
          all(abs(rowSums(result[paste0('PP_H', 0:4)]) - 1) < 1e-10))
inputs <- unique(c(hash_records$path, file.path(base, c('harmonization_summary.json', 'r_validation.json')),
                   'scripts/05_run_real_coloc.R', 'R/coloc_pipeline.R', 'config.R'))
write_tsv(data.frame(path = inputs,
  sha256 = vapply(inputs, function(p) digest::digest(file = p, algo = 'sha256'), character(1))),
  file.path(outdir, 'input_sha256.tsv'))
jsonlite::write_json(list(status = 'coloc_completed', real_data = TRUE,
  coloc_version = as.character(packageVersion('coloc')), groups_analyzed = sum(qc$status == 'analyzed'),
  groups_skipped = sum(qc$status != 'analyzed'), sensitivity_rows = nrow(result),
  completed_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE)),
  file.path(outdir, 'run_summary.json'), pretty = TRUE, auto_unbox = TRUE)
writeLines(outdir, 'results/coloc/latest_real_run.txt')
cat('Completed:', outdir, '\n')
