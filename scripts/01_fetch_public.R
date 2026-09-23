source('R/bootstrap.R'); source('R/utils.R'); source('config.R')
stopifnot(requireNamespace('jsonlite', quietly = TRUE))
offline <- '--offline' %in% commandArgs(trailingOnly = TRUE)
base <- file.path('data/raw', config$study, config$gtex_dataset)
dir.create(base, recursive = TRUE, showWarnings = FALSE)
study_url <- paste0('https://www.ebi.ac.uk/gwas/rest/api/studies/', config$study)
study <- fetch_json(study_url, file.path(base, 'study.json'), offline)
associations <- list(); page <- 0L
u <- paste0(study_url, '/associations?size=500')
repeat {
  z <- fetch_json(u, file.path(base, sprintf('associations_%03d.json', page)), offline)
  a <- z$`_embedded`$associations
  if (is.null(a)) stop('GWAS Catalog association schema missing')
  associations <- c(associations, a)
  next_url <- z$`_links`$`next`$href
  if (is.null(next_url)) break
  u <- next_url; page <- page + 1L
}
rows <- lapply(associations, function(a) {
  # Exclude haplotype / interaction associations from single-variant matching.
  if (isTRUE(a$multiSnpHaplotype) || isTRUE(a$snpInteraction)) return(NULL)
  alleles <- unlist(lapply(a$loci, function(l) vapply(l$strongestRiskAlleles,
    function(r) r$riskAlleleName %||% '', character(1))))
  rs <- unique(sub('-.*$', '', alleles))
  rs <- rs[grepl('^rs[0-9]+$', rs)]
  if (!length(rs)) return(NULL)
  data.frame(rsid = rs, gwas_p = a$pvalue %||% NA_real_,
    reported_allele = paste(alleles, collapse = ';'),
    association_url = a$`_links`$self$href, study = config$study)
})
gwas_all <- bind_rows(rows)
if (!nrow(gwas_all)) stop('No valid Catalog associations')
write_tsv(gwas_all, 'data/processed/gwas_catalog_associations.tsv')
gwas <- gwas_all[is.finite(gwas_all$gwas_p) & gwas_all$gwas_p > 0 &
                 gwas_all$gwas_p <= config$gwas_p, ]
gwas <- gwas[order(gwas$gwas_p, gwas$rsid), ]
gwas <- gwas[!duplicated(gwas$rsid), ]
if (!nrow(gwas)) stop('No genome-wide significant SNPs')
write_tsv(gwas, 'data/processed/gwas_lead_variants.tsv')
mapping <- list()
for (i in seq_len(nrow(gwas))) {
  rs <- gwas$rsid[i]
  message(sprintf('Variant mapping %d/%d: %s', i, nrow(gwas), rs))
  z <- fetch_gtex('dataset/variant', list(snpId = rs, datasetId = config$gtex_dataset),
                  file.path(base, paste0('variant_', rs)), offline)
  mapping[[i]] <- bind_rows(lapply(z, function(v) data.frame(
    rsid = rs, returned_rsid = v$snpId %||% NA_character_,
    variant_id = v$variantId, chr = v$chromosome, pos = v$pos,
    ref = v$ref, alt = v$alt, build = 'GRCh38')))
}
mapping <- unique(bind_rows(mapping))
if (!nrow(mapping)) stop('No GWAS variants mapped in GTEx')
write_tsv(mapping, 'data/processed/variant_mapping.tsv')
write_tsv(data.frame(rsid = setdiff(gwas$rsid, mapping$rsid)),
          'data/processed/unmapped_variants.tsv')
ids <- unique(mapping$variant_id)
chunks <- split(ids, ceiling(seq_along(ids) / 30))
eqtl <- list()
for (tissue in config$tissues) for (i in seq_along(chunks)) {
  message('eQTL query: ', tissue, ' batch ', i, '/', length(chunks))
  z <- fetch_gtex('association/singleTissueEqtl', list(variantId = chunks[[i]],
    tissueSiteDetailId = tissue, datasetId = config$gtex_dataset),
    file.path(base, paste0('eqtl_', tissue, '_batch', i)), offline)
  eqtl[[length(eqtl) + 1L]] <- bind_rows(lapply(z, function(e) data.frame(
    variant_id = e$variantId, gene_id = e$gencodeId, gene_symbol = e$geneSymbol,
    tissue = e$tissueSiteDetailId, eqtl_p = e$pValue, nes = e$nes,
    dataset = e$datasetId)))
}
eqtl <- unique(bind_rows(eqtl))
if (!nrow(eqtl)) eqtl <- data.frame(variant_id = character(), gene_id = character(),
  gene_symbol = character(), tissue = character(), eqtl_p = numeric(), nes = numeric(),
  dataset = character())
if (nrow(eqtl) && (any(!eqtl$variant_id %in% ids) ||
    any(!eqtl$tissue %in% config$tissues) || any(eqtl$dataset != config$gtex_dataset))) {
  stop('API returned rows outside requested filters')
}
write_tsv(eqtl, 'data/processed/gtex_significant_eqtl.tsv')
saveRDS(config, 'data/processed/public_collection_config.rds')
manifest(list.files(base, recursive = TRUE, full.names = TRUE, pattern = '\\.json$'),
         'data/processed/download_manifest.tsv')
writeLines(capture.output(sessionInfo()), 'logs/public_fetch_sessionInfo.txt')
message('Public data collection complete.')
