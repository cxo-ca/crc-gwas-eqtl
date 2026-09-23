validate_stats <- function(d, kind) {
  required <- c('locus', 'chr', 'pos', 'effect_allele', 'other_allele',
                'beta', 'se', 'p', 'build')
  if (kind == 'eqtl') required <- c(required, 'gene_id', 'tissue', 'sdY')
  missing <- setdiff(required, names(d))
  if (length(missing)) stop(kind, ': missing columns: ', paste(missing, collapse = ', '))
  if (!nrow(d)) stop(kind, ': empty input')
  if (anyNA(d[required])) stop(kind, ': missing values in required fields')
  for (nm in c('pos', 'beta', 'se', 'p', if (kind == 'eqtl') 'sdY')) {
    if (!is.numeric(d[[nm]]) || any(!is.finite(d[[nm]]))) stop(kind, ': invalid ', nm)
  }
  if (any(d$pos < 1 | d$pos != floor(d$pos))) stop(kind, ': invalid position')
  if (any(d$se <= 0) || any(!is.finite(d$se^2)) || any(d$se^2 == 0)) stop(kind, ': invalid SE')
  if (any(d$p <= 0 | d$p > 1)) stop(kind, ': p must be in (0,1]; preserve very small p-values')
  if (kind == 'eqtl' && any(d$sdY <= 0)) stop('eqtl: sdY must be positive')
  d$chr <- sub('^chr', '', as.character(d$chr), ignore.case = TRUE)
  if (any(!d$chr %in% c(as.character(1:22), 'X', 'Y'))) stop(kind, ': invalid chromosome')
  d$effect_allele <- toupper(d$effect_allele)
  d$other_allele <- toupper(d$other_allele)
  if (any(d$effect_allele == d$other_allele)) stop(kind, ': identical alleles')
  if (any(!nzchar(d$locus))) stop(kind, ': empty locus')
  chromosomes <- split(d$chr, d$locus)
  if (any(vapply(chromosomes, function(x) length(unique(x)) != 1, logical(1)))) {
    stop(kind, ': each locus must lie on a single chromosome')
  }
  d$key <- paste(d$chr, d$pos, pmin(d$effect_allele, d$other_allele),
                 pmax(d$effect_allele, d$other_allele), sep = ':')
  identity <- if (kind == 'eqtl') paste(d$locus, d$gene_id, d$tissue, d$key) else paste(d$locus, d$key)
  if (anyDuplicated(identity)) stop(kind, ': duplicate variant within analysis group')
  d
}
usable_snvs <- function(d) {
  snv <- d$effect_allele %in% c('A', 'C', 'G', 'T') &
         d$other_allele %in% c('A', 'C', 'G', 'T')
  pair <- paste0(pmin(d$effect_allele, d$other_allele), pmax(d$effect_allele, d$other_allele))
  snv & !pair %in% c('AT', 'CG')
}
harmonize_pair <- function(g, e) {
  # Both inputs must use forward-strand alleles on the same reference build.
  # Ambiguous SNVs and indels are excluded; strand/build conversion is external.
  g <- g[usable_snvs(g), ]; e <- e[usable_snvs(e), ]
  z <- merge(g, e, by = c('locus', 'key', 'chr', 'pos', 'build'), suffixes = c('_gwas', '_eqtl'))
  flip <- z$effect_allele_gwas == z$other_allele_eqtl &
          z$other_allele_gwas == z$effect_allele_eqtl
  same <- z$effect_allele_gwas == z$effect_allele_eqtl &
          z$other_allele_gwas == z$other_allele_eqtl
  if (any(!same & !flip)) stop('Unresolved allele mismatch')
  z$beta_eqtl_aligned <- ifelse(flip, -z$beta_eqtl, z$beta_eqtl)
  z$allele_flipped <- flip
  z[order(z$pos, z$key), ]
}
run_coloc <- function(gwas, eqtl, metadata, config, outdir) {
  if (!requireNamespace('coloc', quietly = TRUE)) stop('Run scripts/setup.R --coloc first')
  if (!isTRUE(metadata$full_region_gwas) || !isTRUE(metadata$full_region_eqtl)) {
    stop('Coloc requires unthresholded regional statistics for both traits')
  }
  if (!identical(metadata$allele_strand, 'forward')) stop('Confirm forward-strand alleles')
  if (!nzchar(metadata$gwas_source %||% '') || !nzchar(metadata$eqtl_source %||% '')) {
    stop('Record both data sources in metadata')
  }
  g <- validate_stats(gwas, 'gwas'); e <- validate_stats(eqtl, 'eqtl')
  if (any(g$build != config$genome_build) || any(e$build != config$genome_build)) {
    stop('Genome build mismatch: explicit liftover and allele checking required')
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  groups <- unique(e[c('locus', 'gene_id', 'tissue')])
  output <- list(); qc <- list(); all_harmonized <- list(); warning_rows <- list()
  g_index <- split(seq_len(nrow(g)), g$locus)
  e_index <- split(seq_len(nrow(e)), paste(e$locus, e$gene_id, e$tissue))
  for (i in seq_len(nrow(groups))) {
    group <- groups[i, ]
    a <- g[g_index[[group$locus]] %||% integer(), ]
    b <- e[e_index[[paste(group$locus, group$gene_id, group$tissue)]], ]
    z <- harmonize_pair(a, b)
    status <- if (nrow(z) < config$min_shared_snps) 'insufficient_shared_snps' else 'analyzed'
    qc[[i]] <- cbind(group, data.frame(n_gwas = nrow(a), n_eqtl = nrow(b),
      n_gwas_excluded_alleles = sum(!usable_snvs(a)), n_eqtl_excluded_alleles = sum(!usable_snvs(b)),
      n_shared = nrow(z), n_flipped = sum(z$allele_flipped),
      fraction_gwas_shared = if (nrow(a)) nrow(z)/nrow(a) else 0,
      fraction_eqtl_shared = nrow(z)/nrow(b), status = status))
    if (status != 'analyzed') next
    if (length(unique(z$sdY)) != 1) stop('sdY must be constant within gene/tissue')
    d1 <- list(snp = z$key, position = z$pos, beta = z$beta_gwas,
               varbeta = z$se_gwas^2, type = 'cc')
    d2 <- list(snp = z$key, position = z$pos, beta = z$beta_eqtl_aligned,
               varbeta = z$se_eqtl^2, type = 'quant', sdY = z$sdY[1])
    for (prior in config$p12) {
      fit_warnings <- character()
      invisible(capture.output(fit <- withCallingHandlers(
        coloc::coloc.abf(d1, d2, p1 = config$p1, p2 = config$p2, p12 = prior),
        warning = function(w) {
          fit_warnings <<- c(fit_warnings, conditionMessage(w))
          invokeRestart('muffleWarning')
        })))
      if (length(fit_warnings)) warning_rows[[length(warning_rows) + 1L]] <- cbind(
        group, data.frame(p12_requested = prior, message = gsub('[\r\n\t]+', ' ', fit_warnings)))
      sm <- fit$summary
      output[[length(output) + 1L]] <- cbind(group, data.frame(
        nsnps = nrow(z), p1 = unname(fit$priors['p1']), p2 = unname(fit$priors['p2']),
        p12 = unname(fit$priors['p12']), p1_requested = config$p1,
        p2_requested = config$p2, p12_requested = prior,
        PP_H0 = unname(sm['PP.H0.abf']), PP_H1 = unname(sm['PP.H1.abf']),
        PP_H2 = unname(sm['PP.H2.abf']), PP_H3 = unname(sm['PP.H3.abf']),
        PP_H4 = unname(sm['PP.H4.abf']),
        passes_pp4 = unname(sm['PP.H4.abf']) >= config$pp4_threshold))
    }
    all_harmonized[[i]] <- z
    png(file.path(outdir, sprintf('locus_%03d.png', i)), width = 1500, height = 850, res = 150)
    par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
    plot(a$pos/1e6, -log10(a$p), pch = 16, col = '#CCCCCC',
      xlab = 'Position (Mb)', ylab = '-log10(GWAS p)', main = as.character(group$locus))
    points(z$pos/1e6, -log10(z$p_gwas), pch = 16, col = '#087E8B')
    legend('topright', c('Shared SNPs', 'Other QC-passing SNPs'), pch = 16,
      col = c('#087E8B', '#CCCCCC'), cex = .7, bty = 'n')
    plot(b$pos/1e6, -log10(b$p), pch = 16, col = '#CCCCCC',
      xlab = 'Position (Mb)', ylab = '-log10(eQTL p)', main = paste(group$gene_id, group$tissue))
    points(z$pos/1e6, -log10(z$p_eqtl), pch = 16, col = '#B45F35')
    dev.off()
    cat(sprintf('[%d/%d] %s %s %s: %d SNPs\n', i, nrow(groups), group$locus,
      group$gene_id, group$tissue, nrow(z)))
  }
  results <- bind_rows(output)
  write_tsv(bind_rows(qc), file.path(outdir, 'qc.tsv'))
  warnings <- if (length(warning_rows)) bind_rows(warning_rows) else data.frame(
    locus = character(), gene_id = character(), tissue = character(),
    p12_requested = numeric(), message = character())
  write_tsv(warnings, file.path(outdir, 'warnings.tsv'))
  if (!nrow(results)) stop('No eligible coloc groups. Inspect QC; no posterior results written.')
  write_tsv(results, file.path(outdir, 'coloc_sensitivity.tsv'))
  write_tsv(bind_rows(all_harmonized), file.path(outdir, 'harmonized.tsv'))
  saveRDS(list(config = config, metadata = metadata), file.path(outdir, 'analysis_settings.rds'))
  writeLines(capture.output(sessionInfo()), file.path(outdir, 'sessionInfo.txt'))
  results
}
