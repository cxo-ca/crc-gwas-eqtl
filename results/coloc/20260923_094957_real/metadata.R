list(full_region_gwas = TRUE, full_region_eqtl = TRUE, allele_strand = "forward", 
    full_region_definition = "All available regional rows after documented QC; no significance threshold; does not imply all variants were measured", 
    gwas_source = "GCST90129505 GRCh38 harmonised summary statistics; see docs/DATA_ACQUISITION.md", 
    eqtl_source = "GTEx v8 original imported Colon_Sigmoid and Colon_Transverse statistics from eQTL Catalogue", 
    gwas_ancestry = "European and East Asian meta-analysis", 
    eqtl_ancestry = "GTEx v8 tissue cohorts; not ancestry-matched to the GWAS meta-analysis", 
    region_definition = "Exploratory merged lead +/-500kb GRCh38 windows in data/processed/acquisition/regions.tsv; not LD-defined", 
    expression_scale_and_sdY_source = "Official GTEx v8 normalized expression matrices; measured sample SD per versioned gene/tissue", 
    limitations = c("Single causal variant per trait per region assumption; no SuSiE or conditional analysis", 
    "Palindromic SNPs and indels excluded; four download lead SNPs removed", 
    "SHROOM2 and GPR143 X-region GWAS absent; three gene/tissue groups skipped", 
    "Shared coverage and ancestry/LD differences limit interpretation; not proof of causal mediation"
    ))
