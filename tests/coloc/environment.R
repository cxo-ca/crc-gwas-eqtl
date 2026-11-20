# CI-only environment; the research environment is recorded separately in logs/.
source('R/bootstrap.R')
coloc_test_library <- normalizePath('.Rlib', mustWork = FALSE)
.libPaths(c(coloc_test_library, .Library), include.site = FALSE)
if (getRversion() != '4.6.1') stop('Coloc smoke test requires R 4.6.1')
coloc_test_packages <- read.delim('tests/coloc/package_versions.tsv', stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(coloc_test_packages$Package))
check_coloc_test_environment <- function() {
  installed <- installed.packages(lib.loc = coloc_test_library)
  actual <- installed[match(coloc_test_packages$Package, installed[, 'Package']), 'Version']
  bad <- is.na(actual) | actual != coloc_test_packages$Version
  if (any(bad)) stop('Restore pinned coloc test packages first: ',
                     paste(coloc_test_packages$Package[bad], collapse = ', '))
}
