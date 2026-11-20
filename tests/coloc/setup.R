# Run from the project root. Install only the recorded closure, never latest CRAN.
source('tests/coloc/environment.R')
snapshot <- 'https://packagemanager.posit.co/cran/'
if (.Platform$OS.type != 'windows') snapshot <- paste0(snapshot, '__linux__/noble/')
snapshot <- paste0(snapshot, '2026-09-23')
options(timeout = 300)
installed <- installed.packages(lib.loc = coloc_test_library)
actual <- installed[match(coloc_test_packages$Package, installed[, 'Package']), 'Version']
needed <- coloc_test_packages$Package[is.na(actual) | actual != coloc_test_packages$Version]
if (length(needed)) {
  install.packages(needed, lib = coloc_test_library, repos = snapshot, dependencies = FALSE)
}
check_coloc_test_environment()
cat('PASS: pinned coloc test environment (', nrow(coloc_test_packages), ' packages)\n', sep = '')
