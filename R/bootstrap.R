# Run scripts from the project root. Keep packages inside this project.
if (.Platform$OS.type == 'windows') {
  try(Sys.setlocale('LC_ALL', 'English_United States.utf8'), silent = TRUE)
}
dir.create('.Rlib', showWarnings = FALSE)
.libPaths(c(normalizePath('.Rlib'), .libPaths()))
options(repos = c(CRAN = 'https://cloud.r-project.org'), timeout = 180)
for (d in c('data/raw', 'data/processed', 'data/input', 'results/public',
            'results/demo', 'results/coloc', 'reports', 'logs')) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
