source('R/bootstrap.R')
args <- commandArgs(trailingOnly = TRUE)
packages <- c('jsonlite', if ('--coloc' %in% args) 'coloc')
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, lib = .libPaths()[1])
stopifnot(all(vapply(packages, requireNamespace, logical(1), quietly = TRUE)))
write.csv(as.data.frame(installed.packages()[, c('Package', 'Version')]),
          'logs/package_versions.csv', row.names = FALSE)
writeLines(capture.output(sessionInfo()), 'logs/setup_sessionInfo.txt')
