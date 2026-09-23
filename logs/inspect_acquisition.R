source('R/bootstrap.R')
cat('data.table:', requireNamespace('data.table', quietly = TRUE), '\n')
s <- jsonlite::fromJSON('data/raw/gtex_openapi.json', simplifyVector = FALSE)
for (nm in grep('eqtl|gene|association', names(s$paths), value = TRUE, ignore.case = TRUE)) {
  z <- s$paths[[nm]]
  cat(nm, '\n')
  for (method in names(z)) {
    cat(method, z[[method]]$description, '\n')
    for (p in z[[method]]$parameters) cat(' ', p$name, ':', p$description, '\n')
  }
}
