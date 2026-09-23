`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = '\t', quote = FALSE, row.names = FALSE,
              na = 'NA', fileEncoding = 'UTF-8')
}
read_tsv <- function(path) read.delim(path, stringsAsFactors = FALSE,
  check.names = FALSE, fileEncoding = 'UTF-8', quote = '', comment.char = '')
bind_rows <- function(rows) {
  rows <- Filter(function(x) !is.null(x) && nrow(x) > 0, rows)
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}
api_url <- function(base, params) {
  paste0(base, '?', paste(unlist(Map(function(k, v) {
    paste0(k, '=', vapply(as.character(v), URLencode, character(1), reserved = TRUE))
  }, names(params), params)), collapse = '&'))
}
fetch_json <- function(url, path, offline = FALSE) {
  source_path <- paste0(path, '.source.txt')
  if (file.exists(path) && file.exists(source_path)) {
    cached_url <- readLines(source_path, n = 1, warn = FALSE)
    if (!identical(cached_url, url)) stop('Cache URL differs from this request: ', path)
  }
  if (!file.exists(path)) {
    if (offline) stop('Missing offline cache: ', path)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    partial <- paste0(path, '.part')
    ok <- FALSE
    for (attempt in 1:3) {
      ok <- tryCatch({
        download.file(url, partial, mode = 'wb', quiet = TRUE, method = 'libcurl')
        jsonlite::fromJSON(partial, simplifyVector = FALSE)
        TRUE
      }, error = function(e) FALSE)
      if (ok) break
      Sys.sleep(attempt)
    }
    if (!ok) stop('Download failed (no result substituted): ', url)
    if (!file.rename(partial, path)) stop('Cannot save cache: ', path)
    writeLines(c(url, format(Sys.time(), tz = 'UTC', usetz = TRUE)), source_path)
  }
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}
fetch_gtex <- function(endpoint, params, prefix, offline = FALSE) {
  result <- list(); page <- 0L
  repeat {
    u <- api_url(paste0('https://gtexportal.org/api/v2/', endpoint),
                 c(params, list(itemsPerPage = 1000, page = page)))
    z <- fetch_json(u, sprintf('%s_%03d.json', prefix, page), offline)
    if (is.null(z$data) || is.null(z$paging_info)) stop('Unexpected GTEx schema: ', u)
    result <- c(result, z$data)
    if (page + 1L >= z$paging_info$numberOfPages) break
    page <- page + 1L
    if (page > 10000L) stop('Unexpected pagination limit')
  }
  result
}
manifest <- function(paths, output) {
  paths <- paths[file.exists(paths)]
  write_tsv(data.frame(path = paths, bytes = file.info(paths)$size,
    md5 = unname(tools::md5sum(paths))), output)
}
html_escape <- function(x) {
  x <- gsub('&', '&amp;', as.character(x), fixed = TRUE)
  x <- gsub('<', '&lt;', x, fixed = TRUE)
  gsub('>', '&gt;', x, fixed = TRUE)
}
html_table <- function(d, n = 30) {
  if (!nrow(d)) return('<p>해당 조건에 맞는 결과가 없습니다.</p>')
  d <- head(d, n)
  paste0('<table><thead><tr>', paste0('<th>', html_escape(names(d)), '</th>', collapse = ''),
    '</tr></thead><tbody>', paste(apply(d, 1, function(r) paste0('<tr>',
      paste0('<td>', html_escape(r), '</td>', collapse = ''), '</tr>')), collapse = ''),
    '</tbody></table>')
}
