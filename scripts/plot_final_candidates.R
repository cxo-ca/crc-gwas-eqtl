source('R/bootstrap.R'); source('R/utils.R')
d <- read_tsv('results/final/candidate_groups.tsv')
d <- d[d$evidence_tier == 'A_conditional_prior_stable', ]
d <- d[order(d$PP_H4_min, decreasing = TRUE), ]
stopifnot(nrow(d) == 11, all(d$PP_H4_min >= .8))
draw <- function() {
  par(mar = c(5, 14, 4, 2))
  y <- rev(seq_len(nrow(d)))
  plot(NA, xlim = c(.75, 1.005), ylim = c(.4, nrow(d) + .7), yaxt = 'n',
       xlab = 'PP.H4: range over three requested p12 values (not a confidence interval)',
       ylab = '', main = 'Prior sensitivity of 11 conditionally supported gene-tissue pairs')
  abline(v = .8, lty = 2, col = '#B45F35')
  abline(h = y, col = '#eeeeee')
  segments(d$PP_H4_min, y, d$PP_H4_max, y, col = '#087E8B', lwd = 3)
  points(d$PP_H4_baseline, y, pch = 19, col = '#243342', cex = 1.1)
  axis(2, at = y, labels = paste(d$gene_symbol, sub('Colon_', '', d$tissue), sep = ' / '), las = 1, tick = FALSE)
  legend('topleft', c('Range across p12 = 1e-6, 1e-5, 1e-4', 'Baseline requested p12 = 1e-5', 'Threshold 0.8'),
         col = c('#087E8B', '#243342', '#B45F35'), lty = c(1, NA, 2),
         pch = c(NA, 19, NA), bty = 'n', cex = .8)
}
png('reports/final_candidate_sensitivity.png', width = 1900, height = 1100, res = 160)
draw(); dev.off()
pdf('reports/final_candidate_sensitivity.pdf', width = 12, height = 7)
draw(); dev.off()
cat('Saved PNG and PDF sensitivity figures.\n')
