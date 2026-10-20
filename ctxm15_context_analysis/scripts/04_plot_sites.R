#!/usr/bin/env Rscript
# Circular map of chromosomal insertion sites: points on the inner circle, outer histogram = genomes per 0.1 min.
# Usage: 04_plot_sites.R results/insertion_sites.tsv results/insertion_sites.pdf [variant] [st]
a <- commandArgs(TRUE)
d <- read.delim(a[1]); d <- d[d$position_min != "ND", ]
if (length(a) >= 3) d <- d[d$gene == a[3], ]
if (length(a) >= 4) d <- d[d$st == a[4], ]
d$position_min <- as.numeric(d$position_min)
h <- table(factor(round(d$position_min, 1), levels = seq(0, 99.9, 0.1)))
ang <- function(m) pi / 2 - 2 * pi * m / 100           # 0 min at top, clockwise
pdf(a[2], 7, 7); par(mar = c(1, 1, 2, 1)); plot(NA, xlim = c(-1.6, 1.6), ylim = c(-1.6, 1.6), asp = 1, axes = FALSE, xlab = "", ylab = "",
  main = paste(c("Chromosomal blaCTX-M insertions", a[3:length(a)][-(1:2)]), collapse = " "))
t <- seq(0, 2 * pi, length.out = 400); lines(cos(t), sin(t))
for (m in seq(0, 90, 10)) { text(1.12 * cos(ang(m)), 1.12 * sin(ang(m)), m, cex = 0.7); segments(cos(ang(m)), sin(ang(m)), 1.04 * cos(ang(m)), 1.04 * sin(ang(m))) }
points(0.9 * cos(ang(d$position_min)), 0.9 * sin(ang(d$position_min)), pch = 16, col = "#00000055")
mx <- max(h, 1)
for (i in which(h > 0)) { m <- as.numeric(names(h)[i]); r <- 1.2 + 0.35 * h[i] / mx
  segments(1.2 * cos(ang(m)), 1.2 * sin(ang(m)), r * cos(ang(m)), r * sin(ang(m)), lwd = 2, col = "#b2182b") }
dev.off()
