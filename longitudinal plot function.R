# fn to plot a single gene over time:
longplot = function(x, use = TRUE, ylab = NULL, col = NULL, labelhdlines = FALSE, yaxt = TRUE) {

  use = use & (!is.na(annot$donation.h) & !is.na(x))
  if (is.null(col)) {
    col = cols$covid[annot$covid]
  }
  if (yaxt) {
    plot(x[use] ~ annot$donation.h[use], col = alpha(col[use], 0.5), pch = 16, 
       xlab = "Days since onset of symptoms", ylab = ylab, cex.lab = 1.3,
       main = NULL)
  } else {
    plot(x[use] ~ annot$donation.h[use], col = alpha(col[use], 0.5), pch = 16, 
         xlab = "Days since onset of symptoms", ylab = ylab, cex.lab = 1.3,
         main = NULL,
         yaxt = "n")
  }
  
  for (id in names(table(annot$Unique.De.ID[use]))[table(annot$Unique.De.ID[use]) > 1]) {
    use2 = (annot$Unique.De.ID == id) & use
    lines(annot[use2, "donation"], x[use2], col = alpha("grey70", 0.5))
  }
  if (sum(!is.na(x[annot$covid == "HD"])) > 1) {
    normalmean = mean(x[annot$covid == "HD"])
    normalquants = quantile(x[annot$covid == "HD"], c(0.1, 0.9))
    abline(h = normalmean, col = cols$covid["HD"])
    abline(h = normalquants, col = cols$covid["HD"], lty = 2)
    if (labelhdlines) {
      text(c(250,250), normalquants, paste0(c(0.1, 0.9), " quantile\n of HD"), 
         col = cols$covid["HD"], cex = 0.7)
    }
  }
  lowessuse = ((annot$covid == "CCD") & !is.na(annot$donation)) &
    ((annot$donation < 225) & (!is.na(x)))
  lines(lowess(annot$donation.h[lowessuse], 
               x[lowessuse]), 
        col = cols$covid["CCD"], lwd = 2)
  if (is.null(col)) {
    legend("topright", pch = c(16, 16, NA), lwd = c(NA, NA, 2),
           col = cols$covid, legend = names(cols$covid))
    
  }
}
