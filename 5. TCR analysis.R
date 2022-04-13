
rm(list = ls())
library(lme4)
library(lmerTest)
library(ComplexHeatmap)
library(viridis)
library(scales)
library(lmtest)
library(lme4)
library(lmerTest)

source("longitudinal plot function.R")

#### load data and DE results --------------------------

annot = read.csv("processed data/sample annotations.csv", row.names = 1)
load("processed data/cols.RData")
load("processed data/perturbation clustering results.RData") # loads "hmdf"

norm = read.csv("processed data/TCR normalized data.csv", row.names = 1)
norm = norm[rownames(annot), ]

cols$pclass = cols$covid
cols$pclass["p1"] = "red"
cols$pclass["p2"] = "darkviolet"


svg("results/fig 3a TCR longitudinal plot.svg")
longplot(annot$tcr, ylab = "TCR score", col = cols$pclass[hmdf$Cluster])
dev.off()


# break samples into time windows:
breaks = c(0, 90, 120, 150, 360)
timegrps = list()
for (i in 1:length(breaks[-1])) {
  timegrps[[i]] = setdiff(rownames(annot)[(annot$donation >= breaks[i]) & (annot$donation < breaks[i+1])], NA)
}
names(timegrps) = paste0(breaks[-length(breaks)], "-", breaks[-1])
sapply(timegrps, length)

timegrp = rep(NA, nrow(annot))
names(timegrp) = rownames(annot)
for (name in names(timegrps)) {
  timegrp[timegrps[[name]]] = name
}

# run ccd vs hd DE within each 
stratests = stratps = stratfdrs = stratci = list()
for (name in names(timegrps)) {
  print(name)
  stratests[[name]] = stratps[[name]] = c()
  
  use = is.element(rownames(annot), timegrps[[name]]) | (annot$covid == 'HD')
  temp = annot[use,]
  
  if (max(table(annot$Unique.De.ID[use])) > 1) {
    mod = lmer(data = temp, tcr ~ covid + Age + Sex + RACE2 + (1|Unique.De.ID))
    stratests[[name]] = -summary(mod)$coefficients["covidHD", "Estimate"]
    stratci[[name]] = -summary(mod)$coefficients["covidHD", "Estimate"] + c(-1,1) * qt(0.975, summary(mod)$coefficients["covidHD", "df"]) * summary(mod)$coefficients["covidHD", "Std. Error"]
    stratps[[name]] = summary(mod)$coefficients["covidHD", "Pr(>|t|)"]
  } else {
    mod = lm(data = temp, tcr ~ covid + Age + Sex + RACE2)
    stratests[[name]] = -summary(mod)$coefficients["covidHD", "Estimate"]
    stratci[[name]] = -summary(mod)$coefficients["covidHD", "Estimate"] + c(-1,1) * qt(0.975, df = mod$df.residual) * summary(mod)$coefficients["covidHD", "Std. Error"]   
    stratps[[name]] = summary(mod)$coefficients["covidHD", "Pr(>|t|)"]
  }
  stratfdrs[[name]] = p.adjust(stratps[[name]], method = "BH")
}

svg("results/fig 3b tcr forestplot.svg", height = 2.5)
par(mar = c(5,8,1,1))
plot(unlist(stratests), 1:4, xlim = c(-0.2, 0.2),
     xlab = "Change from mean of HD",
     yaxt = "n", ylab = "", pch = 16)
axis(2, at = 1:4, labels = paste0(names(stratests), " days"), las = 2)
for (i in 1:length(stratests)) {
  lines(stratci[[i]], rep(i, 2))
}
abline(v = 0, col = "grey60")
dev.off()


#### TCR score vs. perturbation class:
fac = factor(hmdf$Cluster, levels = c("HD", "CCD", "p1", "p2"))
svg("results/fig 3c boxplot of TCR vs perturbation cluster.svg")
par(mar = c(4,5,2,1))
bp = boxplot(annot$tcr ~ fac, outline = F, ylim = range(annot$tcr, na.rm = T), ylab = "TCR score", cex.lab = 1.5, cex.axis = 1.3, xlab = "")
points(jitter(as.numeric(fac)), annot$tcr, col = cols$pclass[hmdf$Cluster], pch = 16)
mod = summary(lm(annot$tcr ~ hmdf$Cluster))$coef
text(3,3.025,paste0("p = ", signif(mod["hmdf$Clusterp1", 4], 3)))
text(4,3.025,paste0("p = ", signif(mod["hmdf$Clusterp2", 4], 3)))
dev.off()



### what TCR genes differ in P1 vs. other CCD? -----------------------------------

ests = ps = c()
use = !is.na(annot$donation)
temp = cbind(annot[use,], norm[use, ])
temp$cluster = hmdf$Cluster[use]

for (gene in colnames(norm)) {
  mod = summary(lmer(data = temp, get(gene) ~ cluster + donation + (1|Unique.De.ID)))$coef
  ests[gene] = mod["clusterp1", "Estimate"]
  ps[gene] = mod["clusterp1", "Pr(>|t|)"]
}

plot(ests, -log10(ps), col = 0)
thresh = max(ps[p.adjust(ps, "BH") < 0.05])
abline(h = -log10(thresh))
text(ests, -log10(ps), names(ests), cex = 0.5)

out = data.frame("log2 fold-change P1 / CCD" = ests,
                 "p-value" = ps,
                 "BH FDR" = p.adjust(ps, "BH"))
write.csv(out, file = "results/TCR gene vs. cluster DE results.csv")

