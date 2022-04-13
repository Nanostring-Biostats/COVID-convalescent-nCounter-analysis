rm(list = ls())
library(scales)
library(readxl)
library(lmtest)
library(pheatmap)

load("processed data/cols.RData")
# load data:
d = read_xlsx("data/12012021 ELISA assay - PD - tidy.xlsx")
# annotation:
annot = as.data.frame(d[,1:2])
rownames(annot) = d$`Sample ID`
# expression:
mat = as.matrix(d[,-(1:3)])
for (i in 1:ncol(mat)) {
  colnames(mat)[i] = substr(colnames(mat)[i], 1, gregexpr(" ", colnames(mat)[i])[[1]][1] - 1)
}
# parse OOR < entries:
min.non0 = c()
for (i in 1:ncol(mat)) {
  if (any(mat[, i] == "OOR <")) {
    min.non0[i] = min(as.numeric(mat[,i]), na.rm = T)
  }
}
mat = apply(mat, 2, as.numeric)
rownames(mat) = d$`Sample ID`


# align to P1/P2 annotations:
load("processed data/perturbation clustering results.RData") # load ("hmdf")
annot$cluster = hmdf$Cluster[match(rownames(annot), gsub("\\.", " ", gsub("X", "", rownames(hmdf))))]
annot$cluster[annot$Condition == 'Healthy'] = "HD"

annot$col = c("red", "darkviolet", "dodgerblue3")[(annot$cluster == "p1") + 
                                                    2 * (annot$cluster == "p2") +
                                                    3 * (annot$cluster == "HD")]

# models for each cytokine:
ps = c()
for (name in colnames(mat)) {
  if (!all(is.na(mat[, name]))) {
    
    y = mat[, name]
    lod = min(y, na.rm = T) / 2
    y[is.na(y)] = lod
    
    #mod = lm(y ~ annot$cluster)
    #ps = rbind(ps, summary(mod)$coef[c("annot$clusterp1", "annot$clusterp2"), "Pr(>|t|)"])
    t1 = t.test(y[annot$cluster == "HD"], y[annot$cluster == "p1"])
    t2 = t.test(y[annot$cluster == "HD"], y[annot$cluster == "p2"])
    ps = rbind(ps, c(t1$p.value, t2$p.value))
    rownames(ps)[nrow(ps)] = name
  }
}

fdrs = ps * NA
for (i in 1:2) {
  fdrs[,i] = p.adjust(ps[,i], "BH")
}


out = cbind(ps, fdrs)
colnames(out) = c("p1 vs. HD p-value", "p2 vs. HD p-value", "p1 vs. HD FDR", "p2 vs. HD FDR")
write.csv(out, file = "results/cytokine p-values.csv")


for (name in rownames(ps)) {
  if (!all(is.na(mat[, name]))) {
    svg(paste0("results/cytokine plots/", name, ".svg"), width = 4, height = 4)
    # get data including below-LOD amounts:
    y = mat[, name]
    lod = min(y, na.rm = T) / 2
    thresh = mean(c(1, 0.5) * min(y, na.rm = T))
    y[is.na(y)] = lod
    
    boxplot(y ~ annot$cluster, outline = F, ylim = c(min(y), max(y) + 0.15 * diff(range(y))), col = 0, ylab = name, xlab = "")
    points(jitter(as.numeric(as.factor(annot$cluster))), y,
           pch = 16, cex = 1.5, col = alpha(annot$col, 0.5))
    abline(h = thresh, lty = 2)
    
    text(2,  max(y) + 0.11 * diff(range(y)), paste0("p = ", signif(ps[name, 1],2)))
    text(3,  max(y) + 0.11 * diff(range(y)), paste0("p = ", signif(ps[name, 2],2)))
    dev.off()
  }
}


