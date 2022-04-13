rm(list = ls())
library(umap)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(mclust)
library(lme4)
library(lmerTest)
library(vioplot)
source("longitudinal plot function.R")
source("genesetplot.R")

#### load data ------------------------------
norm = as.matrix(read.csv("processed data/log2 scale normalized data.csv", row.names = 1))
trvs = as.matrix(read.csv("processed data/TCR normalized data.csv", row.names = 1))
gannot = read.csv("processed data/gene annotations.csv", row.names = 1)
annot = read.csv("processed data/sample annotations.csv", row.names = 1)
annot = annot[rownames(norm), ]
trvs = trvs[match(rownames(norm), rownames(trvs)), ]
load("processed data/cols.RData")
allde = read.csv("results/DE results from host response panel.csv", row.names = 1, check.names = F)
allde = allde[colnames(norm), ]

#### score distance from normals: --------------------------------

# z-transform genes based on normals:
hdmeans = apply(norm[annot$covid == "HD", ], 2, mean)
hdsds = apply(norm[annot$covid == "HD", ], 2, sd)
normz = sweep(sweep(norm, 2, hdmeans, "-"), 2, hdsds, "/")
annot$perturbation = sqrt(rowMeans(normz^2))

# pick a threshold for highly-perturbed samples:
mc = Mclust(annot$perturbation, G = 2, modelNames = "E")

thresh.perturb = min(annot$perturbation[mc$classification == 2])
annot$perturbed = annot$perturbation >= thresh.perturb



if (FALSE) {
  svg("results/perturbation longplot and vioplot.svg")
  layout(mat = matrix(1:2, nrow = 1), widths = c(3.5,7))
  par(mar = c(5,5,2,0))
  vioplot(annot$perturbation ~ annot$covid, col = alpha(rev(cols$covid), 0.5), 
          ylab = "Perturbation score", xlab = "", cex.lab = 1.5)
  par(mar = c(5,0,2,1))
  longplot(x = annot$perturbation, ylab = "Perturbation score", yaxt = F) 
  rect(40,thresh.perturb,160,3.05, col = NA, border = alpha("red", 0.5))
  dev.off()
  
}



#### for all perturbed samples, ID genes most responsible: -------------------------
perturbgenes = list()
for (id in rownames(annot)[annot$perturbed]) {
  perturbgenes[[id]] = colnames(normz)[abs(normz[id, ]) > 5]
}
length(unique(unlist(perturbgenes)))

# classify by cutting tree:
set.seed(0)
htree = hclust(dist(normz[rownames(annot)[annot$perturbed], ]))
cut = cutree(htree, k = 2)
annot$pclass = annot$covid
annot[names(cut), "pclass"] = paste0("p", cut)

cols$pclass = cols$covid
cols$pclass["p1"] = "red"
cols$pclass["p2"] = "darkviolet"

# save results:
hmdf = annot[, "pclass", drop = F]
colnames(hmdf) = "Cluster"
save(hmdf, file = "processed data/perturbation clustering results.RData")


# longitudinal plot showing perturbation clusters
svg("results/fig 2a perturbation longplot with cluster colors.svg", height = 4)
par(mar = c(5,5,0.2,0.1))
longplot(x = annot$perturbation, ylab = "Perturbation score", yaxt = T, col = cols$pclass[annot$pclass]) 
rect(40,thresh.perturb-0.03,160,3.05, col = NA, border = alpha("red", 0.5))
text(95, 2.93, "Samples designated\n\"highly perturbed\"", cex = .9, col = "darkred")
dev.off()




#### heatmap of only the most perturbed genes: -----------------
# get the most pertubed genes:
meanzlist = meandeltalist = list()
pclass.genes = list()
for (name in c("p1", "p2")) {
  use = annot$pclass == name
  meanzlist[[name]] = meanz = colMeans(normz[use, ])
  meandeltalist[[name]] = meanz = colMeans(norm[use, ]) - colMeans(norm[annot$covid == "HD", ])
  meanabsz = colMeans(abs(normz[use, ]))
  #pclass.genes[[name]] = names(meanabsz)[meanabsz > 3]
  pclass.genes[[name]] = c(names(meanabsz)[order(meanz, decreasing = T)[1:10]],
                           names(meanabsz)[order(meanz, decreasing = F)[1:10]])
}
topn = 20
mostperturbedgenes = unique(c(names(meanzlist[[1]])[order(meanzlist[[1]], decreasing = T)[1:topn]],
                              names(meanzlist[[1]])[order(meanzlist[[2]], decreasing = T)[1:topn]],
                              names(meanzlist[[1]])[order(meanzlist[[1]], decreasing = F)[1:topn]],
                              names(meanzlist[[1]])[order(meanzlist[[2]], decreasing = F)[1:topn]]))
hmcols = list(Cluster = c())
hmcols$Cluster["p1"] = "red"
hmcols$Cluster["p2"] = "darkviolet"

# draw heatmap:
svg("results/fig 2b heatmap of most up or down regulated genes in P1 and P2.svg", width = 4, height = 8)
pheatmap(t(normz[rownames(annot)[annot$perturbed], mostperturbedgenes]),
         breaks = seq(-6,6,length.out = 101),
         col = colorRampPalette(c("darkblue", "white", "darkred"))(100),
         show_rownames = T, fontsize_row = 7, show_colnames = F, annotation_legend = F,
         annotation_col = hmdf[, "Cluster", drop = F], annotation_colors = hmcols)
dev.off()

# legend
svg("results/clusters color legend.svg", width = 1.2, height = 1.5)
par(mar = c(0,0,0,0))
frame()
legend("center", pch = 16,  col = cols$pclass, legend = names(cols$pclass), cex = 1.5)
dev.off()

#### longitudinal plots of most interesting genes --------------------------

showgenes = c("PLAU", "IL1B", #"CXCL3", "VEGFA", # up in P1
              "NFKB1", "PLEK","LCP2", # up in p1, dn in p2
              #"SAMHD1", "ATG3", # dn in both 
              #"CXCL2", "OSM", "ICOSLG", # up in both
              "IRF3", # dn in p1, up in p2
              "MTOR", "IL18BP", "RACK1", 'TGFB1') # up in p2

svg("results/fig 2c perturbed genes longitudinal plots.svg", height = 10, width = 6)
layout(mat = matrix(c(1:5,11,6:10,11), 6), heights = c(rep(1,10), 0.1, 0.1))
par(mar = c(2,4.2,3,1))
for (gene in showgenes) {
  longplot(x = norm[, gene], ylab = gene, col = cols$pclass[annot$pclass]) 
}
par(mar = c(0,0,0,0))
frame()
legend("top", "Days since onset of symptoms", bty = "n", cex = 1.5)
dev.off()




if (FALSE) {  # explorations outside of manuscript below
  
  #### TCR score vs. perturbation class:
  longplot(annot$tcr, col = cols$pclass[hmdf$Cluster], ylab = "TCR score")
  fac = factor(hmdf$Cluster, levels = c("HD", "CCD", "p1", "p2"))
  bp = boxplot(annot$tcr ~ fac, outline = F, ylim = range(annot$tcr, na.rm = T), ylab = "TCR score", cex.lab = 1.5, cex.axis = 1.3, xlab = "")
  points(jitter(as.numeric(fac)), annot$tcr, col = cols$pclass[hmdf$Cluster], pch = 16)
  mod = summary(lm(annot$tcr ~ hmdf$Cluster))$coef
  text(3,3.025,paste0("p = ", signif(mod["hmdf$Clusterp1", 4], 3)))
  text(4,3.025,paste0("p = ", signif(mod["hmdf$Clusterp2", 4], 3)))
  
  
  
  #### antibodies vs. perturbation class: --------------------
  temp = lm(as.numeric(annot$IgG.Anti.SARS.CoV.2..ORTHO.) ~ hmdf$Cluster + annot$Sex + annot$RACE2 + annot$Age)
  use = as.numeric(names(temp$fitted.values))
  mod1 = lmer(as.numeric(annot$IgG.Anti.SARS.CoV.2..ORTHO.)[use] ~ hmdf$Cluster[use] + annot$Sex[use] + annot$RACE2[use] + annot$Age[use] + (1|annot$Unique.De.ID[use]))
  summary(mod1)$coef
  mod2 = lmer(as.numeric(annot$Total.Anti.SARS.CoV.2..ORTHO.)[use] ~ hmdf$Cluster[use] + annot$Sex[use] + annot$RACE2[use] + annot$Age[use] + (1|annot$Unique.De.ID[use]))
  summary(mod2)$coef
  
  mod1 = lm(as.numeric(annot$IgG.Anti.SARS.CoV.2..ORTHO.) ~ hmdf$Cluster + annot$donation.h)
  summary(mod1)$coef
  mod2 = lm(as.numeric(annot$Total.Anti.SARS.CoV.2..ORTHO.) ~ hmdf$Cluster + annot$donation.h)
  summary(mod2)$coef
  
  plot(as.numeric(annot$IgG.Anti.SARS.CoV.2..ORTHO.) ~ annot$donation.h, col = cols$pclass[hmdf$Cluster],pch = 16)
  plot(as.numeric(annot$Total.Anti.SARS.CoV.2..ORTHO.) ~ annot$donation.h, col = cols$pclass[hmdf$Cluster],pch = 16)
  
  # heatmap of all after subclustering perturbed samples:
  pheatmap(t(normz[, unique(unlist(perturbgenes))]),
           breaks = seq(-6,6,length.out = 101),
           col = colorRampPalette(c("darkblue", "white", "darkred"))(100),
           show_rownames = T, fontsize_row = 4, show_colnames = F,
           annotation_col = annot[, "pclass", drop = F], annotation_colors = cols)
  
  #### id gene sets in perturbed samples:
  meanzlist = meandeltalist = list()
  pclass.genes = list()
  for (name in c("p1", "p2")) {
    use = annot$pclass == name
    meanzlist[[name]] = meanz = colMeans(normz[use, ])
    meandeltalist[[name]] = meanz = colMeans(norm[use, ]) - colMeans(norm[annot$covid == "HD", ])
    meanabsz = colMeans(abs(normz[use, ]))
    #pclass.genes[[name]] = names(meanabsz)[meanabsz > 3]
    pclass.genes[[name]] = c(names(meanabsz)[order(meanz, decreasing = T)[1:10]],
                             names(meanabsz)[order(meanz, decreasing = F)[1:10]])
  }
  
  
  
  showgenes = (abs(meanzlist[[1]]) + abs(meanzlist[[2]]) > 5)
  plot(meanzlist[[1]], meanzlist[[2]], 
       xlab = "Change in average Z-score from HD", ylab = "Mean delta z in p2", cex.lab = 1.25,
       col = c(alpha("dodgerblue4", 0.5), "white")[1 + showgenes], pch = 16, cex = 0.75, asp = 1)
  text(meanzlist[[1]][showgenes], meanzlist[[2]][showgenes], names(meanzlist[[1]])[showgenes], cex = 0.5)
  abline(h = 0, col = alpha("grey50", 0.75))
  abline(v = 0, col = alpha("grey50", 0.75))
  
  
  showgenes = (abs(meandeltalist[[1]]) + abs(meandeltalist[[2]]) > 2)
  plot(meandeltalist[[1]], meandeltalist[[2]], 
       xlab = "Cluster P1 mean log2 fold-change in from HD", 
       ylab = "Cluster P1 mean log2 fold-change in from HD", 
       cex.lab = 1.25,
       col = c(alpha("dodgerblue4", 0.5), "white")[1 + showgenes], pch = 16, cex = 0.75, asp = 1)
  text(meandeltalist[[1]][showgenes], meandeltalist[[2]][showgenes], names(meandeltalist[[1]])[showgenes], cex = 0.5)
  abline(h = 0, col = alpha("grey50", 0.75))
  abline(v = 0, col = alpha("grey50", 0.75))
  
  
}




