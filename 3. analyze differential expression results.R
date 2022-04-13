
rm(list = ls())

library(ComplexHeatmap)
library(circlize)
library(viridis)
library(scales)
library(lmtest)
library(lme4)
library(lmerTest)
library(pheatmap)
library(np)
source("longitudinal plot function.R")


#### load data and DE results --------------------------

norm = as.matrix(read.csv("processed data/log2 scale normalized data.csv", row.names = 1))
gannot = read.csv("processed data/gene annotations.csv", row.names = 1)
annot = read.csv("processed data/sample annotations.csv", row.names = 1)
load("processed data/cols.RData")
allde = read.csv("results/DE results from host response panel.csv", row.names = 1, check.names = F)
allde = allde[colnames(norm), ]
probeannot = read.csv("data/ProbeAnnotations_NS_Hs_HostResponse_v1.0.csv")

#### fig 1a: volcano plots from each timepoint ---------------------------------
timenames = c("0-90d", "90-120d", "120-150d", "150-360d")
plottimenames = c("0-89 days", "90-119 days", "120-149 days", "150-241 days")
names(plottimenames) = timenames

svg("results/volcano plots per time window.svg", height = 4, width = 10)
layout(mat = matrix(1:4,1), widths = c(2.65,2,2,2))
par(mar = c(4,4,2,1))
for (name in timenames) {
  
  ests = allde[[paste0(name, " CCD/HD log2 fold-change")]]
  ps = allde[[paste0(name, " CCD/HD p-value")]]
  fdrs = allde[[paste0(name, " CCD/HD FDR")]]
  text.thresh = ps[order(ps)[20]]
  p.thresh = max(ps[fdrs <= 0.05])
  plot(ests, -log10(ps), xlab = "log2 fold-change", main = plottimenames[name], 
       ylim = c(0, -log10(min(allde[, paste0(timenames, " CCD/HD p-value")]))),
       xlim = range(allde[, paste0(timenames, " CCD/HD log2 fold-change")]),
       ylab = "-log10(p-value)", cex.lab = 1.5, pch = 16,
       yaxt = "n", 
       col = c(alpha("grey30", 0.5), alpha("white", 0.1))[1 + (ps < text.thresh)])
  abline(h = -log10(p.thresh), lty = 2)
  text(ests[ps < text.thresh], -log10(ps)[ps < text.thresh], rownames(allde)[ps < text.thresh], cex = 0.75,
       col = c("darkblue", "darkred")[1 + (ests[ps < text.thresh] > 0)])
  if (name == "0-90d") {
    axis(2)
  }
  par(mar = c(4,0,2,1))
}
legend("topright", lty = 2, legend = "FDR = 0.05")
dev.off()


# summarize extent of DE:
numsig = c()
for (name in timenames) {
  ests = allde[[paste0(name, " CCD/HD log2 fold-change")]]
  ps = allde[[paste0(name, " CCD/HD p-value")]]
  fdrs = allde[[paste0(name, " CCD/HD FDR")]]
  
  print(name)
  print(sum((abs(ests) > log2(1.2)) & (fdrs < 0.05)))
  numsig[name] = sum((abs(ests) > log2(1.2)) & (fdrs < 0.05))
}
write.csv(numsig, "results/number of significant genes per time window.csv")



#### fig 1b: heatmap of most DE genes -------------------------------

topgenes = c()
for (name in timenames) {
  ests = allde[[paste0(name, " CCD/HD log2 fold-change")]]
  ps = allde[[paste0(name, " CCD/HD p-value")]]
  fdrs = allde[[paste0(name, " CCD/HD FDR")]]
  
  print(name)
  print(paste0("up: ", paste0(rownames(allde)[(ests > log2(1.5)) & (fdrs < 0.05)], collapse = ", ")))
  print(paste0("down: ", paste0(rownames(allde)[(ests< -log2(1.5)) & (fdrs < 0.05)], collapse = ", ")))
  topgenes = unique(c(topgenes, rownames(allde)[(abs(ests) > log2(1.5)) & (fdrs < 0.05)]))
}
topgenes

# make summary table to top genes:
df = data.frame(gene = topgenes)
rownames(df) = df$gene
df$gene = gsub("\\.", "/", df$gene)
df$geneset = probeannot$Probe.Annotation[match(df$gene, probeannot$Probe.Label)]
# record which genesets are among our top genes:
genesets = c()
for (i in 1:nrow(df)) {
  genesets = unique(c(genesets, strsplit(df[i, "geneset"], ";")[[1]]))
}
# how often does each gene set appear in our gene list?
appearances = c()
for (name in genesets) {
  appearances[name] = sum(grepl(name, df$geneset))
}
genesets = names(appearances)[appearances >= 1]
# remove vague genesets:
genesets = setdiff(genesets, c("Chemokine Signaling",
                               "DNA Sensing", "RNA Sensing", "Other Interleukin Signaling"))
# add to df:
for (name in genesets) {
  df[[name]] = 1 * grepl(name, df$geneset)
}
# now record gene's DE results:
for (name in timenames) {
  ests = allde[[paste0(name, " CCD/HD log2 fold-change")]]
  fdrs = allde[[paste0(name, " CCD/HD FDR")]]
  
  #df[[name]] = (ests * (fdrs < 0.05))[match(rownames(df), rownames(allde))]
  df[[name]] = (ests)[match(rownames(df), rownames(allde))]
}

p0 = pheatmap(df[, genesets])
annotcols = list()
for (name in genesets) {
  annotcols[[name]] = c("0" = "white", "1" = "darkviolet")
}

selected.genesets = c(
  "Myeloid Activation", "TNF Signaling","NF-kappaB Signaling","JAK-STAT Signaling","Immune Exhaustion",
  "Type I Interferon Signaling","T-cell Costimulation","Complement System",
  "IL-1 Signaling", "IL-2 Signaling", "IL-6 Signaling", "Virus-Host Interaction"
)

dev.off()
svg("results/fig 1b heatmap.svg", width = 5.5, height = 6)
pheatmap(df[p0$tree_row$order, timenames],
         cluster_cols = F, cluster_rows = T,
         col = colorRampPalette(c("darkblue", "white", "darkred"))(100),
         breaks = seq(-1.5,1.5,length.out=101), 
         annotation_legend = F, annotation_colors = annotcols,
         annotation_row = df[, selected.genesets])
dev.off()



#### fig 1c: plots of selected genes ---------------------------

selectedgenes = c("CTLA4", "CXCR4",   # uniform drop
                  "OSM",    # delayed drop
                  "CXCL2", "CCL3.L1.L3", "IL1B", # many high, some normal
                  "IFNA6", 'HERC5')#,      # stays low

svg("results/selected genes over time.svg", height = 10, width = 6)
layout(mat = t(matrix(c(1:9,9), 2)), heights = c(rep(1,10), 0.1))
par(mar = c(2,4,3,1))
for (gene in selectedgenes) {
  longplot(x = norm[, gene], ylab = gene) 
}
par(mar = c(0,0,0,0))
frame()
legend("top", "Days since onset of symptoms", bty = "n", cex = 1.5)
dev.off()


for (gene in selectedgenes) {
  svg(paste0("results/expression vs. time - ", gene, ".svg"), height = 4)
  par(mar = c(4,4,1,1))
  longplot(x = norm[, gene], ylab = gene) 
  dev.off()
}




#### heatmap of loess fit on residuals from demographics ---------------------------------------------

if (FALSE) {  # these results were not in the manuscript
  # record how many genes were ever significant:
  everlowfdr = c()
  for (name in timenames) {
    everlowfdr = unique(c(everlowfdr, rownames(allde)[
      allde[[paste0(name, " CCD/HD FDR")]] < 0.05]))
  }
  
  
  # choose genes to show: anything meeting a significance threshold at any time window:
  eversig = c()
  for (name in timenames) {
    eversig = unique(c(eversig, rownames(allde)[
      (allde[[paste0(name, " CCD/HD FDR")]] < 0.05) & (abs(allde[[paste0(name, " CCD/HD log2 fold-change")]]) > log2(1.33))]))
  }
  
  # first, fit a model of only demographics data and get resids:
  resids = norm[, eversig] * NA
  for (gene in eversig) {
    mod = lm(norm[, gene] ~ annot$Sex + annot$RACE2 + annot$Age + annot$covid + annot$donation.h)
    # get fitted value based on demographics only:
    demovars = paste0("annot$", c("SexM", "RACE2Black", "RACE2Caucasian", "RACE2Hispanic", "RACE2Other", "Age"))
    X.demo = model.matrix(mod)[, demovars]
    yhat.demo = X.demo %*% mod$coefficients[demovars]
    tempresids = norm[rownames(yhat.demo), gene] - yhat.demo
    # center to have mean = 0 in HD:
    hdnames = rownames(annot)[annot$covid == "HD"]
    tempresids = tempresids - mean(tempresids[is.element(rownames(tempresids), hdnames), ])
    resids[rownames(tempresids), gene] = tempresids
  }
  
  # now fit loess models to residuals from demographics:
  use = setdiff(rownames(resids), hdnames)
  fitted = c()
  issig = c()
  times = seq(60,200,5)
  for (gene in colnames(resids)) {
    
    low = loess(resids[use, gene] ~ annot[use, "donation"])
    pr = predict(low, newdata = times, se = T)  
    
    issig.up = (pr$fit - qt(0.975,pr$df) * pr$se.fit > 0)
    issig.dn = (pr$fit + qt(0.975,pr$df) * pr$se.fit < 0) 
    
    fitted = rbind(fitted, pr$fit)
    issig = rbind(issig, issig.up | issig.dn)
    
  }
  rownames(fitted) = rownames(issig) = colnames(resids)
  colnames(fitted) = colnames(issig) = times
  
  table(issig[,"180"])
  table(issig[,"200"])
  
  # get genesets to annotate heatmap:
  genesets = c()
  for (i in 1:nrow(df)) {
    genesets = unique(c(genesets, strsplit(df[i, "geneset"], ";")[[1]]))
  }
  #selected.genesets = setdiff(genesets, NA)
  selected.genesets = c(
    "Myeloid Activation", "TNF Signaling","NF-kappaB Signaling","JAK-STAT Signaling","Immune Exhaustion",
    "Type I Interferon Signaling","T-cell Costimulation","Complement System",
    "IL-1 Signaling", "IL-2 Signaling", "IL-6 Signaling", "Virus-Host Interaction",
    "Lysosome","MHC Class II Antigen Presentation"
  )
  df = data.frame(gene = eversig)
  rownames(df) = df$gene
  df$gene = gsub("\\.", "/", df$gene)
  df$geneset = probeannot$Probe.Annotation[match(df$gene, probeannot$Probe.Label)]
  for (name in selected.genesets) {
    df[[name]] = 1 * grepl(name, probeannot[match(rownames(df), make.names(probeannot$Probe.Label)), "Probe.Annotation"])
  }
  annotcols = list()
  for (name in genesets) {
    annotcols[[name]] = c("0" = "white", "1" = "darkviolet")
  }
  
  # shorten name:
  colnames(df)[colnames(df) == "MHC Class II Antigen Presentation"] = "MHCII Antigen Presentation"
  annotcols[["MHCII Antigen Presentation"]] = c("0" = "white", "1" = "darkviolet")
  annotcols[["MHC Class II Antigen Presentation"]] = NULL
  selected.genesets = setdiff(c(selected.genesets, "MHCII Antigen Presentation"), "MHC Class II Antigen Presentation")
  
  p1 = pheatmap((fitted * issig)[eversig, ])
  dev.off()
  svg("results/loess heatmap - resids - eversig.svg", width = 6.5, height = 8) #, units = "in", res = 600)
  pheatmap((fitted * issig)[eversig, ][p1$tree_row$order, ], cluster_rows = F, 
           col = colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(100),
           breaks = seq(-2, 2, length.out = 101), 
           cluster_cols = F, fontsize_row = 6, fontsize_col = 8,
           #main = "log2 fold-change from HD",
           annotation_row = df[, selected.genesets], annotation_colors = annotcols, annotation_legend = F)
  dev.off()
  
  
}
