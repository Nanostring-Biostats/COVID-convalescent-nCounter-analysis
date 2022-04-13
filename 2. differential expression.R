
rm(list = ls())

library(lmtest)
library(lme4)
library(lmerTest)

source("differential expression pipeline function.R")


#### load data ------------------------------
norm = as.matrix(read.csv("processed data/log2 scale normalized data.csv", row.names = 1))
gannot = read.csv("processed data/gene annotations.csv", row.names = 1)
annot = read.csv("processed data/sample annotations.csv", row.names = 1)
trvs = as.matrix(read.csv("processed data/TCR normalized data.csv", row.names = 1))


#### run global DE function --------------------------

# DE of normalized data
de.norm = runDE(mat = norm, annot = annot, lookattcr = TRUE)
write.csv(de.norm, file = "results/DE results from host response panel.csv", row.names = F)
# condensed results of just analyses presented in manuscript:
shownames = c('gene', colnames(de.norm)[grepl("d CCD/HD", colnames(de.norm))])
write.csv(de.norm[, shownames], file = "results/Supp Table 3 - DE results from host response panel.csv", row.names = F)


# DE of TCR diversity score:
de.tcr = runDE(mat = annot[, c("tcr", "tcr")], annot = annot, lookattcr = FALSE)
write.csv(de.tcr[1, , drop = F], file = "results/DE results from TCR score.csv", row.names = F)

de.trvs = runDE(mat = trvs, annot = annot, lookattcr = TRUE)
write.csv(de.trvs, file = "results/DE results from TCR panel.csv", row.names = F)
