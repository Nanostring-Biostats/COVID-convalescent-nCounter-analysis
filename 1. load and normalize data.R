rm(list = ls())
library(NormqPCR)
library(pheatmap)
library(readxl)

#### data loading ---------------------------
source("RCC_gather_n_collect.r")

## load and parse RCCs:
dat = collate.lanes.2.dataset(allRCC.path = "data/Host Response panel RCCs/")
dat = dat$datasets$ds1
raw = dat$count
# sample annotations:
annot = dat$sample.annot
# gene annotations:
gannot = dat$gannot
rownames(gannot) = colnames(raw)
rm(dat)


## load clinical annot (covid and healthy in different sheets)
clin.covid = read.csv("data/NanoString IDs 07272021 CCD formatted as dates.csv")
clin.healthy = read.csv("data/Nanostring IDs 07272021 HD.csv")

# harmonize variable names:
clin.healthy$Sex = clin.healthy$Gender
clin.healthy$RACE = clin.healthy$Race

## merge the clinical datasets:
# shared variables only:
clin = rbind(clin.covid[, intersect(colnames(clin.covid), colnames(clin.healthy))], 
             clin.healthy[, intersect(colnames(clin.covid), colnames(clin.healthy))])
# now add other variables available for covid samples:
covidonlycolumns = setdiff(colnames(clin.covid), colnames(clin.healthy))
coviddf = rbind(clin.covid[, covidonlycolumns], 
                matrix(NA, nrow = nrow(clin.healthy), ncol = length(covidonlycolumns), 
                       dimnames = list(rownames(clin.healthy), covidonlycolumns)))
clin = cbind(clin, coviddf) 
clin = clin[match(rownames(annot), make.names(clin$Nanostring_ID)), ]

# parse dates:
clin$Date.of.Blood.Collection = as.Date(clin$Date.of.Blood.Collection, format = "%m/%d/%y")
clin$Date.of.Diagnostic = as.Date(clin$Date.of.Diagnostic, format = "%m/%d/%y")
clin$Date.of.symptoms.onset = as.Date(clin$Date.of.symptoms.onset, format = "%m/%d/%y")
clin$Date.of.symptoms.resolution = as.Date(clin$Date.of.symptoms.resolution, format = "%m/%d/%y")


# report lists of patients with missing data:
print(paste0(setdiff(clin$Unique.De.ID[is.na(clin$Date.of.symptoms.onset)], NA), collapse = ", "))
print(paste0(setdiff(clin$Unique.De.ID[is.na(clin$Date.of.Diagnostic)], NA), collapse = ", "))

# get time variables:
clin$diagnosis = as.numeric(clin$Date.of.Diagnostic - clin$Date.of.symptoms.onset)
clin$donation = as.numeric(clin$Date.of.Blood.Collection - clin$Date.of.symptoms.onset)
clin$resolution = as.numeric(clin$Date.of.symptoms.resolution - clin$Date.of.symptoms.onset)

# define CCD vs. HD:
clin$covid = NA
clin$covid[is.element(clin$Nanostring_ID, clin.covid$Nanostring_ID)] = "CCD"
clin$covid[is.element(clin$Nanostring_ID, clin.healthy$Nanostring_ID)] = "HD"

#stopifnot(all(!is.na(clin$donation[clin$covid == "CCD"])))
table(is.na(clin$donation), clin$covid)

# for the purpose of plotting, a new time-since-donation variable with values for healthy:
clin$donation.h = clin$donation
clin$donation.h[clin$covid == "HD"] = 275

# condense the race variable:
clin$RACE2 = clin$RACE
clin$RACE2[clin$RACE2 == "White"] = "Caucasian"
clin$RACE2[clin$RACE2 == "Declined"] = "Other"
clin$RACE2[clin$RACE2 == "Mixed Race"] = "Other"
clin$RACE2[clin$RACE2 == "Unknown"] = "Other"
clin$RACE2[is.na(clin$RACE2)] = "Other"
table(clin$RACE2)

clin$Clin2 = clin$Clinical.Course
clin$Clin2[(clin$Clin2 == "MILD ") & !is.na(clin$Clin2)] = "MILD"
clin$Clin2[(clin$Clin2 == "ASYM") & !is.na(clin$Clin2)] = "MILD"
table(clin$Clin2)

annot = cbind(annot, clin)
annot$Unique.De.ID[annot$covid == "HD"] = paste0("HD", 1:sum(annot$covid == "HD"))

# define coloring:
cols = list()
cols$covid = c("HD" = "dodgerblue3", "CCD" = "orange")


#### add tcr results to annot: -------------------------


# load TCR data:
tcr = read.csv("data/TCRscores.csv", row.names = 1)
# align to annot:
m = match(rownames(annot),
          substr(make.names(rownames(tcr)), 1, nchar(make.names(rownames(tcr))) - 4))
annot$tcr = tcr$TCRScore[m]


#### QC --------------------------------

# calc signal factors:
hks = rownames(gannot)[gannot$CodeClass == "Housekeeping"]
annot$hkgeomean = exp(rowMeans(log(raw[, hks])))
annot$posgeomean = exp(rowMeans(log(raw[, rownames(gannot)[gannot$CodeClass == "Positive"]])))
annot$negmean = rowMeans(raw[, rownames(gannot)[gannot$CodeClass == "Negative"]])


# QC samples for low signal:
plot(annot$hkgeomean, annot$negmean)
plot(annot$posgeomean, annot$negmean)

hist(annot$hkgeomean, breaks = 40, xlim = c(0,max(annot$hkgeomean)), 
     xlab = "Housekeeper geomean", main = "", ylab = "Number of samples", cex.lab = 1.5)
thresh.hks = 200
abline(v = thresh.hks, lty = 2, col = 2)
legend("topleft", lty = 2, col = 2, legend = "Threshold for removal")

keep.samples = rownames(raw)[annot$hkgeomean > thresh.hks]
print(setdiff(rownames(raw), keep.samples))

# QC genes for low signal:
signal.to.bg = sweep(raw, 1, annot$negmean, "/")
hist(log2(apply(signal.to.bg, 2, mean)))
hist(log2(apply(signal.to.bg, 2, max)))
paste0(names(which(apply(signal.to.bg, 2, max) < 2)), collapse = ", ")
# a handful of genes have low signal, but we're keeping all of them

# keep only endogenous genes:
keep.genes = rownames(gannot)[gannot$CodeClass == "Endogenous"]

# remove bad RCCs:
annot = annot[keep.samples, ]
raw = raw[keep.samples, keep.genes]



#### Normalize ---------------------------------

# linear-scale normalized data:
normlin = sweep(raw, 1, annot$hkgeomean, "/") * 500
# log2-scale normalized data:
norm = log2(pmax(normlin, 1))


#### load TCR expression data (pre-normalized as part of NanoString analysis service): -----------------------------

# load TCR data:
trvs = read.csv("data/TSCOL_0223_NIHTCR_norm_data.csv", row.names = 1)
trvs = trvs[grepl("Endogenous", rownames(trvs)), ]
rownames(trvs) = gsub("Endogenous_", "", rownames(trvs))
trvs = t(trvs)
rownames(trvs) = gsub(".RCC", "", rownames(trvs))

# remove non-TCR genes:
trvs = trvs[, colnames(trvs)[grepl("TR", colnames(trvs))]]

# save all outputs:
write.csv(raw, file = "processed data/raw count data.csv")
write.csv(norm, file = "processed data/log2 scale normalized data.csv")
write.csv(normlin, file = "processed data/linear scale normalized data.csv")
write.csv(trvs, file = "processed data/TCR normalized data.csv")
write.csv(annot, file = "processed data/sample annotations.csv")
write.csv(gannot, file = "processed data/gene annotations.csv")
save(cols, file = "processed data/cols.RData")

dev.off()


