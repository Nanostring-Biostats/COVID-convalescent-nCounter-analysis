#' Run all DE analyses for a normalized data matrix (e.g. host response panel data, or TCR panel data)
#' @param norm Matrix of normalized data, to be analyzed on log2-scale
#' @param annot Annotation data frame
runDE = function(mat, annot, lookattcr = TRUE) {
  
  # harmonize IDs:
  sharedids = intersect(rownames(mat), rownames(annot))
  annot = annot[sharedids, ]
  mat = mat[sharedids, ]

  # identify each patient's earliest samples:
  firstsample = c()
  for (pt in unique(annot$Unique.De.ID)) {
    firstsample = c(
      firstsample,
      which(annot$Unique.De.ID == pt)[
        order(annot$donation[annot$Unique.De.ID == pt])[1]
      ])
  }
  
  # prep output:
  allde = data.frame(gene = colnames(mat))
  
  #### Differential expression: HD vs CCD: ----------------------------------------
  
  ests = ps = c()
  temp = cbind(annot[firstsample, ], mat[firstsample, ])
  for (gene in colnames(mat)) {
    mod = lm(data = temp, get(gene) ~ covid + Age + Sex + RACE2)
    ests[gene] = -summary(mod)$coefficients["covidHD", "Estimate"]
    ps[gene] = summary(mod)$coefficients[2, "Pr(>|t|)"]
  }
  fdrs = p.adjust(ps, method = "BH")
  
  # append to DE results:
  allde$'CCD/HD log2 fold-change' = ests[allde$gene]
  allde$'CCD/HD p-value' = ps[allde$gene]
  allde$'CCD/HD FDR' = fdrs[allde$gene]
  

  
  
  #### DE: HD vs CCD, by time window: -----------------------------------
  
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
  stratests = stratps = stratfdrs = list()
  for (name in names(timegrps)) {
    print(name)
    stratests[[name]] = stratps[[name]] = c()
    
    use = is.element(rownames(annot), timegrps[[name]]) | (annot$covid == 'HD')
    temp = cbind(annot[use,], mat[use, ])
    
    for (gene in colnames(mat)) {
      if (max(table(annot$Unique.De.ID[use])) > 1) {
        mod = lmer(data = temp, get(gene) ~ covid + Age + Sex + RACE2 + (1|Unique.De.ID))
        stratests[[name]][gene] = -summary(mod)$coefficients["covidHD", "Estimate"]
        stratps[[name]][gene] = summary(mod)$coefficients["covidHD", "Pr(>|t|)"]
      } else {
        mod = lm(data = temp, get(gene) ~ covid + Age + Sex + RACE2)
        stratests[[name]][gene] = -summary(mod)$coefficients["covidHD", "Estimate"]
        stratps[[name]][gene] = summary(mod)$coefficients["covidHD", "Pr(>|t|)"]
      }
    }
    stratfdrs[[name]] = p.adjust(stratps[[name]], method = "BH")
  }
  
  # append to allde:
  for (name in names(timegrps)) {
    allde[[paste0(name, "d CCD/HD log2 fold-change")]] = stratests[[name]][allde$gene]
    allde[[paste0(name, "d CCD/HD p-value")]] = stratps[[name]][allde$gene]
    allde[[paste0(name, "d CCD/HD FDR")]] = stratfdrs[[name]][allde$gene]
  }
  
  
  #### Differential expression: vs time: ----------------------------------------
  
  tempp = c()
  tempest = c()
  use = (annot$covid == "CCD") & !is.na(annot$donation)
  tempannot = annot[use, ]
  tempmat = mat[use, ]
  temp = cbind(annot[use, ], mat[use, ])
  for (gene in colnames(mat)) {
    mod = lmer(data = temp, get(gene) ~ donation + Age + Sex + RACE2 + (1 | Unique.De.ID))
    #mod = lmer(tempmat[, gene] ~ tempannot$donation + tempannot$Age + tempannot$Sex + tempannot$RACE2 + (1 | tempannot$Unique.De.ID))
    tempest[gene] = summary(mod)$coefficients[2,"Estimate"]
    tempp[gene] = summary(mod)$coefficients[2,"Pr(>|t|)"]
  }
  tempfdrs = p.adjust(tempp, method = "BH")
  min(tempfdrs)
  
  
  allde$'time log2 fold-change' = tempest[allde$gene]
  allde$'time p-value' = tempp[allde$gene]
  allde$'time FDR' = tempfdrs[allde$gene]
  
  
#### Differential expression vs. disease severity -----------------------
  
  tempannot = annot[firstsample, ]
  tempmat = mat[firstsample, ]
 
  # DE mods:
  est.sev.m = p.sev.m = est.sev.s = p.sev.s = c()
  for (gene in colnames(mat)) {
    mod = lm(tempmat[, gene] ~ tempannot$Clin2 + tempannot$Age + tempannot$Sex + tempannot$RACE2)
    est.sev.m[gene] = summary(mod)$coefficients["tempannot$Clin2MODERATE", "Estimate"]
    p.sev.m[gene] = summary(mod)$coefficients["tempannot$Clin2MODERATE", "Pr(>|t|)"]
    est.sev.s[gene] = summary(mod)$coefficients["tempannot$Clin2SEVERE", "Estimate"]
    p.sev.s[gene] = summary(mod)$coefficients["tempannot$Clin2SEVERE", "Pr(>|t|)"]
  }
  fdr.sev.m = p.adjust(p.sev.m, method = "BH")
  fdr.sev.s = p.adjust(p.sev.s, method = "BH")
  
  allde$'moderate/mild log2 fold change' = est.sev.m
  allde$'moderate/mild p-value' = p.sev.m
  allde$'moderate/mild fdr' = fdr.sev.m
  allde$'severe/mild log2 fold change' = est.sev.s
  allde$'severe/mild p-value' = p.sev.s
  allde$'severe/mild fdr' = fdr.sev.s
  

  #### DE of disease severity vs. healthy --------------------------------------
  
  ests.sevh = ps.sevh = fdrs.sevh = list()
  for (severity in c("MILD", "MODERATE", "SEVERE")) {
    use = (annot$Clin2 == severity)
    use[firstsample] = TRUE
    use[annot$covid == "HD"] = TRUE
    tempannot = annot[use, ]
    tempmat = mat[use, ]
    ests.sevh[[severity]] = ps.sevh[[severity]] = fdrs.sevh[[severity]]= c()
    for (gene in colnames(mat)) {
      mod = lm(tempmat[, gene] ~ tempannot$covid + tempannot$Age + tempannot$Sex + tempannot$RACE2)
      ests.sevh[[severity]][gene] = -summary(mod)$coefficients["tempannot$covidHD", "Estimate"]
      ps.sevh[[severity]][gene] = summary(mod)$coefficients["tempannot$covidHD", "Pr(>|t|)"]
    }
    fdrs.sevh[[severity]] = p.adjust(ps.sevh[[severity]], method = "BH")
  }
  
  for (severity in c("MILD", "MODERATE", "SEVERE")) {
    allde[[paste0(severity, "/healthy log2 fold change")]] = ests.sevh[[severity]]
    allde[[paste0(severity, "/healthy p-value")]] = ps.sevh[[severity]]
    allde[[paste0(severity, "/healthy FDR")]] = fdrs.sevh[[severity]]
  }
  
  
  #### association with demographic variables -----------------------------------
  
  demovars = c("Age","Sex","racecat")
  demops = matrix(NA, ncol(mat), length(demovars),
                  dimnames = list(colnames(mat), demovars))
  annot$racecat = as.factor(annot$RACE)
  annot$racecat = relevel(annot$racecat, "Caucasian")
  raceests = c()
  est.age = est.male = c()
  for (name in demovars) {
    for (gene in colnames(mat)) {
      mod = lm(mat[, gene] ~ annot[, name])
      demops[gene, name] = lrtest(mod)$Pr[2]   
      if (name == "Age") {
        est.age[gene] = mod$coefficients[2]
      }
      if (name == "Sex") {
        est.male[gene] = mod$coefficients[2]
      }
      if (name == "racecat") {
        raceests = rbind(raceests, mod$coefficients[-1])
      }
    }
  }
  
  demofdrs = demops*NA
  for (name in demovars) {
    demofdrs[, name] = p.adjust(demops[, name], method = "BH")
  }

  allde$'sex M/F log2 fold-change' = est.male
  allde$'sex p-value' = demops[, "Sex"]
  allde$'sex FDR' = demofdrs[, "Sex"]
  allde$'Age log2 fold change per year' = est.age
  allde$'Age p-value' = demops[, "Age"]
  allde$'Age FDR' = demofdrs[, "Age"]
  allde$'Race p-value' = demops[, "racecat"]
  allde$'Race FDR' = demofdrs[, "racecat"]


  #### correlation w ab titer -------------------
  
  # vs neutralizing antibody:
  x = NA
  x[is.element(annot$Neutralizing.Antibody..Holbrook., c("NONE", "None"))] = 0
  x[is.element(annot$Neutralizing.Antibody..Holbrook., c("<1:40","1:160","1:320","1:40","1:640","1:80"))] = 1
  use = !is.na(x)
  temp = cbind(annot[use, ], mat[use, ])
  temp$x = x[use]
  
  est.ab = p.ab = c()
  for (gene in colnames(mat)) {
    mod = lmer(data = temp, get(gene) ~ x + (1 | Unique.De.ID) + Age + Sex + RACE2)
    #mod = lmer(tempmat[, gene] ~ tempannot$x + (1 | tempannot$Unique.De.ID) + tempannot$Age + tempannot$Sex + tempannot$RACE2)
    est.ab[gene] = summary(mod)$coef[2,1]
    p.ab[gene] = summary(mod)$coef[2,"Pr(>|t|)"]
  }
  fdr.ab = p.adjust(p.ab, method = "BH")
  
  allde$'neutralizing ab > 0 log2 fold-change' = est.ab
  allde$'neutralizing ab > 0 p-value' = p.ab
  allde$'neutralizing ab > 0 FDR' = fdr.ab
  

  
  #### DE vs. annot$IgG.Anti.SARS.CoV.2..ORTHO. ------------------------------
  x = as.numeric(annot$IgG.Anti.SARS.CoV.2..ORTHO.) > 12
  use = !is.na(x)
  temp = cbind(annot[use, ], mat[use, ])
  temp$x = x[use]
  
  est.ab = p.ab = c()
  for (gene in colnames(mat)) {
    mod = lmer(data = temp, get(gene) ~ x + (1 | Unique.De.ID) + Age + Sex + RACE2)
    est.ab[gene] = summary(mod)$coef[2,1]
    p.ab[gene] = summary(mod)$coef[2,"Pr(>|t|)"]
  }
  fdr.ab = p.adjust(p.ab, method = "BH")
  
  
  allde$'ortho ab > 12 log2 fold-change' = est.ab
  allde$'ortho ab > 12 p-value' = p.ab
  allde$'ortho ab > 12 FDR' = fdr.ab
  
  
  #### DE vs. TCR diversity -----------------------------
  
  if (lookattcr) {
    
    use = (annot$covid == "CCD") & (!is.na(annot$tcr))
    temp = cbind(annot[use, ], mat[use, ])
    
    # TCR DE:
    est.tcr = p.tcr = c()
    for (gene in colnames(mat)) {
      mod = lmer(data = temp, get(gene) ~ tcr + (1 | Unique.De.ID) + Age + Sex + RACE2)
      est.tcr[gene] = summary(mod)$coef[2,1]
      p.tcr[gene] = summary(mod)$coef[2,"Pr(>|t|)"]
    }
    fdr.tcr = p.adjust(p.tcr, method = "BH")
    
    allde$'TCR log2 fold-change' = est.tcr
    allde$'TCR p-value' = p.tcr
    allde$'TCR FDR' = fdr.tcr
    
  }
  
  return(allde)
}

