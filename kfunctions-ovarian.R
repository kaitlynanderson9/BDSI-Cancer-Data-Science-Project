library(dplyr)
library(tidyverse)
library(MItools)
library(forcats)


ovKData <- ovarianData$ovarian_metadata

# create columns and calculate proportions of the cell types in each image

ovKData <- ovKData %>% left_join(
  ovarianData$ovarian_cells %>% group_by(sample_id) %>% 
    dplyr::summarize(
      totalCell = n(),
      p_CK = sum(pheno == "Tumor")/n(),
      p_CD8 = sum(pheno == "Cytotoxic T")/n(),
      p_CD19 = sum(pheno == "B Cell")/n(),
      p_CD14 = sum(pheno == "Macrophage")/n(),
      p_CD4 = sum(pheno == "T Helper")/n(),
    ), by = "sample_id"
)

# create columns for k and k cross values

ovKData$k_CK <- NA
ovKData$k_CD8 <- NA 
ovKData$k_CD19 <- NA
ovKData$k_CD14 <- NA
ovKData$k_CD4 <- NA
ovKData$k_CK_CD8 <- NA
ovKData$k_CK_CD19 <- NA
ovKData$k_CK_CD14 <- NA
ovKData$k_CK_CD4 <- NA
ovKData$k_CD19_CD8 <- NA
ovKData$k_CD19_CD14 <- NA
ovKData$k_CD19_CD4 <- NA
ovKData$k_CD8_CD14 <- NA
ovKData$k_CD8_CD4 <- NA
ovKData$k_CD14_CD4 <- NA

# calculate k tumor, k cytotoxic T, and k cross tumor x cytotoxic T

for(i in 1:nrow(ovKData)){
  
  tempData <- ovarianData$ovarian_cells %>% filter(
    sample_id == ovKData$sample_id[i]
  )
  
  pppObj <- spatstat.geom::ppp(
    x = tempData$x, y = tempData$y, 
    xrange = range(tempData$x), yrange = range(tempData$y),
    marks = factor(tempData$pheno)
  )
  
  if (ovKData$p_CK[i] > 0) {
    
    kEstCK_temp <- MItools::permute.envelope(pppObj, fun = Kest, funargs = list(i = "Tumor", correction = "isotropic"))
    
    ovKData$k_CK[i] <- (kEstCK_temp$obs - kEstCK_temp$mmean)[
      which.min(abs(kEstCK_temp$r - 40))
    ]
    
  } else {
    ovKData$k_CK[i] <- NA
  }
  
  if (ovKData$p_CD8[i] > 0) {
    
    kEstCD8_temp <- MItools::permute.envelope(pppObj, fun = Kest, funargs = list(i = "Cytotoxic T", correction = "isotropic"))
    
    ovKData$k_CD8[i] <- (kEstCD8_temp$obs - kEstCD8_temp$mmean)[
      which.min(abs(kEstCD8_temp$r - 40))
    ]
    
  } else {
    ovKData$k_CD8[i] <- NA
  }
  
  if (ovKData$p_CD8[i] > 0 & ovKData$p_CK[i] > 0) {
    
    kCross_temp <- MItools::permute.envelope(pppObj, fun = Kcross, funargs = list(i = "Tumor", j = "Cytotoxic T", correction = "isotropic"))
    
    ovKData$k_CK_CD8[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 40))
    ]
    
  } else {
    ovKData$k_CK_CD8[i] <- NA
  }
}

# calculate k B Cell, k Macrophage, and k T Helper

for(i in 6:6){
  
  tempData <- ovarianData$ovarian_cells %>% filter(
    sample_id == ovKData$sample_id[i]
  )
  
  pppObj <- spatstat.geom::ppp(
    x = tempData$x, y = tempData$y, 
    xrange = range(tempData$x), yrange = range(tempData$y),
    marks = factor(tempData$pheno)
  )
  
  if(ovKData$p_CD19[i] > 0){
    
    kEstCD19_temp <- MItools::permute.envelope(pppObj, fun = Kest, funargs = list(i = "B Cell", correction = "isotropic"))
    
    ovKData$k_CD19[i] <- (kEstCD19_temp$obs - kEstCD19_temp$mmean)[
      which.min(abs(kEstCD19_temp$r - 40))
    ]
  }
  else{
    
    ovKData$k_CD19[i] <- NA
  }
  
  if(ovKData$p_CD14[i] > 0){
    
    kEstCD14_temp <- MItools::permute.envelope(pppObj, fun = Kest, funargs = list(i = "Macrophage", correction = "isotropic"))
    
    ovKData$k_CD14[i] <- (kEstCD14_temp$obs - kEstCD14_temp$mmean)[
      which.min(abs(kEstCD14_temp$r - 40))
    ]
  }
  else{
    ovKData$k_CD14[i] <- NA
  }
  
  if(ovKData$p_CD4[i] > 0){
    
    kEstCD4_temp <- MItools::permute.envelope(pppObj, fun = Kest, funargs = list(i = "T Helper", correction = "isotropic"))
    
    ovKData$k_CD4[i] <- (kEstCD4_temp$obs - kEstCD4_temp$mmean)[
      which.min(abs(kEstCD4_temp$r - 40))
    ]
  }
  else{
    ovKData$k_CD4[i] <- NA
  }
  
}

# calculate k cross tumor x B cell, k cross tumor x macrophage, 
# and k cross tumor x T helper


for(i in 1:nrow(ovKData)){
  
  tempData <- ovarianData$ovarian_cells %>% filter(
    sample_id == ovKData$sample_id[i]
  )
  
  pppObj <- spatstat.geom::ppp(
    x = tempData$x, y = tempData$y, 
    xrange = range(tempData$x), yrange = range(tempData$y),
    marks = factor(tempData$pheno)
  )
  
  if(ovKData$p_CD19[i] > 0 & ovKData$p_CK[i] > 0){
    
    kCross_temp <- MItools::permute.envelope(pppObj, fun=Kcross, funargs=list(i="Tumor", j="B Cell", correction = "isotropic")
    )
    
    ovKData$k_CK_CD19[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 40))
    ]
  }
  else{
    
    ovKData$k_CK_CD19[i] <- NA
  }
  
  if(ovKData$p_CD14[i] > 0 & ovKData$p_CK[i] > 0){
    
    kCross_temp <- MItools::permute.envelope(pppObj, fun=Kcross, funargs=list(i="Tumor", j="Macrophage", correction = "isotropic")
    )
    
    ovKData$k_CK_CD14[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 40))
    ]
  }
  else{
    
    ovKData$k_CK_CD14[i] <- NA
  }
  
  if(ovKData$p_CD4[i] > 0 & ovKData$p_CK[i] > 0){
    
    kCross_temp <- MItools::permute.envelope(pppObj, fun=Kcross, funargs=list(i="Tumor", j="T Helper", correction = "isotropic")
    )
    
    ovKData$k_CK_CD4[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 40))
    ]
  }
  else{
    
    ovKData$k_CK_CD4[i] <- NA
  }
}

# calculate k cross cytotoxic T x B cell, k cross b cell x macrophage, 
# and k cross b cell x t helper

for(i in 1:nrow(ovKData)){
  
  tempData <- ovarianData$ovarian_cells %>% filter(
    sample_id == ovKData$sample_id[i]
  )
  
  pppObj <- spatstat.geom::ppp(
    x = tempData$x, y = tempData$y, 
    xrange = range(tempData$x), yrange = range(tempData$y),
    marks = factor(tempData$pheno)
  )
  
  if(ovKData$p_CD19[i] > 0 & ovKData$p_CD8[i] > 0){
    
    kCross_temp <- MItools::permute.envelope(pppObj, fun=Kcross, funargs=list(i="Cytotoxic T", j="B Cell", correction = "isotropic")
    )
    
    ovKData$k_CD19_CD8[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 40))
    ]
  }
  else{
    
    ovKData$k_CD19_CD8[i] <- NA
  }
  
  if(ovKData$p_CD19[i] > 0 & ovKData$p_CD14[i] > 0){
    
    kCross_temp <- MItools::permute.envelope(pppObj, fun=Kcross, funargs=list(i="B Cell", j="Macrophage", correction = "isotropic")
    )
    
    ovKData$k_CD19_CD14[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 40))
    ]
  }
  else{
    
    ovKData$k_CD19_CD14[i] <- NA
  }
  
  if(ovKData$p_CD4[i] > 0 & ovKData$p_CD19[i] > 0){
    
    kCross_temp <- MItools::permute.envelope(pppObj, fun=Kcross, funargs=list(i="B Cell", j="T Helper", correction = "isotropic")
    )
    
    ovKData$k_CD19_CD4[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 40))
    ]
  }
  else{
    
    ovKData$k_CD19_CD4[i] <- NA
  }
}

# calculate k cross cytotoxic t x macrophage, k cross cytotoxic t x t helper, 
# and k cross macrophage x t helper


for(i in 1:nrow(ovKData)){
  
  tempData <- ovarianData$ovarian_cells %>% filter(
    sample_id == ovKData$sample_id[i]
  )
  
  pppObj <- spatstat.geom::ppp(
    x = tempData$x, y = tempData$y, 
    xrange = range(tempData$x), yrange = range(tempData$y),
    marks = factor(tempData$pheno)
  )
  
  if(ovKData$p_CD14[i] > 0 & ovKData$p_CD8[i] > 0){
    
    kCross_temp <- MItools::permute.envelope(pppObj, fun=Kcross, funargs=list(i="Cytotoxic T", j="Macrophage", correction = "isotropic")
    )
    
    ovKData$k_CD8_CD14[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 40))
    ]
  }
  else{
    
    ovKData$k_CD8_CD14[i] <- NA
  }
  
  if(ovKData$p_CD8[i] > 0 & ovKData$p_CD4[i] > 0){
    
    kCross_temp <- MItools::permute.envelope(pppObj, fun=Kcross, funargs=list(i="Cytotoxic T", j="T Helper", correction = "isotropic")
    )
    
    ovKData$k_CD8_CD4[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 40))
    ]
  }
  else{
    
    ovKData$k_CD8_CD4[i] <- NA
  }
  
  if(ovKData$p_CD4[i] > 0 & ovKData$p_CD14[i] > 0){
    
    kCross_temp <- MItools::permute.envelope(pppObj, fun=Kcross, funargs=list(i="Macrophage", j="T Helper", correction = "isotropic")
    )
    
    ovKData$k_CD14_CD4[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 40))
    ]
  }
  else{
    
    ovKData$k_CD14_CD4[i] <- NA
  }
}

 




