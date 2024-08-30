
library(dplyr)
library(tidyverse)
library(MItools)
library(forcats)


lungKData <- lungData$lung_metadata

# create columns and calculate proportions of the cell types in each image

lungKData <- lungKData %>% left_join(
  lungData$lung_cells %>% group_by(image_id) %>% 
    dplyr::summarize(
      totalCell = n(),
      p_CK = sum(pheno == "CK")/n(),
      p_CD8 = sum(pheno == "CD8")/n(),
      p_CD19 = sum(pheno == "CD19")/n(),
      p_CD14 = sum(pheno == "CD14")/n(),
      p_CD4 = sum(pheno == "CD4")/n(),
    ), by = "image_id"
)

# create columns for k and k cross values

lungKData$k_CK <- NA
lungKData$k_CD8 <- NA 
lungKData$k_CD19 <- NA
lungKData$k_CD14 <- NA
lungKData$k_CD4 <- NA
lungKData$k_CK_CD8 <- NA
lungKData$k_CK_CD19 <- NA
lungKData$k_CK_CD14 <- NA
lungKData$k_CK_CD4 <- NA
lungKData$k_C19_CD8 <- NA
lungKData$k_C19_CD14 <- NA
lungKData$k_C19_CD4 <- NA
lungKData$k_C8_CD14 <- NA
lungKData$k_C8_CD4 <- NA
lungKData$k_C14_CD4 <- NA


# calculate k tumor, k cytotoxic T, and k cross tumor x cytotoxic T


for(i in 1:nrow(lungKData)){

  tempData <- lungData$lung_cells %>% filter(
    image_id == lungKData$image_id[i]
  )
  
  pppObj <- spatstat.geom::ppp(
    x = tempData$x, y = tempData$y, 
    xrange = range(tempData$x), yrange = range(tempData$y),
    marks = factor(tempData$pheno)
  )
  
  if(lungKData$p_CK[i] > 0){
    
    kEstCK_temp <- MItools::permute.envelope(pppObj, fun = Kest, funargs = list(i = "CK", correction = "isotropic"))
    
    lungKData$k_CK[i] <- (kEstCK_temp$obs - kEstCK_temp$mmean)[
      which.min(abs(kEstCK_temp$r - 40))
    ]
  }
  else{

    lungKData$k_CK[i] <- NA
  }
  
  if(lungKData$p_CD8[i] > 0){
    
    kEstCD8_temp <- MItools::permute.envelope(pppObj, fun = Kest, funargs = list(i = "CD8", correction = "isotropic"))
    
    lungKData$k_CD8[i] <- (kEstCD8_temp$obs - kEstCD8_temp$mmean)[
      which.min(abs(kEstCD8_temp$r - 40))
    ]
  }
  else{
    lungKData$k_CD8[i] <- NA
  }
  
  if(lungKData$p_CD8[i] > 0 & lungKData$p_CK[i] > 0){

    kCross_temp <- MItools::permute.envelope(pppObj, fun=Kcross, funargs=list(i="CK", j="CD8", correction = "isotropic")
    )
    
    lungKData$k_CK_CD8[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 40))
    ]
  }
  else{

    lungKData$k_CK_CD8[i] <- NA
  }
}

# calculate k B Cell, k macrophage, and k T Helper 


for(i in 1:nrow(lungKData)){
  
  tempData <- lungData$lung_cells %>% filter(
    image_id == lungKData$image_id[i]
  )
  
  pppObj <- spatstat.geom::ppp(
    x = tempData$x, y = tempData$y, 
    xrange = range(tempData$x), yrange = range(tempData$y),
    marks = factor(tempData$pheno)
  )
  
  if(lungKData$p_CD19[i] > 0){
    
    kEstCD19_temp <- MItools::permute.envelope(pppObj, fun = Kest, funargs = list(i = "CD19", correction = "isotropic"))
    
    lungKData$k_CD19[i] <- (kEstCD19_temp$obs - kEstCD19_temp$mmean)[
      which.min(abs(kEstCD19_temp$r - 40))
    ]
  }
  else{
    
    lungKData$k_CD19[i] <- NA
  }
  
  if(lungKData$p_CD14[i] > 0){
    
    kEstCD14_temp <- MItools::permute.envelope(pppObj, fun = Kest, funargs = list(i = "CD14", correction = "isotropic"))
    
    lungKData$k_CD14[i] <- (kEstCD14_temp$obs - kEstCD14_temp$mmean)[
      which.min(abs(kEstCD14_temp$r - 40))
    ]
  }
  else{
    lungKData$k_CD14[i] <- NA
  }
  
  if(lungKData$p_CD4[i] > 0){
    
    kEstCD4_temp <- MItools::permute.envelope(pppObj, fun = Kest, funargs = list(i = "CD4", correction = "isotropic"))
    
    lungKData$k_CD4[i] <- (kEstCD4_temp$obs - kEstCD4_temp$mmean)[
      which.min(abs(kEstCD4_temp$r - 40))
    ]
  }
  else{
    lungKData$k_CD4[i] <- NA
  }
  
}

# calculate k cross tumor x B cell, k cross tumor x macrophage, 
# and k cross tumor x T helper


for(i in 1:nrow(lungKData)){
  
  tempData <- lungData$lung_cells %>% filter(
    image_id == lungKData$image_id[i]
  )
  
  pppObj <- spatstat.geom::ppp(
    x = tempData$x, y = tempData$y, 
    xrange = range(tempData$x), yrange = range(tempData$y),
    marks = factor(tempData$pheno)
  )
  
  if(lungKData$p_CD19[i] > 0 & lungKData$p_CK[i] > 0){
    
    kCross_temp <- MItools::permute.envelope(pppObj, fun=Kcross, funargs=list(i="CK", j="CD19", correction = "isotropic")
    )
    
    lungKData$k_CK_CD19[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 40))
    ]
  }
  else{
    
    lungKData$k_CK_CD19[i] <- NA
  }
  
  if(lungKData$p_CD14[i] > 0 & lungKData$p_CK[i] > 0){
    
    kCross_temp <- MItools::permute.envelope(pppObj, fun=Kcross, funargs=list(i="CK", j="CD14", correction = "isotropic")
    )
    
    lungKData$k_CK_CD14[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 40))
    ]
  }
  else{
    
    lungKData$k_CK_CD14[i] <- NA
  }
  
  if(lungKData$p_CD4[i] > 0 & lungKData$p_CK[i] > 0){
    
    kCross_temp <- MItools::permute.envelope(pppObj, fun=Kcross, funargs=list(i="CK", j="CD4", correction = "isotropic")
    )
    
    lungKData$k_CK_CD4[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 40))
    ]
  }
  else{
    
    lungKData$k_CK_CD4[i] <- NA
  }
}

# calculate k cross cytotoxic T x B cell, k cross b cell x macrophage, 
# and k cross b cell x t helper


lungKData$k_CD19_CD8 <- NA
lungKData$k_CD19_CD14 <- NA
lungKData$k_CD19_CD4 <- NA
lungKData$k_CD8_CD14 <- NA
lungKData$k_CD8_CD4 <- NA
lungKData$k_CD14_CD4 <- NA

for(i in 1:nrow(lungKData)){
  
  tempData <- lungData$lung_cells %>% filter(
    image_id == lungKData$image_id[i]
  )
  
  pppObj <- spatstat.geom::ppp(
    x = tempData$x, y = tempData$y, 
    xrange = range(tempData$x), yrange = range(tempData$y),
    marks = factor(tempData$pheno)
  )
  
  if(lungKData$p_CD19[i] > 0 & lungKData$p_CD8[i] > 0){
    
    kCross_temp <- MItools::permute.envelope(pppObj, fun=Kcross, funargs=list(i="CD8", j="CD19", correction = "isotropic")
    )
    
    lungKData$k_CD19_CD8[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 40))
    ]
  }
  else{
    
    lungKData$k_C19_CD8[i] <- NA
  }
  
  if(lungKData$p_CD19[i] > 0 & lungKData$p_CD14[i] > 0){
    
    kCross_temp <- MItools::permute.envelope(pppObj, fun=Kcross, funargs=list(i="CD19", j="CD14", correction = "isotropic")
    )
    
    lungKData$k_CD19_CD14[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 40))
    ]
  }
  else{
    
    lungKData$k_CD19_CD14[i] <- NA
  }
  
  if(lungKData$p_CD4[i] > 0 & lungKData$p_CD19[i] > 0){
    
    kCross_temp <- MItools::permute.envelope(pppObj, fun=Kcross, funargs=list(i="CD19", j="CD4", correction = "isotropic")
    )
    
    lungKData$k_CD19_CD4[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 40))
    ]
  }
  else{
    
    lungKData$k_CD19_CD4[i] <- NA
  }
}

# calculate k cross cytotoxic t x macrophage, k cross cytotoxic t x t helper, 
# and k cross macrophage x t helper


for(i in 1:nrow(lungKData)){
  
  tempData <- lungData$lung_cells %>% filter(
    image_id == lungKData$image_id[i]
  )
  
  pppObj <- spatstat.geom::ppp(
    x = tempData$x, y = tempData$y, 
    xrange = range(tempData$x), yrange = range(tempData$y),
    marks = factor(tempData$pheno)
  )
  
  if(lungKData$p_CD14[i] > 0 & lungKData$p_CD8[i] > 0){
    
    kCross_temp <- MItools::permute.envelope(pppObj, fun=Kcross, funargs=list(i="CD8", j="CD14", correction = "isotropic")
    )
    
    lungKData$k_CD8_CD14[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 40))
    ]
  }
  else{
    
    lungKData$k_CD8_CD14[i] <- NA
  }
  
  if(lungKData$p_CD8[i] > 0 & lungKData$p_CD4[i] > 0){
    
    kCross_temp <- MItools::permute.envelope(pppObj, fun=Kcross, funargs=list(i="CD8", j="CD4", correction = "isotropic")
    )
    
    lungKData$k_CD8_CD4[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 40))
    ]
  }
  else{
    
    lungKData$k_CD8_CD4[i] <- NA
  }
  
  if(lungKData$p_CD4[i] > 0 & lungKData$p_CD14[i] > 0){
    
    kCross_temp <- MItools::permute.envelope(pppObj, fun=Kcross, funargs=list(i="CD14", j="CD4", correction = "isotropic")
    )
    
    lungKData$k_CD14_CD4[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 40))
    ]
  }
  else{
    
    lungKData$k_CD14_CD4[i] <- NA
  }
}

columns_to_remove <- c("k_C19_CD8", "k_C19_CD14", "k_C19_CD4", 
                       "k_C8_CD14", "k_C8_CD4", "k_C14_CD4")

lungKData <- lungKData[, !names(lungKData) %in% columns_to_remove]













