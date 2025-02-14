
##### Set working directory
#setwd("/Users/dianaochoa/Desktop/data for R/CoEx_AGL jun 2023")
#setwd("G:/Mi unidad/STRI_CTPA/Colaborations_STRI/DianaOchoa/CoEx_2025_01")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   

# LOAD COEX FUNCTIONS
#source("/Users/dianaochoa/Desktop/data for R/CoEx_AGL jun 2023/Coexistence_functions_REV.R")
source("Coexistence_functions_REV.R")


# Load Packages
library(readxl)
library(BIEN)
library(raster)
library(KernSmooth)
library(geodata)
library(viridis)
library(ggplot2)
#library(tidyverse)
library(gridExtra)
library(reshape2)



# # # # Load Diana's data

# Taxa affinities
affs <- as.data.frame(read_excel("affinitiesMasterTable_Ochoa_only pollen.xlsx", sheet = 1))

# New modern records
modern <- as.data.frame(read_excel("1 Master Flora file ap2023.xlsx", sheet = 1))### Master file refers to DataFile S7 - Published  data
spl1 <- as.data.frame(read_excel("1 Master Flora file ap2023.xlsx", sheet = 2)) ### Master file refers to DataFile S7 - Herbaria data

# Species Link records
spl <- read.csv("speciesLink-20230503115019-0005596 rev.txt", sep = "\t")

# fossil counts data and formating

conteo <- as.data.frame(read_excel("AGL Coex counts_2023 only pollen.xlsx" ,col_names = FALSE)) 
# conteo <- as.data.frame(read_excel("AGL Coex counts_2023 only pollen noParacas.xlsx" ,col_names = FALSE)) 
conteo1 <- as.data.frame(t(conteo[-(1:2),-1]))
colnames(conteo1)<- conteo[-(1:2),1]
rownames(conteo1)<- conteo[2,-1]
time <- (conteo[1,-1])

fossil <- as.data.frame(lapply(conteo1,function (x) as.numeric(as.character(x))))
colnames(fossil)<- conteo[-(1:2),1]
rownames(fossil)<- conteo[2,-1]
# NA's = zero
fossil[is.na(fossil)] <- 0

# check fossil taxa in affinities

all(colnames(fossil) %in% affs$morphotype)
# must be TRUE

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # #  Coexistence analyses (step by step) # # # #

# Modern taxa affinities names
   mtx <- as.character(unique(affs$affinity))	

# Load BIEN occurrences (from downloaded data - the faster way -)
  load("BIEN_occurrences.RData")

# Compare which taxa are in the BIEN db
  newtx <- setdiff(mtx,names(BIEN_occs))
  bien_ls <- vector("list", length(newtx))
  names(bien_ls) <- newtx

# = = =# = = =# = = =# = = =# = = =
# Search new taxa in the BIEN web db ( just incase there are new records, this can be skiped to save time
# = = =# = = =# = = =# = = =
 for (j in which(lapply(bien_ls,length)== 0)){

	a <- BIEN_occurrence_genus(newtx[j],
  		cultivated = F ,
  		new.world = T,
  		all.taxonomy = F,
  		native.status = F,
  		natives.only = T,
  		observation.type = F,
  		political.boundaries = F ,
  		collection.info = F
		)
	if(nrow(a) > 0 ){	
		bien_ls[[j]] <- a[,c("latitude","longitude")]
	} 
    } #j
 # = = =# = = =# = = =# = = =
 
  # Eliminate empty taxa

   bien_ls[which(unlist(lapply(bien_ls,is.null)))] <- NULL

  # Keep only taxa with 3 o more records
   bien_ls <- bien_ls[lapply(bien_ls,nrow) >= 3]
  # combine data into the same object
   BIEN_occs <- c(BIEN_occs,bien_ls)

# = = =

#  include info from species link & records from Master Flora and 
# combine with BIEN data into modern occurrences from the modern affinities

maffs <- sort(unique(affs$affinity))
moccs <- list()

for(afi in 1:length(maffs)){
  #m0 <- modern[which(modern$Genus == maffs[afi]),c("Longitude","Latitude")]
    m0 <- modern[which(modern$Genus == maffs[afi] | modern$Family == maffs[afi]), c("Longitude","Latitude")]
  names(m0) <- c("longitude","latitude")
  
   h0 <- spl1[which(spl1$genus == maffs[afi] | spl$family == maffs[afi]),  c("longitude","latitude")] # h0 for Herbaria

  #l0 <- spl[which(spl$genus == maffs[afi]),c("longitude","latitude")]
  l0<- spl[which(spl$genus == maffs[afi] | spl$family == maffs[afi]),c("longitude","latitude")]
  
  b01 <- which(names(BIEN_occs) == maffs[afi])
  if(length(b01) == 0){
	b0 <- data.frame(longitude = NA, latitude= NA)
  } else{
      b02 <- BIEN_occs[[b01]]
      b0 <- na.omit(b02)
  }
  mb0 <- rbind(m0,h0, l0,b0)
  moccs[[afi]] <- na.omit(mb0)
 
  names(moccs)[[afi]] <- maffs[afi]
}

length(which(lapply(moccs,nrow) == 0))


which(lapply(moccs,nrow) < 3 & lapply(moccs,nrow) != 0)

# Eliminate records with less than 3 observations
moccs[which(lapply(moccs,nrow) < 3)] <- NULL

# 188 taxa with at least 3 observations
length(moccs) ### 48 ? - ?


# # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Extracting data from World Clim (geodata package)
  
#  bioclim BIO18  
wc.bio <- worldclim_global(var = "bio", res = 5, path = getwd())

#BIO1 = Annual Mean Temperature
#BIO5 = Max Temperature of Warmest Month
#BIO6 = Min Temperature of Coldest Month
#BIO7 = Temperature Annual Range (BIO5-BIO6)
#BIO8 = Mean Temperature of Wettest Quarter
#BIO9 = Mean Temperature of Driest Quarter
#BIO10 = Mean Temperature of Warmest Quarter
#BIO11 = Mean Temperature of Coldest Quarter
#BIO12 = Annual Precipitation
#BIO13 = Precipitation of Wettest Month
#BIO14 = Precipitation of Driest Month
#BIO15 = Precipitation Seasonality (Coefficient of Variation)##a statistical measure of the dispersion of data points around the mean
#BIO16 = Precipitation of Wettest Quarter
#BIO17 = Precipitation of Driest Quarter
#BIO18 = Precipitation of Warmest Quarter
#BIO19 = Precipitation of Coldest Quarter


# set climatic parameters

Var1 <- wc.bio[[1]] 	# Annual Mean Temperature
Var2 <- wc.bio[[12]] 	# Annual Precipitation

nameX <- "MAT"
nameY <- "MAP"

range_x <- c(0,30)
range_y <- c(0,1000)   ### PP restricted to  see the lower limit

 # select a bandwidth for  kernel distributions
bw1 <- abs(max(range_x)- min(range_x)) /20
bw2 <- abs(max(range_y)- min(range_y)) /20

###########################################################
# Afinidad climática por taxon moderno

mod.taxa.2d <- list()

for(moi in 1:length(moccs)){

  # extract complete climatic data for spatial occurrences
  moccsi <- moccs[[moi]][c("longitude","latitude")]
  moccsi$v1 <- extract(Var1, moccsi[,c("longitude","latitude")])[ ,2]
  moccsi$v2 <- extract(Var2, moccsi[,c("longitude","latitude")])[ ,2]

  cvars <- na.omit(moccsi[c("v1","v2")])
 
  mod.taxa.2d[[moi]] <- bkde2D(cvars,bandwidth=c(bw1,bw2),
    gridsize = c(101L, 101L), range.x = list(range_x, range_y))
  names(mod.taxa.2d[[moi]]) <-  c(nameX,nameY,"Prob")
  names(mod.taxa.2d)[moi] <- names(moccs)[moi]

  # partial distibution list
  nr <- nrow(cvars)
  samp50 <- list()
  
  for(j in 1:100){		# 100 k repetitions
    # subsample 50% of occurences per taxa
    s50 <- sample(1:nr, ceiling(nr * 0.50))
    varx <- cvars[s50,]

    # kernel 2D distribution
    samp50[[j]] <- bkde2D(varx , bandwidth = c(bw1,bw2),
      gridsize = c(101L, 101L), range.x = list(range_x, range_y))$fhat
  } # j

  # Average of subsamples and replace distribution
  mod.taxa.2d[[moi]][[3]] <- Reduce("+", samp50)/100
} # i

# print average distribution per taxa
pdf("Modern_2D_2025feb7.pdf") #### modern pdf ####### change date
  for(moi in 1:length(mod.taxa.2d)){
    image(mod.taxa.2d[[moi]][[1]],mod.taxa.2d[[moi]][[2]],mod.taxa.2d[[moi]][[3]],
	main= names(mod.taxa.2d)[moi],
   	xlab = nameX, ylab = nameY, col = viridis(50))
    box()
  }
dev.off()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# FOSSIL SAMPLES	# DEBUG (ed)

# Ensure fos_affs is correctly aggregated
fos_affs <- aggregate(affs$affinity ~ affs$morphotype, FUN = c)

fossil.affs <- list()

for (fi in 1:nrow(fossil)) {
  # fossil species presence
  fossp <- names(fossil[fi, which(fossil[fi, ] > 0)])
  fossp.dist <- list()
  
  print(paste("Processing fossil row:", fi))
  print(paste("Fossil species present:", paste(fossp, collapse = ", ")))
  
  for (fj in 1:length(fossp)) {
    # extract modern affinities for each fossil taxa
    modaff <- unlist(fos_affs[which(fos_affs$`affs$morphotype` == fossp[fj]), 2])
    
    print(paste("Processing fossil species:", fossp[fj]))
    print(paste("Modern affinities:", paste(modaff, collapse = ", ")))
    
    # locate in the modern 2d distributions 
    mods <- mod.taxa.2d[which(names(mod.taxa.2d) %in% modaff)]
    mods.z <- lapply(mods, function (ff) ff[[3]])
    
    if (length(mods.z) > 0) {
      # combine records if they are different taxa
      fosspjj <- Reduce("+", mods.z)  
      
      # make the distribution sum = 1
      fossp.dist[[fj]] <- fosspjj / sum(fosspjj)
    } else {
      fossp.dist[[fj]] <- matrix(0, nrow = nrow(mod.taxa.2d[[1]][[3]]), ncol = ncol(mod.taxa.2d[[1]][[3]]))
    }
  }
  
  if (length(fossp.dist) > 0) {
    fossil.affs[[fi]] <- Reduce("+", fossp.dist)
  } else {
    fossil.affs[[fi]] <- matrix(0, nrow = nrow(mod.taxa.2d[[1]][[3]]), ncol = ncol(mod.taxa.2d[[1]][[3]]))
  }
}

# Print the final fossil.affs for verification
print("Final fossil.affs:")
print(fossil.affs)

# PLOT
axisx <- mod.taxa.2d[[1]][[1]]
axisy <- mod.taxa.2d[[1]][[2]]

pdf("FOSSIL_2D_MATMAP_2025feb7.pdf") ####### change date
for (foi in 1:length(fossil.affs)) {
  image(axisx, axisy, fossil.affs[[foi]],
        main = rownames(fossil)[foi],
        xlab = nameX, ylab = nameY, col = viridis(50))
  box()
}
dev.off()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # Estimates

  ci <- c(0.025,0.975)
  quant <- 0.95

 estimates <- matrix(nrow = nrow(fossil),ncol = 10)
 colnames(estimates) <- c("Mean var1","Median var1","SD var1","LL var1","UL var1",
                "Mean var2","Median var2","SD var2","LL var2","UL var2")
 rownames(estimates ) <- rownames(fossil)

  for(i in 1:nrow(estimates)){
     stims <- which(fossil.affs[[i]] > quantile(fossil.affs[[i]],quant),arr.ind=TRUE)
     stims.x <- axisx[stims[,1]]
     stims.y <- axisy[stims[,2]]

     estimates[i,1] <- mean(stims.x)
     estimates[i,2] <- quantile(stims.x,0.5)
     estimates[i,3] <- sd(stims.x)
     estimates[i,4] <- quantile(stims.x,ci[1])
     estimates[i,5] <- quantile(stims.x,ci[2])

     estimates[i,6] <- mean(stims.y)
     estimates[i,7] <- quantile(stims.y,0.5)
     estimates[i,8] <- sd(stims.y)
     estimates[i,9] <- quantile(stims.y,ci[1])
     estimates[i,10] <- quantile(stims.y,ci[2])
} 
estimates
write.csv(estimates,"Estimates_2025feb7_MAT_MAP.csv") ####### change date 
estimates1_TPA<-estimates

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# ahora para 
#BIO17 = Precipitation of Driest Quarter
#BIO16 = Precipitation of Wettest Quarter


Var1 <- wc.bio[[17]] 	# PDQ 0-1543
Var2 <- wc.bio[[16]] 	# PWQ 0-5757

nameX <- "PDQ"
nameY <- "PWQ"

range_x <- c(0,1000)
range_y <- c(0,1500)   ### PP restricted to  see the lower limit

 # select a bandwidth for  kernel distributions
bw1 <- abs(max(range_x)- min(range_x)) /20
bw2 <- abs(max(range_y)- min(range_y)) /20

###########################################################
# Afinidad climática por taxon moderno

mod.taxa.2d <- list()

for(moi in 1:length(moccs)){

  # extract complete climatic data for spatial occurrences
  moccsi <- moccs[[moi]][c("longitude","latitude")]
  moccsi$v1 <- extract(Var1, moccsi[,c("longitude","latitude")])[ ,2]
  moccsi$v2 <- extract(Var2, moccsi[,c("longitude","latitude")])[ ,2]

  cvars <- na.omit(moccsi[c("v1","v2")])
 
  mod.taxa.2d[[moi]] <- bkde2D(cvars,bandwidth=c(bw1,bw2),
    gridsize = c(101L, 101L), range.x = list(range_x, range_y))
  names(mod.taxa.2d[[moi]]) <-  c(nameX,nameY,"Prob")
  names(mod.taxa.2d)[moi] <- names(moccs)[moi]

  # partial distibution list
  nr <- nrow(cvars)
  samp50 <- list()
  
  for(j in 1:100){		# 100 k repetitions
    # subsample 50% of occurences per taxa
    s50 <- sample(1:nr, ceiling(nr * 0.50))
    varx <- cvars[s50,]

    # kernel 2D distribution
    samp50[[j]] <- bkde2D(varx , bandwidth = c(bw1,bw2),
      gridsize = c(101L, 101L), range.x = list(range_x, range_y))$fhat
  } # j

  # Average of subsamples and replace distribution
  mod.taxa.2d[[moi]][[3]] <- Reduce("+", samp50)/100
} # i

# print average distribution per taxa
pdf("Modern_2D_2025feb7_PDQ_PWQ.pdf") #### modern pdf   ####### change date
  for(moi in 1:length(mod.taxa.2d)){
    image(mod.taxa.2d[[moi]][[1]],mod.taxa.2d[[moi]][[2]],mod.taxa.2d[[moi]][[3]],
	main= names(mod.taxa.2d)[moi],
   	xlab = nameX, ylab = nameY, col = viridis(50))
    box()
  }
dev.off()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# FOSSIL SAMPLES	# DEBUG (ed)

# Ensure fos_affs is correctly aggregated
fos_affs <- aggregate(affs$affinity ~ affs$morphotype, FUN = c)

fossil.affs <- list()

for (fi in 1:nrow(fossil)) {
  # fossil species presence
  fossp <- names(fossil[fi, which(fossil[fi, ] > 0)])
  fossp.dist <- list()
  
  print(paste("Processing fossil row:", fi))
  print(paste("Fossil species present:", paste(fossp, collapse = ", ")))
  
  for (fj in 1:length(fossp)) {
    # extract modern affinities for each fossil taxa
    modaff <- unlist(fos_affs[which(fos_affs$`affs$morphotype` == fossp[fj]), 2])
    
    print(paste("Processing fossil species:", fossp[fj]))
    print(paste("Modern affinities:", paste(modaff, collapse = ", ")))
    
    # locate in the modern 2d distributions 
    mods <- mod.taxa.2d[which(names(mod.taxa.2d) %in% modaff)]
    mods.z <- lapply(mods, function (ff) ff[[3]])
    
    if (length(mods.z) > 0) {
      # combine records if they are different taxa
      fosspjj <- Reduce("+", mods.z)  
      
      # make the distribution sum = 1
      fossp.dist[[fj]] <- fosspjj / sum(fosspjj)
    } else {
      fossp.dist[[fj]] <- matrix(0, nrow = nrow(mod.taxa.2d[[1]][[3]]), ncol = ncol(mod.taxa.2d[[1]][[3]]))
    }
  }
  
  if (length(fossp.dist) > 0) {
    fossil.affs[[fi]] <- Reduce("+", fossp.dist)
  } else {
    fossil.affs[[fi]] <- matrix(0, nrow = nrow(mod.taxa.2d[[1]][[3]]), ncol = ncol(mod.taxa.2d[[1]][[3]]))
  }
}

# Print the final fossil.affs for verification
print("Final fossil.affs:")
# print(fossil.affs)

# PLOT
axisx <- mod.taxa.2d[[1]][[1]]
axisy <- mod.taxa.2d[[1]][[2]]

pdf("FOSSIL_2D_PDQ_PWQ_2025feb7.pdf")       ####### change date
for (foi in 1:length(fossil.affs)) {
  image(axisx, axisy, fossil.affs[[foi]],
        main = rownames(fossil)[foi],
        xlab = nameX, ylab = nameY, col = viridis(50))
  box()
}
dev.off()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # Estimates

  ci <- c(0.025,0.975)
  quant <- 0.95

 estimates <- matrix(nrow = nrow(fossil),ncol = 10)
 colnames(estimates) <- c("Mean var1","Median var1","SD var1","LL var1","UL var1",
                "Mean var2","Median var2","SD var2","LL var2","UL var2")
 rownames(estimates ) <- rownames(fossil)

  for(i in 1:nrow(estimates)){
     stims <- which(fossil.affs[[i]] > quantile(fossil.affs[[i]],quant),arr.ind=TRUE)
     stims.x <- axisx[stims[,1]]
     stims.y <- axisy[stims[,2]]

     estimates[i,1] <- mean(stims.x)
     estimates[i,2] <- quantile(stims.x,0.5)
     estimates[i,3] <- sd(stims.x)
     estimates[i,4] <- quantile(stims.x,ci[1])
     estimates[i,5] <- quantile(stims.x,ci[2])

     estimates[i,6] <- mean(stims.y)
     estimates[i,7] <- quantile(stims.y,0.5)
     estimates[i,8] <- sd(stims.y)
     estimates[i,9] <- quantile(stims.y,ci[1])
     estimates[i,10] <- quantile(stims.y,ci[2])
} 
estimates
write.csv(estimates,"Estimates_2025feb7_PDQ_PWQ.csv")  ####### change date
estimates1_WDQ<-estimates

#######

estimates2_TPA <- read.csv("v estimates 2024_MAT_MAP.csv") ## data from first analyses without assessing uncertainty  for sampling or spatial distribution
colnames(estimates2_TPA) <- gsub("X", "", colnames(estimates2_TPA))  # Remove 'X'
colnames(estimates2_TPA) <- gsub("\\.", " ", colnames(estimates2_TPA))  # Replace '.' with a space
print(colnames(estimates2_TPA))

estimates2_WDQ <- read.csv("v estimates 2024_WetDryQ.csv") ## data from first analyses without assessing uncertainty  for sampling or spatial distribution
colnames(estimates2_WDQ) <- gsub("X", "", colnames(estimates2_WDQ))  # Remove 'X'
colnames(estimates2_WDQ) <- gsub("\\.", " ", colnames(estimates2_WDQ))  # Replace '.' with a space
print(colnames(estimates2_WDQ))


################################################
###########                                 MAT_MAP
################################################
########    MIOCENE Comparisons  (rows 1:9)   MAT_MAP
################################################

## MIOCENE MAT
group1_m1_1 <- estimates1_TPA[1:9, "Mean var1"]  ## after  subsampling_2025
group2_m1_1 <- estimates2_TPA[1:9, "Mean var1"] ## prior  subsampling_2024


# MAT Miocene Density Plot  
d1_m1 <- density(group1_m1_1)
d2_m1 <- density(group2_m1_1)
plot(d1_m1, col = "red", lwd = 2,
     main = "Density Plot:  MAT (Miocene data)",
     xlab = "MAT (ºC)")
lines(d2_m1, col = "blue", lwd = 2)
legend("topright", legend = c("after  subsampling_2025", "initial data 2024"),
       col = c("red", "blue"), lwd = 2)

# MAT Miocene Boxplot 
data_m1_1 <- data.frame(
  value = c(group1_m1_1, group2_m1_1),
  Source = rep(c("after  subsampling_2025", "initial data 2024"), each = length(group1_m1_1))
)
boxplot(value ~ Source, data = data_m1_1,
        main = "Boxplot: MAT (Miocene data)",
        ylab = "MAT (ºC)")
stripchart(value ~ Source, data = data_m1_1,
           vertical = TRUE, method = "jitter",
           pch = 20, col = "darkgray", add = TRUE)


# MAT Miocene t-Test   
#  No statistically significant difference in means between the two groups    p-value = 0.2468641
ttest_m1_1 <- t.test(group1_m1_1, group2_m1_1)
print(ttest_m1_1)
p_value1 <- ttest_m1_1$p.value

if (p_value1 < 0.05) {
  cat("There is a statistically significant difference in means between the two groups (p-value =", p_value1, ").\n")
} else {
  cat("There is no statistically significant difference in means between the two groups (p-value =", p_value1, ").\n")
}


# MIOCENE MAP Comparison  (rows 1:9)
group1_m2_1 <- estimates1_TPA[1:9, "Mean var2"]  ## after  subsampling_2025
group2_m2_1 <- estimates2_TPA[1:9, "Mean var2"] ## prior  subsampling_2024

# MAP Miocene Density Plot  
d1_m2 <- density(group1_m2_1)
d2_m2 <- density(group2_m2_1)
plot(d1_m2, col = "red", lwd = 2, ylim=c(0,0.006),
     main = "Density Plot: MAP (Miocene data)",
     xlab = "MAP (mm)")
lines(d2_m2, col = "blue", lwd = 2)
legend("topright", legend = c("after  subsampling_2025", "initial data 2024"),
       col = c("red", "blue"), lwd = 2)

# MAP Miocene Boxplot 
data_m2_1 <- data.frame(
  value = c(group1_m2_1, group2_m2_1),
  Source = rep(c("after  subsampling_2025", "initial data 2024"), each = length(group1_m2_1))
)
boxplot(value ~ Source, data = data_m2_1,
        main = "Boxplot: MAP (Miocene data)",
        ylab = "MAP (mm)")
stripchart(value ~ Source, data = data_m2_1,
           vertical = TRUE, method = "jitter",
           pch = 20, col = "darkgray", add = TRUE)
           
           
# MAP Miocene t-Test   
#  No statistically significant difference in means between the two groups    p-value = 0.8040923
ttest_m2_1 <- t.test(group1_m2_1, group2_m2_1)
print(ttest_m2_1)
p_value2 <- ttest_m2_1 $p.value

if (p_value2 < 0.05) {
  cat("There is a statistically significant difference in means between the two groups (p-value =", p_value2,").\n")
} else {
  cat("There is no statistically significant difference in means between the two groups (p-value =", p_value2,").\n")
}
    
           
################################################
########    HOLOCENE Comparisons  (rows 10:18)
################################################      

## HOLOCENE MAT
group1_m1_2 <- estimates1_TPA[10:18, "Mean var1"] ## after  subsampling_2025
group2_m1_2 <- estimates2_TPA[10:18, "Mean var1"]  ## prior  subsampling_2024

# MAT Holocene Density Plot  
d1_m1_2 <- density(group1_m1_2)
d2_m1_2 <- density(group2_m1_2)
plot(d1_m1_2, col = "red", lwd = 2, ylim=c(0,0.13),
     main = "Density Plot: MAT (Holocene data)",
     xlab = "MAT (ºC)")
lines(d2_m1_2, col = "blue", lwd = 2)
legend("topright", legend = c("after  subsampling_2025", "initial data 2024"),
       col = c("red", "blue"), lwd = 2)

# MAT Miocene Boxplot 
data_m1_2 <- data.frame(
  value = c(group1_m1_2, group2_m1_2),
  Source = rep(c("after  subsampling_2025", "initial data 2024"), each = length(group1_m1_2))
)
boxplot(value ~ Source, data = data_m1_2,
        main = "Boxplot: MAT (Holocene data)",
        ylab = "MAT (ºC)")
stripchart(value ~ Source, data = data_m1_2,
           vertical = TRUE, method = "jitter",
           pch = 20, col = "darkgray", add = TRUE)

# MAT Holocene t-Test   
#  No statistically significant difference in means between the two groups    p-value = 0.7071261
ttest_m1_2 <- t.test(group1_m1_2, group2_m1_2)
print(ttest_m1_2)
p_value3 <- ttest_m1_2 $p.value

if (p_value3 < 0.05) {
  cat("There is a statistically significant difference in means between the two groups (p-value =", p_value3,").\n")
} else {
  cat("There is no statistically significant difference in means between the two groups (p-value =", p_value3,").\n")
}
      

# HOLOCENE MAP Comparison 
group1_m2_2 <- estimates1_TPA[10:18, "Mean var2"]  ## after  subsampling_2025
group2_m2_2 <- estimates2_TPA[10:18, "Mean var2"]   ## prior  subsampling_2024

# MAP Holocene Density Plot  
d1_m2_2 <- density(group1_m2_2)
d2_m2_2 <- density(group2_m2_2)
plot(d1_m2_2, col = "red", lwd = 2, ylim=c(0,0.007),
     main = "Density Plot: MAP (Holocene data)",
     xlab = "MAP (mm)")
lines(d2_m2_2, col = "blue", lwd = 2)
legend("topright", legend = c("after  subsampling_2025", "initial data 2024"),
       col = c("red", "blue"), lwd = 2)

# MAT Miocene Boxplot 
data_m2_2 <- data.frame(
  value = c(group1_m2_2, group2_m2_2),
  Source = rep(c("after subsampling_2025", "initial data 2024"), each = length(group1_m2_2))
)
boxplot(value ~ Source, data = data_m2_2,
        main = "Boxplot: MAP (Holocene data)",
        ylab = "MAP (mm)")
stripchart(value ~ Source, data = data_m2_2,
           vertical = TRUE, method = "jitter",
           pch = 20, col = "darkgray", add = TRUE)

# MAT Holocene t-Test   
#  No statistically significant difference in means between the two groups    p-value = 0.7813607
ttest_m2_2 <- t.test(group1_m2_2, group2_m2_2)
print(ttest_m2_2)
p_value4 <- ttest_m2_2 $p.value

if (p_value4 < 0.05) {
  cat("There is a statistically significant difference in means between the two groups (p-value =", p_value4,").\n")
} else {
  cat("There is no statistically significant difference in means between the two groups (p-value =", p_value4,").\n")
}

golden_yellow <- "#af8e02"
deep_blue <- "#131384"

######### Compiled graphs:  
par(mfrow = c(2, 2)) 

# Boxplot for Mean var1 (Rows 1-9)
boxplot(value ~ Source, data = data_m1_1,
        main = "Miocene Mean Annual Temperature - MAT (ºC)",
        ylab = "MAT", col = c(golden_yellow, golden_yellow))
stripchart(value ~ Source, data = data_m1_1,
           vertical = TRUE, method = "jitter",
           pch = 20, col = "darkgray", add = TRUE)

# Boxplot for Mean var2 (Rows 1-9)
boxplot(value ~ Source, data = data_m2_1,
        main = "Miocene Mean Annual Precipitation - MAP (mm)",
        ylab = "MAP", col = c(golden_yellow, golden_yellow))
stripchart(value ~ Source, data = data_m2_1,
           vertical = TRUE, method = "jitter",
           pch = 20, col = "darkgray", add = TRUE)

# Boxplot for Mean var1 (Rows 10-18)
boxplot(value ~ Source, data = data_m1_2,
        main = "Holocene Mean Annual Temperature - MAT (ºC)",
        ylab = "MAT", col = c(deep_blue, deep_blue))
stripchart(value ~ Source, data = data_m1_2,
           vertical = TRUE, method = "jitter",
           pch = 20, col = "darkgray", add = TRUE)

# Boxplot for Mean var2 (Rows 10-18)
boxplot(value ~ Source, data = data_m2_2,
        main = "Holocene Mean Annual Precipitation - MAP (mm)",
        ylab = "MAP", col = c(deep_blue, deep_blue))
stripchart(value ~ Source, data = data_m2_2,
           vertical = TRUE, method = "jitter",
           pch = 20, col = "darkgray", add = TRUE)


# Density plot
par(mfrow = c(2, 2))

# Density plot for Mean var1 (Rows 1-9)
plot(d1_m1, col = "red", lwd = 2,
     main = "Miocene Mean Annual Temperature - MAT (ºC)",
     xlab = "MAT")
lines(d2_m1, col = "blue", lwd = 2)
legend("topright", legend = c("after  subsampling_2025", "initial data 2024"),
       col = c("red", "blue"), lwd = 2)

# Density plot for Mean var2 (Rows 1-9)
plot(d1_m2, col = "red", lwd = 2,ylim=c(0,0.006),
     main = "Miocene Mean Annual Precipitation - MAP (mm)",
     xlab = "MAP")
lines(d2_m2, col = "blue", lwd = 2)
legend("topright", legend = c("after  subsampling_2025", "initial data 2024"),
       col = c("red", "blue"), lwd = 2)

# Density plot for Mean var1 (Rows 10-18)
plot(d1_m1_2, col = "red", lwd = 2, ylim=c(0,0.13),
     main = "Holocene Mean Annual Temperature - MAT (ºC)",
     xlab = "MAT")
lines(d2_m1_2, col = "blue", lwd = 2)
legend("topright", legend = c("after  subsampling_2025", "initial data 2024"),
       col = c("red", "blue"), lwd = 2)

# Density plot for Mean var2 (Rows 10-18)
plot(d1_m2_2, col = "red", lwd = 2, ylim=c(0,0.007),
     main = "Holocene Mean Annual Precipitation - MAP (mm)",
     xlab = "MAP")
lines(d2_m2_2, col = "blue", lwd = 2)
legend("topright", legend = c("after  subsampling_2025", "initial data 2024"),
       col = c("red", "blue"), lwd = 2)
       

################################################
###########                                 WetDryQ
################################################
########    MIOCENE Comparisons  (rows 1:9)   WetDryQ
################################################


## MIOCENE PDQ
group1_m1_1 <- estimates1_WDQ[1:9, "Mean var1"]  ## after  subsampling_2025  PDQ
group2_m1_1 <- estimates2_WDQ[1:9, "Mean var1"] ## prior  subsampling_2024 PDQ


# PDQ Miocene Density Plot  
d1_m1 <- density(group1_m1_1)
d2_m1 <- density(group2_m1_1)
plot(d1_m1, col = "red", lwd = 2, xlim=c(0,100), ylim=c(0,0.2),
     main = "Density Plot:  PDQ (Miocene data)",
     xlab = "PDQ (mm)")
lines(d2_m1, col = "blue", lwd = 2)
legend("topright", legend = c("after  subsampling_2025", "initial data 2024"),
       col = c("red", "blue"), lwd = 2)

# PDQ Miocene Boxplot 
data_m1_1 <- data.frame(
  value = c(group1_m1_1, group2_m1_1),
  Source = rep(c("after  subsampling_2025", "initial data 2024"), each = length(group1_m1_1))
)
boxplot(value ~ Source, data = data_m1_1,
        main = "Boxplot: PDQ (Miocene data)",
        ylab = "PDQ (mm)")
stripchart(value ~ Source, data = data_m1_1,
           vertical = TRUE, method = "jitter",
           pch = 20, col = "darkgray", add = TRUE)

# PDQ Miocene t-Test   
#  There is a statistically significant difference in means between the two groups    p-value =  2.079846e-08
ttest_m1_1 <- t.test(group1_m1_1, group2_m1_1)
print(ttest_m1_1)
p_value5 <- ttest_m1_1$p.value

if (p_value5 < 0.05) {
  cat("There is a statistically significant difference in means between the two groups (p-value =", p_value5, ").\n")
} else {
  cat("There is no statistically significant difference in means between the two groups (p-value =", p_value5, ").\n")
}



# MIOCENE PWQ Comparison  (rows 1:9)
group1_m2_1 <- estimates1_WDQ[1:9, "Mean var2"]  ## after  subsampling_2025 PWQ
group2_m2_1 <- estimates2_WDQ[1:9, "Mean var2"] ## prior  subsampling_2024 PWQ
 
# PWQ Miocene Density Plot  
d1_m2 <- density(group1_m2_1)
d2_m2 <- density(group2_m2_1)
plot(d1_m2, col = "red", lwd = 2, xlim=c(150,400), ylim=c(0,0.022),
     main = "Density Plot: PWQ (Miocene data)",
     xlab = "PWQ (mm)")
lines(d2_m2, col = "blue", lwd = 2)
legend("topright", legend = c("after  subsampling_2025", "initial data 2024"),
       col = c("red", "blue"), lwd = 2)

# PWQ Miocene Boxplot 
data_m2_1 <- data.frame(
  value = c(group1_m2_1, group2_m2_1),
  Source = rep(c("after  subsampling_2025", "initial data 2024"), each = length(group1_m2_1))
)
boxplot(value ~ Source, data = data_m2_1,
        main = "Boxplot: PWQ (Miocene data)",
        ylab = "PWQ (mm)")
stripchart(value ~ Source, data = data_m2_1,
           vertical = TRUE, method = "jitter",
           pch = 20, col = "darkgray", add = TRUE)
           
           
# PWQ Miocene t-Test   
#  There is a statistically significant difference in means between the two groups      p-value = 5.018901e-05
ttest_m2_1 <- t.test(group1_m2_1, group2_m2_1)
print(ttest_m2_1)
p_value6 <- ttest_m2_1 $p.value

if (p_value6 < 0.05) {
  cat("There is a statistically significant difference in means between the two groups (p-value =", p_value6,").\n")
} else {
  cat("There is no statistically significant difference in means between the two groups (p-value =", p_value6,").\n")
}
    



           
################################################
########    HOLOCENE Comparisons  (rows 10:18)
################################################      

## HOLOCENE PDQ
group1_m1_2 <- estimates1_WDQ[10:18, "Mean var1"] ## after  subsampling_2025  PDQ
group2_m1_2 <- estimates2_WDQ[10:18, "Mean var1"]  ## prior  subsampling_2024  PDQ

# PDQ Holocene Density Plot  
d1_m1_2 <- density(group1_m1_2)
d2_m1_2 <- density(group2_m1_2)
plot(d1_m1_2, col = "red", lwd = 2, ylim=c(0,0.05), xlim=c(0,150),
     main = "Density Plot: PDQ (Holocene data)",
     xlab = "PDQ (mm)")
lines(d2_m1_2, col = "blue", lwd = 2)
legend("topright", legend = c("after  subsampling_2025", "initial data 2024"),
       col = c("red", "blue"), lwd = 2)

# PDQ Miocene Boxplot 
data_m1_2 <- data.frame(
  value = c(group1_m1_2, group2_m1_2),
  Source = rep(c("after  subsampling_2025", "initial data 2024"), each = length(group1_m1_2))
)
boxplot(value ~ Source, data = data_m1_2,
        main = "Boxplot: MAT (Holocene data)",
        ylab = "PDQ (mm)")
stripchart(value ~ Source, data = data_m1_2,
           vertical = TRUE, method = "jitter",
           pch = 20, col = "darkgray", add = TRUE)

# PDQ Holocene t-Test   
#  There is a statistically significant difference in means between the two groups       p-value  = 2.44097e-07
ttest_m1_2 <- t.test(group1_m1_2, group2_m1_2)
print(ttest_m1_2)
p_value7 <- ttest_m1_2 $p.value

if (p_value7 < 0.05) {
  cat("There is a statistically significant difference in means between the two groups (p-value =", p_value7,").\n")
} else {
  cat("There is no statistically significant difference in means between the two groups (p-value =", p_value7,").\n")
}



# HOLOCENE PWQ Comparison 
group1_m2_2 <- estimates1_WDQ[10:18, "Mean var2"] #  after  subsampling_2025 PWQ
group2_m2_2 <- estimates2_WDQ[10:18, "Mean var2"] #  prior  subsampling_2024  PWQ

# PWQ Holocene Density Plot  
d1_m2_2 <- density(group1_m2_2)
d2_m2_2 <- density(group2_m2_2)
plot(d1_m2_2, col = "red", lwd = 2, ylim=c(0,0.035), xlim=c(150,350),
     main = "Density Plot: PWQ (Holocene data)",
     xlab = "PWQ (mm)")
lines(d2_m2_2, col = "blue", lwd = 2)
legend("topright", legend = c("after  subsampling_2025", "initial data 2024"),
       col = c("red", "blue"), lwd = 2)

# PWQ Miocene Boxplot 
data_m2_2 <- data.frame(
  value = c(group1_m2_2, group2_m2_2),
  Source = rep(c("after subsampling_2025", "initial data 2024"), each = length(group1_m2_2))
)
boxplot(value ~ Source, data = data_m2_2,
        main = "Boxplot: PWQ (Holocene data)",
        ylab = "PWQ (mm)")
stripchart(value ~ Source, data = data_m2_2,
           vertical = TRUE, method = "jitter",
           pch = 20, col = "darkgray", add = TRUE)

# PWQ Holocene t-Test   
#  There is a statistically significant difference in means between the two groups        p-value = 1.670438e-09 
ttest_m2_2 <- t.test(group1_m2_2, group2_m2_2)
print(ttest_m2_2)
p_value8 <- ttest_m2_2 $p.value

if (p_value8 < 0.05) {
  cat("There is a statistically significant difference in means between the two groups (p-value =", p_value8,").\n")
} else {
  cat("There is no statistically significant difference in means between the two groups (p-value =", p_value8,").\n")
}



############ Compiled graphs:  
par(mfrow = c(2, 2)) 

# Boxplot for Mean var1 (Rows 1-9)
boxplot(value ~ Source, data = data_m1_1,
        main = "Miocene Precipitation Dryest Quarter - PDQ (mm)",
        ylab = "PDQ", col = c(golden_yellow, golden_yellow))
stripchart(value ~ Source, data = data_m1_1,
           vertical = TRUE, method = "jitter",
           pch = 20, col = "darkgray", add = TRUE)

# Boxplot for Mean var2 (Rows 1-9)
boxplot(value ~ Source, data = data_m2_1,
        main = "Miocene Precipitation Wettest Quarter - PWQ (mm)",
        ylab = "PWQ", col = c(golden_yellow, golden_yellow))
stripchart(value ~ Source, data = data_m2_1,
           vertical = TRUE, method = "jitter",
           pch = 20, col = "darkgray", add = TRUE)

# Boxplot for Mean var1 (Rows 10-18)
boxplot(value ~ Source, data = data_m1_2,
        main = "Holocene Precipitation Dryest Quarter - PDQ (mm)",
        ylab = "PDQ", col = c(deep_blue, deep_blue))
stripchart(value ~ Source, data = data_m1_2,
           vertical = TRUE, method = "jitter",
           pch = 20, col = "darkgray", add = TRUE)

# Boxplot for Mean var2 (Rows 10-18)
boxplot(value ~ Source, data = data_m2_2,
        main = "Holocene Precipitation Wettest Quarter - PWQ (mm)",
        ylab = "PWQ", col = c(deep_blue, deep_blue))
stripchart(value ~ Source, data = data_m2_2,
           vertical = TRUE, method = "jitter",
           pch = 20, col = "darkgray", add = TRUE)


# Density plot
par(mfrow = c(2, 2))

# Density plot for Mean var1 (Rows 1-9)
plot(d1_m1, col = "red", lwd = 2, xlim=c(0,150), ylim=c(0,0.18),
     main = "Miocene Precipitation Dryest Quarter- PDQ (mm)",
     xlab = "PDQ")
lines(d2_m1, col = "blue", lwd = 2)
legend("topright", legend = c("after  subsampling_2025", "initial data 2024"),
       col = c("red", "blue"), lwd = 2)

# Density plot for Mean var2 (Rows 1-9)
plot(d1_m2, col = "red", lwd = 2, xlim=c(150,400), ylim=c(0,0.022),
     main = "Miocene Precipitation Wettest Quarter- PWQ (mm)",
     xlab = "PWQ")
lines(d2_m2, col = "blue", lwd = 2)
legend("topright", legend = c("after  subsampling_2025", "initial data 2024"),
       col = c("red", "blue"), lwd = 2)

# Density plot for Mean var1 (Rows 10-18)
plot(d1_m1_2, col = "red", lwd = 2, ylim=c(0,0.05), xlim=c(0,150),
     main = "Holocene Precipitation Dryest Quarter- PDQ (mm)",
     xlab = "PDQ")
lines(d2_m1_2, col = "blue", lwd = 2)
legend("topright", legend = c("after  subsampling_2025", "initial data 2024"),
       col = c("red", "blue"), lwd = 2)

# Density plot for Mean var2 (Rows 10-18)
plot(d1_m2_2, col = "red", lwd = 2, ylim=c(0,0.035), xlim=c(150,400),
     main = "Holocene Precipitation Wettest Quarter- PWQ (mm)",
     xlab = "PWQ")
lines(d2_m2_2, col = "blue", lwd = 2)
legend("topright", legend = c("after  subsampling_2025", "initial data 2024"),
       col = c("red", "blue"), lwd = 2)




################################################
###########         Differences between Holocene and Miocene  after Resampling  : MAP, PDQ, and PWQ
################################################
########    MAP  
################################################




