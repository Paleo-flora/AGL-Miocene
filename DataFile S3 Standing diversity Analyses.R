###
## Floristic Analysis: Diversity and Composition   -   Peruvian Miocene paleoflora  [Ochoa et al., xx]

library(fossil) 
library(ggplot2)
library(vegan)
library(readxl)
library(ggplot2)
library(dplyr)
library(reshape2)
library(gridExtra)
library(fossil) 
library(BIEN)
library(raster)
library(KernSmooth)
library(readxl)
library(CommEcol)
library(cluster)
library(rioja)

setwd("/Users/dianaochoa/Desktop/data for R/HP3 Paleoflora dic2023")

######################################### 
#DIVERSITY - DATA  WITH  ALL COUNTS 
######################################### 

# counts data
counts <- read_excel("AGL Coex counts_2023.xlsx" ,col_names = FALSE) ## All counts [species*samples], excluding pollenites/spore indet

pf.counts<-as.data.frame(counts[3:nrow(counts),-1]) ## Pollen counts start in row 3 - file without Dinocyst and Fungi
dim(pf.counts)###252 18
pf.counts[is.na(pf.counts)] <- 0 # NA's = zero
colnames(pf.counts)<- as.character(counts[2,2:ncol(counts)])
#rownames(pf.counts)<- as.character(counts[3:nrow(counts),1])
pf.counts <- as.data.frame(lapply(pf.counts,function (x) as.numeric(as.character(x))))
total_col<-apply(pf.counts, 2, sum) #sum per sample

sites <- ncol(pf.counts) #determine number of sites
species<-nrow(pf.counts) #determine number of species
pctgs <- data.frame(matrix(ncol = sites, nrow = species))

gr <- c(rep("Miocene",9),rep("Holocene",9)) ###Changes in proportions calculated using complete database 
Mio.sites<-c(1:9) ## depending in gr
Modern.sites<-c(10:18) ## depending in gr

all.sites<-c(1:ncol(pf.counts))
age.type <-(all.sites>9)*1 ##### ### 1 is modern and 0 is fossil-AGL
age.type<-as.vector(age.type) ####  1 is modern and 0 is fossil-AGL

abundancedata <- pf.counts 

rar.50<-rarefy(t(abundancedata),50,se= TRUE)
rar.50test=rar.50[1,which(rar.50[2,]>0)]
G1.DUC.50<-all.sites[which(rar.50[2,]>0)]
age.type50<-age.type[which(rar.50[2,]>0)]

rar.100<-rarefy(t(abundancedata),100,se= TRUE)
rar.100test=rar.100[1,which(rar.100[2,]>0)]
G1.DUC.100<-all.sites[which(rar.100[2,]>0)]
age.type100<-age.type[which(rar.100[2,]>0)]

rar.150<-rarefy(t(abundancedata),150,se= TRUE)
rar.150test=rar.150[1,which(rar.150[2,]>0)]
G1.DUC.150<-all.sites[which(rar.150[2,]>0)]
age.type150<-age.type[which(rar.150[2,]>0)]

rar.200 <-rarefy(t(abundancedata),200,se= TRUE)
rar.200test=rar.200[1,which(rar.200[2,]>0)]
G1.DUC.200<-all.sites[which(rar.200[2,]>0)]
age.type200<-age.type[which(rar.200[2,]>0)]


######Tproof:
## rar100
test.rar100<-t.test(rar.100test ~ age.type100)
test.rar100

## Welch Two Sample t-test
## data:  rar.100test by age.type100
## t = 1.8993, df = 13.842, p-value = 0.07856
## alternative hypothesis: true difference in means between group 0 and group 1 is not equal to 0
## 95 percent confidence interval:
## -0.5900276      9.6367403
## sample estimates:
## mean in group 0 mean in group 1 
##       28.11041        23.58705 


## rar150
test.rar150<-t.test(rar.150test~age.type150)
test.rar150

## Welch Two Sample t-test
## data:  rar.150test by age.type150
## t = 3.8831, df = 9.2414, p-value = 0.003531
## alternative hypothesis: true difference in means between group 0 and group 1 is not equal to 0
## 95 percent confidence interval:
##   4.594483 17.296967
## sample estimates:
## mean in group 0 mean in group 1 
##        35.32526        24.37953 

## rar200
test.rar200<-t.test(rar.200test ~ age.type200)
test.rar200

## Welch Two Sample t-test
## data:  rar.200test by age.type200
## t = 3.6743, df = 6.5658, p-value = 0.008866
## alternative hypothesis: true difference in means between group 0 and group 1 is not equal to 0
## 95 percent confidence interval:
##   5.044152    23.969629
## sample estimates:
## mean in group 0 mean in group 1 
##      42.81658        28.30969 



plot(G1.DUC.100, rar.100test, xlab="Miocene (sites 1 to 8) versus Holocene (sites 9 to 15)", ylab="Number of Species", pch=19, main="Rarefied richness", ylim=c(13,33))
points(G1.DUC.50, rar.50test, col="blue", pch=21)
legend("topright", legend= c("Rarefaction at 100","Rarefaction at 50"), pch=c(19,21), col=c("black", "blue"),  bg="white")


par(mfrow=c(1,3))
#boxplot(rar.50test ~ age.type50 ,col = "lightgray", 
#main = "Richness estimation Rar.50", ylim=c(15,45), ylab = "Estimated Richness",
#xlab = "Miocene (0) versus Holocene (1)")

boxplot(rar.100test ~ age.type100, col = "lightgray", 
main = "Richness estimation Rar.100",ylim=c(15,45),
xlab = "Miocene (0) versus Holocene (1)")

boxplot(rar.150test ~ age.type150, col = "lightgray", 
main = "Richness estimation Rar.150", ylim=c(15,45),
xlab = "Miocene (0) versus Holocene (1)")

boxplot(rar.200test ~ age.type200, col = "lightgray", 
main = "Richness estimation Rar.200", ylim=c(20,60),
xlab = "Miocene (0) versus Holocene (1)")


#######################################################
## Rarefaction Boxplots   - Within-sample diversity  ## All counts [species*samples] 
##################################################################

# Define USGS colors 
colors <- c("gold", "#FCD5B5") # Create a vector of colors for the Tortonian and Holocene

# Convert the data into a dataframe
rar.1<-cbind(rar.100test,age.type100)
rar.1_df <- data.frame(rar.100test = rar.1[,1],
  age.type100 = factor(rar.1[,2], levels = c(0,1), labels = c("Miocene", "Holocene")))
  
rar.2<-cbind(rar.150test,age.type150)
rar.2_df <- data.frame(rar.150test = rar.2[,1],
  age.type150 = factor(rar.2[,2], levels = c(0,1), labels = c("Miocene", "Holocene")))
  
rar.3<-cbind(rar.200test,age.type200)
rar.3_df <- data.frame(rar.200test = rar.3[,1],
  age.type200 = factor(rar.3[,2], levels = c(0,1), labels = c("Miocene", "Holocene")))


# Create the plot
plot.rar1 <-ggplot(rar.1_df, aes(x = age.type100, y = rar.100test)) +
  geom_boxplot(alpha = 0.5, size = 0.5,outlier.shape = NA, aes(fill = age.type100)) + 
  scale_fill_manual(values = colors)+
  geom_jitter(alpha = 0.5, size = 2, width = 0.2) +
  labs(x = "Age", y = "Rarefied diversity") +
  theme_classic()

plot.rar2 <-ggplot(rar.2_df, aes(x = age.type150, y = rar.150test))  +
  geom_boxplot(alpha = 0.5, size = 0.5, outlier.shape = NA,  aes(fill = age.type150)) +  
  scale_fill_manual(values = colors)+
  geom_jitter(alpha = 0.5, size = 2, width = 0.2) +
  labs(x = "Age", y = "Rarefied diversity") +
  theme_classic()

plot.rar3<-ggplot(rar.3_df, aes(x = age.type200, y = rar.200test))+
  geom_boxplot(alpha = 0.5, size = 0.5, outlier.shape = NA, aes(fill = age.type200)) +  
  scale_fill_manual(values = colors)+
  geom_jitter(alpha = 0.5, size = 2, width = 0.2) +
  labs(x = "Age", y = "Rarefied diversity") +
  theme_classic()


# Add titles and remove axis labels
plot.rar1 <- plot.rar1 + ggtitle("Rarefied diversity at 100") + ylab(NULL) + xlab(NULL)
plot.rar2 <- plot.rar2 + ggtitle("Rarefied diversity at 150") + ylab(NULL) + xlab(NULL)
plot.rar3 <- plot.rar3 + ggtitle("Rarefied diversity at 200") + ylab(NULL) + xlab(NULL)

# Remove legend
plot.rar1 <- plot.rar1 + guides(fill = FALSE)
plot.rar2 <- plot.rar2 + guides(fill = FALSE)
plot.rar3 <- plot.rar3 + guides(fill = FALSE)

# Combine plots, remove grid lines in strips and save
Rarefaction.plot<-grid.arrange(plot.rar1, plot.rar2, plot.rar3,  nrow = 1, top = "Rarefied diversity")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank())

#ggsave("Rarefied diversity.pdf", plot= Rarefaction.plot, dpi = 500)

  
##########################################################################
## Accumulation curve ## All counts [species*samples] 
#############################################################################
 pf.counts.acc<-t(pf.counts)
 Mio.counts <-(pf.counts.acc[1:9,])
 dim(Mio.counts)
 Hol.counts <-(pf.counts.acc[10:nrow(pf.counts.acc),])
 dim(Hol.counts)
acc_cur.Mio <- specaccum(Mio.counts)
sp2.M <- specaccum(Mio.counts, "random")
acc_cur.Hol <- specaccum(Hol.counts)
sp2.H <- specaccum(Hol.counts, "random")
par(mfrow=c(1,2))
plot(acc_cur.Mio, ci.type="poly", col="gold", main =("Miocene accumulation curve"), ylim=c(0,200), xlim=c(0,10), lwd=3, ci.lty=0, ci.col="light yellow")
boxplot(sp2.M, col="light gray", add=TRUE, pch="+")
plot(acc_cur.Hol,  ci.type="poly", col="#FCD5B5", main =("Holocene accumulation curve"), ylim=c(0,200),xlim=c(0,10), lwd=3, ci.lty=0, ci.col="mistyrose")
boxplot(sp2.H, col="light gray", add=TRUE, pch="+")

#####################################################################
## Standing diversity (after range-through assumption) ## All counts [species*samples] 
######################################################################

counts <- read_excel("AGL Coex counts_2023.xlsx" ,col_names = FALSE) ## All counts [species*samples], excluding pollenites/spore indet

conteo.all <-as.data.frame(t(counts[3:nrow(counts),-1]))## Pollen counts start in row 3 - file without Dinocyst and Fungi
dim(conteo.all)###252 18
conteo.all[is.na(conteo.all )] <- 0 # NA's = zero
colnames(conteo.all)<- counts[-(1:2),1]
rownames(conteo.all)<- counts[2,-1]

fossil <- as.data.frame(lapply(conteo.all,function (x) as.numeric(as.character(x))))
colnames(fossil)<- counts[-(1:2),1]
rownames(fossil)<- counts[2,-1]
# NA's = zero
fossil[is.na(fossil)] <- 0

fossil.t<-t(fossil)##species in rows

Miocene.data<-fossil.t[, 1:9]
Holocene.data<-fossil.t[, 10:ncol(fossil.t)]
dim(Miocene.data) ##should be #sps=80 by 9 samples
dim(Holocene.data) #should be #sps=80 by 9 samples
samples<-c(1:9) #same for Holocene and Miocene, as similar number of samples


##Rarefaction Miocene vs Holocene 
########## code for range-through analysis
fill.occur=function(sp)
{
occur=which(sp>0)  ##array of row numbers where sp>0
fad=occur[1] #row of first number in array
numboccur=length(occur)# how many numbers in the array
lad=occur[numboccur]#row number of the last position in the array
alloccur=rep(0,length(sp))## produces a matrix of zeroes that is similar in size to original
i=1:length(sp)## all positions in the matrix are now labeled i
alloccur[i>=fad & i<=lad]=1## replaces all i in between fad and lad by 1
return(alloccur)#gives the matrix out
}


Mio.Rthrough<-t(apply(Miocene.data,1,fill.occur))## matrix after range through, samples in columns
Mio.total<-apply(Mio.Rthrough,2,sum)##diversity sum

Hol.Rthrough <-t(apply(Holocene.data,1,fill.occur))## matrix after range through, samples in columns
Hol.total<-apply(Hol.Rthrough,2,sum)##diversity sum

leg.txt<-c("Miocene", "Holocene")
plot(samples,  Mio.total, main="Range-through standing diversity", col="dark blue", pch=20, xlab="observed diverstity", ylab="Samples")
points(samples, Hol.total, col="brown", pch=20)
legend("topright", leg.txt, pch =c(20), col = c("dark blue", "brown"), cex=0.7)


## DROPING OUT SINGLETONES 

Mio.summa1<-apply(Miocene.data,1,sum)
Hol.summa1<-apply(Holocene.data,1,sum)

Mio.uniq<-Miocene.data[which(Mio.summa1>1),] #44 sp
Hol.uniq<-Holocene.data[which(Hol.summa1>1),] #24 sp

MioRT.uniq<-t(apply(Mio.uniq,1,fill.occur))## matrix after range through, samples in columns
Mio.unique<-apply(MioRT.uniq,2,sum)##diversity sum data without singles

HolRT.uniq <-t(apply(Hol.uniq,1,fill.occur))## matrix after range through, samples in columns
Holo.unique<-apply(HolRT.uniq,2,sum)##diversity sum data without singles

leg.txt2<-c("Miocene", "Holocene", "Miocene No singletones", "Holocene No singletones")
plot(samples,  Mio.total, main="Range-through diversity", col="dark blue", pch=20, xlab="observed diverstity", ylab="Samples", ylim=c(15, 95))
points(samples, Hol.total, col="brown", pch=20)
points(samples,  Mio.unique, col="dark blue", pch=18, )
points(samples, Holo.unique, col="brown", pch=18)
legend("topright", leg.txt2, pch =c(20, 20, 18, 18), col = c("dark blue", "brown", "dark blue", "brown"), cex=0.7)


#####################################################################
## Simpson index ## All counts [species*samples] using  samples with counts above 80 grains
######################################################################
age.type <- c(1:ncol(fossil.t))
age.type <-(age.type>9)*1 ##### ### 1 is modern and 0 is fossil-AGL
age.type<-as.vector(age.type) ####  1 is modern and 0 is fossil-AGL

species_data <- fossil[rowSums(fossil) > 80, ]
age.type80<-age.type[rowSums(fossil) > 80]
simpson_index <- diversity(species_data, index = "simpson")

# Print Simpson index
print(simpson_index)

Simpson.test<-t.test(simpson_index ~ age.type80)
Simpson.test

#	Welch Two Sample t-test
# data:  simpson_index by age.type80
# t = 0.94116, df = 11.291, p-value = 0.3663
# alternative hypothesis: true difference in means between group 0 and group 1 is not equal to 0
# 95 percent confidence interval:
# -0.01252323      0.03133740
# sample estimates:
# mean in group 0 mean in group 1 
 #      0.9018213       0.8924142 

#p_value < 0.05; therefore there are no significant differences within groups.

#####################################################################
##  Detrended Correspondence Analyses (DCA) ## All counts [species*samples]
######################################################################
abundancedata <- fossil.t 
decorana(t(abundancedata))->deco

plot(deco)
plot(deco, display="sites")####samples
plot(deco, display="species")####species
plot(deco, display="both")####to see both.. samples and species
deco
# Call: decorana(veg = t(abundancedata)) 
# Detrended correspondence analysis with 26 segments.
# Rescaling of axes with 4 iterations.
# Total inertia (scaled Chi-square): 3.0744

#                      				DCA1   		DCA2   		DCA3    		DCA4
#  Eigenvalues          	0.6311 		0.2628 		0.1291 	0.13845
#  Additive Eigenvalues 0.6311 	0.2453 		0.1197 	0.13134
#  Decorana values      0.6327 		0.2856 		0.1298 	0.06404
#  Axis lengths         		2.9694 		2.3861 		1.4768 	1.56032

scores(deco)->values
DCA1<-values[,1] ## DCA1 0.6311
DCA2<-values[,2] ## DCA2 0.2628
all.samples<-c(1:18)
## DCA versus samples
plot(DCA1[1:9], all.samples[1:9], main="DCA Miocene versus Holocene", col="dark blue", pch=20, xlim=c(-2,2), ylim=c(0,18), xlab="DCA  axis 1", ylab="Samples") ## Miocene samples
points(DCA1[10:18], all.samples[10:18], main="DCA Miocene versus Holocene", col="brown", pch=18) ## Holocene samples
legend("topleft", legend= c("Miocene","Holocene"), pch=c(20,18), col=c("dark blue", "brown"),  bg="white")

values
  #                   		  		DCA1       		 DCA2        			DCA3        			DCA4
# AGL1          		  -1.12184918 		-0.42833427  	0.07194585 		-0.62733911
# AGL2      		  -1.23588026  		1.17458431 		-0.17879764  	0.22079153
# AGL2b       		  -0.91460577  		0.31622926 		-0.32777662  	0.17076256
# AGL3         		  -1.08805480 		-0.30362759 	-0.13971536  	0.93297937
# AGL4            	-1.07959183 			-1.21154002 		-0.06304717 	-0.32239766
# AGL5            	-0.91816684  		0.04730828  	0.53131992 		-0.11299118
# AGL6            	-0.47405196 		-0.07268660  	0.22374558  		0.06044278
# AGL7            	-0.79704351 			-0.53651785  	0.17727187  		0.25199746
# AGL8            	-1.13963854 			-0.74566175  	0.25011977 		-0.06144807
# Paracas1        -0.01192517 			-0.15893207  	0.14539943  		0.07954916
# Paracas2         0.28728540 		-0.12240146  	0.10699929  	0.27832201
# PiscoHolocene_1  1.41159299  0.36137901  	0.88714793  		0.13083500
# PiscoHolocene_2  1.48186019  0.22028032 -0.47937223 	-0.05447134
# PiscoHolocene_3  1.47944315  0.25797974 	-0.53055362 	-0.01236784
# PiscoHolocene_4  1.43066704  0.03825223 -0.58962868 	-0.05809486
# PiscoHolocene_5  1.54636134  0.15819634 	-0.20697306 	-0.37366338
# PiscoHolocene_6  1.73349258  0.69489705  0.07691685  	0.04934570
# PiscoHolocene_7  1.58068164  0.13646431 	-0.07420633  	0.31013997

#####################################################################
## Chao Dissimilarity Index (CDI)   ## All counts [species*samples]
######################################################################

# colors for plots 
jet.colors <- colorRampPalette(c("#0C0C38", "#0000FF", "#00CED1", "#CAFF70", "#FFD700", "#FFA500", "#FF0000", "#CD2626"))
clrs <- colorRampPalette(c("#440154FF","#3B528BFF","#21908CFF","#5DC863FF","greenyellow"),alpha=T)
# transparent colors
t_col <- function(color, percent = 50, name = NULL) {
	rgb.val <- col2rgb(color)
	t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
             max = 255, alpha = (100 - percent) * 255 / 100, names = name)
	invisible(t.col)
}   


chaod <- dis.chao(fossil, index='sorensen') 
mchaod <- as.matrix(chaod)

# plot distances
image(mchaod,xaxt="n",yaxt="n", col= clrs(50) )
axis(side=1, at= seq(0,1,length.out=length(colnames(mchaod))), labels=colnames(mchaod), cex.axis= 0.5,las=2) 
axis(side=2, at= seq(0,1,length.out=length(colnames(mchaod))), labels=colnames(mchaod), cex.axis= 0.5,las=2) 


# group by cluster
chclust(chaod, method = "coniss")-> conisclust
plot(conisclust, hang= -0.1, cex.axis=0.7,las=1, cex=0.7)
mtext("Depth (m)", side = 1, line = 3)

# Separate culled data in the TWO main groups
gr <- ifelse(cutree(conisclust,k=2)==1,"HOLOCENE", "MIOCENE")
g1 <- gr == "HOLOCENE"
g2 <- gr == "MIOCENE"

mchaod_g1 <- mchaod[g1,g1]
cdi_g1 <- mchaod_g1[upper.tri(mchaod_g1)]

mchaod_g2 <- mchaod[g2,g2]
cdi_g2 <- mchaod_g2[upper.tri(mchaod_g2)]

cdi_g1g2 <- as.numeric(mchaod[g1,g2])

cdi.g1 = cbind(cdi_g1,rep("G1",length(cdi_g1)))
cdi.g2 = cbind(cdi_g2,rep("G2",length(cdi_g2)))
cdi.1_2 = cbind(cdi_g1g2,rep("G1G2",length(cdi_g1g2)))
cdi = as.data.frame(rbind(cdi.g1,cdi.g2,cdi.1_2))

# Evaluate diferences between groups
cdi$chaod <- as.numeric(as.character(cdi[,1]))
anova_cdi = aov(cdi$chaod ~ cdi[,2])
summary(anova_cdi)

       ##    Df Sum Sq Mean Sq F value Pr(>F)    
##  cdi[, 2]      2  3.506  1.7530   108.1 <2e-16 ***
##  Residuals   150  2.432  0.0162                   
##  ---
##  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_cdi <- TukeyHSD(anova_cdi)
tukey_cdi

  ##    Tukey multiple comparisons of means
  ##      95% family-wise confidence level
  ##  Fit: aov(formula = cdi$chaod ~ cdi[, 2])
  ##  
  ##  $`cdi[, 2]`
  ##                 diff        lwr         upr     p adj
  ##  G1G2-G1  0.28037927  0.2200045  0.34075401 0.0000000    ##  
  ##  G2-G1   -0.04293161 -0.1139744  0.02811116 0.3278639  ##  
  ##  G2-G1G2 -0.32331088 -0.3836856 -0.26293614 0.0000000  ##  

#plot Distances within and between groups

plot(0, xlim=c(-0.1,1),ylim=c(0,6.5),t="n", xlab= "Chao dissimilarity index", ylab="frequency",las=1,cex.axis=0.5)
polygon(density(cdi_g1,from=-0.3,to=1.2), col=t_col("navy"),border=NA)
polygon(density(cdi_g1g2,from=-0.3,to=1.2), col=t_col("orange"),border=NA)
polygon(density(cdi_g2,from=-0.3,to=1.2), col=t_col("darkcyan"),border=NA)
legend("topright",c("Holocene","Miocene","between groups"),col=c(t_col("navy"),t_col("darkcyan"),t_col("orange")),bty="n", pch=15)


#####################################################################
## habit growth analyses  ## 
######################################################################

# Define USGS colors 
colors <- c("#FCD5B5", "gold") # Create a vector of colors for the Tortonian and Holocene

gr <- c(rep("Miocene",9),rep("Holocene",9)) ###Changes in proportions calculated using complete database 
#gr <- c(rep("Miocene",8),rep("Holocene",9)) #excluding samples with counts <100 grains/sample 


habit <- as.data.frame(read_excel("summary modern data.xlsx", sheet = 2))
habit.ptg <- as.data.frame(t(habit[-(1:16),2:19]))
colnames(habit.ptg)<- habit[-(1:16),1]
rownames(habit.ptg)<- habit[1,2:19]
# NA's = zero
habit.ptg[is.na(habit.ptg)] <- 0
chars <- sapply(habit.ptg, is.character)
habit.ptg[ , chars] <- as.data.frame(apply(habit.ptg[ , chars], 2, as.numeric)) #convert all character columns to numeric
data.habit<-cbind(habit.ptg,gr)


#######################################################
#Plot by taxonomic grps
##################################################################

# Set up plots
plot1 <- ggplot(data.habit, aes(x = gr, y = habit.ptg$Angiosperm, fill = gr)) +
  geom_boxplot(alpha = 0.5, size = 0.5, outlier.shape = NA) +  geom_point() + ylim(0, 85)+
    scale_fill_manual(values = colors)+
  geom_jitter(alpha = 0.5, size = 2, width = 0.2) 

plot2 <- ggplot(data.habit, aes(x = gr, y = habit.ptg$Gymnosperm, fill = gr)) +
  geom_boxplot(alpha = 0.5, size = 0.5, outlier.shape = NA) +
      geom_point() + ylim(0, 10)+
    scale_fill_manual(values = colors)+
  geom_jitter(alpha = 0.5, size = 2, width = 0.2, shape= 3) 

plot3 <- ggplot(data.habit, aes(x = gr, y = habit.ptg$FERNS, fill = gr)) +
  geom_boxplot(alpha = 0.5, size = 0.5, outlier.shape = NA) +
      geom_point() + ylim(0, 40)+
    scale_fill_manual(values = colors)+
  geom_jitter(alpha = 0.5, size = 2, width = 0.2, shape= 10) 


plot4 <- ggplot(data.habit, aes(x = gr, y = habit.ptg$unknwn.gp, fill = gr)) +
  geom_boxplot(alpha = 0.5, size = 0.5, outlier.shape = NA) +
    geom_point() + ylim(0, 40)+
    scale_fill_manual(values = colors)+
  geom_jitter(alpha = 0.5, size = 2, width = 0.2, shape= 20) 


# Add titles and remove axis labels
plot1 <- plot1 + ggtitle("Angiosperms") + ylab(NULL) + xlab(NULL)
plot2 <- plot2 + ggtitle("Gymnosperms") + ylab(NULL) + xlab(NULL)
plot3 <- plot3 + ggtitle("Pteridophytes") + ylab(NULL) + xlab(NULL)
plot4 <- plot4 + ggtitle("Unknown taxonomic group") + ylab(NULL) + xlab(NULL)

# Remove legend
plot1 <- plot1 + guides(fill = FALSE)
plot2 <- plot2 + guides(fill = FALSE)
plot3 <- plot3 + guides(fill = FALSE)
plot4 <- plot4 + guides(fill = FALSE)

# Combine plots, remove grid lines in strips and save
tx_gp_plot<-grid.arrange(plot1, plot2, plot3, plot4, top = "Comparison by Taxonomic Group", nrow=1) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank())

#ggsave("by taxon gps.pdf", plot=tx_gp_plot, dpi = 500)

#######################################################
##Plot by habit
##################################################################

# Set up plots
plot.h1 <- ggplot(data.habit, aes(x = gr, y = habit.ptg$Herbs.EXC, fill = gr)) +
  geom_boxplot(alpha = 0.5, size = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = colors)+
  geom_jitter(alpha = 0.5, size = 1.5, width = 0.1, col="black", shape=3) +


plot.h2 <- ggplot(data.habit, aes(x = gr, y = habit.ptg$Shrubs.EXC, fill = gr)) +
  geom_boxplot(alpha = 0.5, size = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = colors)+
      ylim(0, 10) +
  geom_jitter(alpha = 0.5, size = 1.5, width = 0.1, col="black", shape=3) +
  

plot.h3 <- ggplot(data.habit, aes(x = gr, y = habit.ptg$Trees.EXC, fill = gr)) +
  geom_boxplot(alpha = 0.5, size = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = colors)+
    ylim(0, 20) +
  geom_jitter(alpha = 0.5, size = 1.5, width = 0.1, col="black", shape=3) +


plot.h4 <- ggplot(data.habit, aes(x = gr, y = habit.ptg$Unknown.Habit, fill = gr)) +
  geom_boxplot(alpha = 0.5, size = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = colors)+
  geom_jitter(alpha = 0.5, size = 2, width = 0.2) +


plot.h5 <- ggplot(data.habit, aes(x = gr, y = habit.ptg$Herbs, fill = gr)) +
  geom_boxplot(alpha = 0.5, size = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = colors)+
  geom_jitter(alpha = 0.5, size = 2, width = 0.2) +


plot.h6 <- ggplot(data.habit, aes(x = gr, y = habit.ptg$Shrubs, fill = gr)) +
  geom_boxplot(alpha = 0.5, size = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = colors)+
  geom_jitter(alpha = 0.5, size = 2, width = 0.2) +

  
plot.h7 <- ggplot(data.habit, aes(x = gr, y = habit.ptg$Trees, fill = gr)) +
  geom_boxplot(alpha = 0.5, size = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = colors)+
  geom_jitter(alpha = 0.5, size = 2, width = 0.2) +



# Add titles and remove axis labels
plot.h1 <- plot.h1 + ggtitle("Herbs (exclusive)") + ylab(NULL) + xlab(NULL)
plot.h2 <- plot.h2 + ggtitle("Shrubs (exclusive)") + ylab(NULL) + xlab(NULL)
plot.h3 <- plot.h3 + ggtitle("Trees (exclusive)") + ylab(NULL) + xlab(NULL)
plot.h4 <- plot.h4 + ggtitle("Unknown habit") + ylab(NULL) + xlab(NULL)
plot.h5<- plot.h5 + ggtitle("Herbs") + ylab(NULL) + xlab(NULL)
plot.h6 <- plot.h6 + ggtitle("Shrubs") + ylab(NULL) + xlab(NULL)
plot.h7 <- plot.h7 + ggtitle("Trees") + ylab(NULL) + xlab(NULL)

# Remove legend
plot.h1 <- plot.h1 + guides(fill = FALSE)
plot.h2 <- plot.h2 + guides(fill = FALSE)
plot.h3 <- plot.h3 + guides(fill = FALSE)
plot.h4 <- plot.h4 + guides(fill = FALSE)
plot.h5 <- plot.h5 + guides(fill = FALSE)
plot.h6 <- plot.h6 + guides(fill = FALSE)
plot.h7 <- plot.h7 + guides(fill = FALSE)

# Combine plots, remove grid lines in strips and save
Habit_plot<-grid.arrange(plot.h1, plot.h2, plot.h3, plot.h4, plot.h5, plot.h6, plot.h7, plot.h4, 
top = "Comparison by Habit", heights = c(1, 1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank())

#ggsave("habit.pdf", plot=tx_gp_plot, dpi = 500)

######################################################################
####   Analyss using only counts of species with known affinity
######################################################################

conteo <- as.data.frame(read_excel("AGL Coex counts_2023 only pollen.xlsx" ,col_names = FALSE)) ## counts [species*samples] using  only species wiht known affinity
conteo1 <- as.data.frame(t(conteo[-(1:2),-1]))
colnames(conteo1)<- conteo[-(1:2),1]
rownames(conteo1)<- conteo[2,-1]
time <- (conteo[1,-1])

fossil <- as.data.frame(lapply(conteo1,function (x) as.numeric(as.character(x))))
colnames(fossil)<- conteo[-(1:2),1]
rownames(fossil)<- conteo[2,-1]
# NA's = zero
fossil[is.na(fossil)] <- 0

fossil.t<-t(fossil)##species in rows

Miocene.data<-fossil.t[, 1:9]
Holocene.data<-fossil.t[, 10:ncol(fossil.t)]
dim(Miocene.data) ##should be #sps=80 by 9 samples
dim(Holocene.data) #should be #sps=80 by 9 samples
samples<-c(1:9) #same for Holocene and Miocene, as similar number of samples


##Rarefaction Miocene vs Holocene 
########## code for range-through analysis
fill.occur=function(sp)
{
occur=which(sp>0)  ##array of row numbers where sp>0
fad=occur[1] #row of first number in array
numboccur=length(occur)# how many numbers in the array
lad=occur[numboccur]#row number of the last position in the array
alloccur=rep(0,length(sp))## produces a matrix of zeroes that is similar in size to original
i=1:length(sp)## all positions in the matrix are now labeled i
alloccur[i>=fad & i<=lad]=1## replaces all i in between fad and lad by 1
return(alloccur)#gives the matrix out
}


Mio.Rthrough<-t(apply(Miocene.data,1,fill.occur))## matrix after range through, samples in columns
Mio.total<-apply(Mio.Rthrough,2,sum)##diversity sum

Hol.Rthrough <-t(apply(Holocene.data,1,fill.occur))## matrix after range through, samples in columns
Hol.total<-apply(Hol.Rthrough,2,sum)##diversity sum

leg.txt<-c("Miocene", "Holocene")
plot(samples,  Mio.total, main="Range-through standing diversity", col="dark blue", pch=20, xlab="observed diverstity", ylab="Samples")
points(samples, Hol.total, col="brown", pch=20)
legend("topright", leg.txt, pch =c(20), col = c("dark blue", "brown"), cex=0.7)


## DROPING OUT SINGLETONES 

Mio.summa1<-apply(Miocene.data,1,sum)
Hol.summa1<-apply(Holocene.data,1,sum)

Mio.uniq<-Miocene.data[which(Mio.summa1>1),] #44 sp
Hol.uniq<-Holocene.data[which(Hol.summa1>1),] #24 sp

MioRT.uniq<-t(apply(Mio.uniq,1,fill.occur))## matrix after range through, samples in columns
Mio.unique<-apply(MioRT.uniq,2,sum)##diversity sum data without singles

HolRT.uniq <-t(apply(Hol.uniq,1,fill.occur))## matrix after range through, samples in columns
Holo.unique<-apply(HolRT.uniq,2,sum)##diversity sum data without singles

leg.txt2<-c("Miocene", "Holocene", "Miocene No singletones", "Holocene No singletones")
plot(samples,  Mio.total, main="Range-through diversity", col="dark blue", pch=20, xlab="observed diverstity", ylab="Samples", ylim=c(0, 45))
points(samples, Hol.total, col="brown", pch=20)
points(samples,  Mio.unique, col="dark blue", pch=18, )
points(samples, Holo.unique, col="brown", pch=18)
legend("topright", leg.txt2, pch =c(20, 20, 18, 18), col = c("dark blue", "brown", "dark blue", "brown"), cex=0.7)

### selecting most common elements per Flora (above 80 grains total)
Mio.comm<-Miocene.data[which(Mio.summa1>80),] #44 sp
Hol.comm<-Holocene.data[which(Hol.summa1>80),] #24 sp

###Within-Sample Diversity (Rarefaction) using  only species with known affinity

abundancedata <- fossil.t 
age.type <- c(1:ncol(fossil.t))
age.type <-(age.type>9)*1 ##### ### 1 is modern and 0 is fossil-AGL
age.type<-as.vector(age.type) ####  1 is modern and 0 is fossil-AGL

rar.50<-rarefy(t(abundancedata),50,se= TRUE)
rar.50test=rar.50[1,which(rar.50[2,]>0)]
Samp.50<-samples[which(rar.50[2,]>0)]
age.type50<-age.type[which(rar.50[2,]>0)]

rar.100<-rarefy(t(abundancedata),100,se= TRUE)
rar.100test=rar.100[1,which(rar.100[2,]>0)]
Samp.100<-samples[which(rar.100[2,]>0)]
age.type100<-age.type[which(rar.100[2,]>0)]

rar.150<-rarefy(t(abundancedata),150,se= TRUE)
rar.150test=rar.150[1,which(rar.150[2,]>0)]
Samp.150<-samples[which(rar.150[2,]>0)]
age.type150<-age.type[which(rar.150[2,]>0)]

plot(Samp.100, rar.100test, xlab="Miocene (sites 1 to 8) versus Holocene (sites 9 to 15)", ylab="Number of Species", pch=19, main="Rarefied richness", ylim=c(13,33))
points(Samp.50, rar.50test, col="blue", pch=21)
legend("topright", legend= c("Rarefaction at 100","Rarefaction at 50"), pch=c(19,21), col=c("black", "blue"),  bg="white")

par(mfrow=c(1,2))
boxplot(rar.100test ~ age.type100, col = "lightgray", 
main = "Richness estimation Rar.100",ylim=c(0,30),
xlab = "Miocene (0) versus Holocene (1)")

boxplot(rar.150test ~ age.type150, col = "lightgray", 
main = "Richness estimation Rar.150", ylim=c(0,30),
xlab = "Miocene (0) versus Holocene (1)")



######Tproof:

test.rar50<-t.test(rar.50test ~ age.type50)
test.rar50

#Welch Two Sample t-test
# data:  rar.50test by age.type50
# t = 1.8821, df = 13.692, p-value = 0.08126
# alternative hypothesis: true difference in means between group 0 and group 1 is not equal to 0
# 95 percent confidence interval:
#  -0.320652  4.836760
# sample estimates:
# mean in group 0 mean in group 1 
#       12.044857        9.786803 


test.rar100<-t.test(rar.100test~age.type100)
test.rar100

# Welch Two Sample t-test
# data:  rar.100test by age.type100
# t = 3.4371, df = 5.3825, p-value = 0.01645
# alternative hypothesis: true difference in means between group 0 and group 1 is not equal to 0
# 95 percent confidence interval:
#   1.632819 10.560989
# sample estimates:
# mean in group 0 mean in group 1 
#        16.26906        10.17216 


test.rar150<-t.test(rar.150test~age.type150)
test.rar150

# Welch Two Sample t-test

# data:  rar.150test by age.type150
# t = 3.2267, df = 4.1651, p-value = 0.03027
# alternative hypothesis: true difference in means between group 0 and group 1 is not equal to 0
# 95 percent confidence interval:
#  1.457667 17.620190
# sample estimates:
# mean in group 0 mean in group 1 
#       20.46798        10.92905 

# Extract the p-value
p_value <- test.rar150 $p.value

# Interpret the results
if (p_value < 0.05) {
  print("There are significant differences within the group.")
} else {
  print("There are no significant differences within the group.")
}


##   Simpson index (samples above 80 grains)

species_data <- fossil[rowSums(fossil) > 80, ]
age.type80<-age.type[rowSums(fossil) > 80]
simpson_index <- diversity(species_data, index = "simpson")

# Print the Simpson index
print(simpson_index)

Simpson.test<-t.test(simpson_index ~ age.type80)
Simpson.test

#	Welch Two Sample t-test
# data:  simpson_index by age.type80
# t = 1.3358, df = 9.5646, p-value = 0.2125
# alternative hypothesis: true difference in means between group 0 and group 1 is not equal to 0
# 95 percent confidence interval:
# -0.02421204  0.09559675
# sample estimates:
# mean in group 0 mean in group 1 
 #      0.8172509       0.7815586 

# Extract the p-value
p_value <- Simpson.test $p.value

# Interpret the results
if (p_value < 0.05) {
  print("There are significant differences within the group.")
} else {
  print("There are no significant differences within the group.")
}


###DCA analyses using all dataset

decorana(t(abundancedata))->deco

plot(deco)
plot(deco, display="sites")####samples
plot(deco, display="species")####species
plot(deco, display="both")####to see both.. samples and species

# Call: decorana(veg = t(abundancedata)) 
# Detrended correspondence analysis with 26 segments.
# Rescaling of axes with 4 iterations.
# Total inertia (scaled Chi-square): 2.1016 

 #                       					DCA1   	DCA2   	DCA3    	DCA4
# Eigenvalues         			 0.7140 	0.2009 	0.1600 	0.09635
# Additive Eigenvalues 	0.7140 	0.1967 	0.1588 	0.08134
# Decorana values      		0.7264 	0.1937 	0.1282 	0.04761
# Axis lengths         			3.5170 	1.4286 	1.6216 	1.06109

scores(deco)->values
DCA1<-values[,1] ## DCA1 0.714
DCA2<-values[,2] ## DCA2 0.20
all.samples<-c(1:18)
## DCA versus samples
plot(DCA1[1:9], all.samples[1:9], main="DCA Miocene versus Holocene", col="dark blue", pch=20, xlim=c(-2,2), ylim=c(0,18), xlab="DCA  axis 1", ylab="Samples") ## Miocene samples
points(DCA1[10:18], all.samples[10:18], main="DCA Miocene versus Holocene", col="brown", pch=18) ## Holocene samples
legend("topleft", legend= c("Miocene","Holocene"), pch=c(20,18), col=c("dark blue", "brown"),  bg="white")

##selecting species driginc DCA ordination
max(deco$cproj[,1]) ##max socre in DCA1
min(deco$cproj[,1]) ##min socre in DCA1
highest.scores<-deco$cproj>=2
table(highest.scores[,1])
hs.species <- rownames(deco$cproj[deco$rproj[, 1] >= 2, ]) ##high scored sp
hs.scores <- deco$cproj[hs.species, ]

ls.species <- rownames(deco$cproj[deco$rproj[, 1] <= 0.5, ]) ##high scored sp
ls.scores <- deco$cproj[ls.species, ]

######