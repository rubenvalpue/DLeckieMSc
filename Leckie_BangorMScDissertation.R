#13/05/2021
#Duncan Leckie
#M.S.c Thesis Analysis
#How does CWD Characteristics affect fungal Diversity?
#What moderators have a significant affect on fungal Diversity?
#Following Steps from http://environmentalcomputing.net/meta-analysis-1/
#Install the metafor package 
#install.packages("metafor")
#Install the dplyr package 
#install.packages("dplyr")
#Install plotrix package
#install.packages("plotrix")
#insatll read excel package
#install.packages("readxl")
##############################################################################
#Start of Script when opening again

#clean memory
rm(list=ls())


#load Packages 
library(metafor)
library(dplyr)
library(plotrix)
library(readxl)
library(ape)

# load the MA1 data
library(readxl)
dat3 <- MA1_CWD_Volume <- read_excel(paste0(getwd(),"/Data/MA1_CWD_Volume.xlsx"), 
                             col_types = c("numeric", "text", "numeric", 
                                           "text", "text", "text", "numeric", 
                                           "text", "text", "numeric", "numeric", 
                                           "numeric", "numeric", "numeric", 
                                           "text", "text", "text", "text"))
View(dat3)

# ANALYSIS
# VISUALISATION OF EFFECT SIZES PER STUDY
# Visualize data$yi effect sizes of the different studies
# following Q7 and fig.5d of https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-017-0357-7#Fig5
# Visualisation forest plot
# "Forest plot visualizes point estimates (effect sizes) and their 95% confidence intervals (CI) 
# (based on sampling error variance) by using the forest function"

# Create a forest plot
forest(dat3$Zr, dat3$SEV, main="Standarised Effect Size MA1", slab=paste(dat3$In_Ref))

# Create a basic robust model

# Adding weights to account for publication bias
# vi = variance of the effect sizes per study
# We make weights, which is 1/vi and stick that into the argument, weights
dat3[, "wi"] <- 1/dat3$SEV
View(dat3)

# This model uses the weights (make it robust) and a multilevel model (paper bias)
# Model without moderators
robustml_m <- rma.mv(yi = Zr, V = SEV, W = wi, random = list(~1 | Paper, ~1 | ID), method = "REML", data = dat3)
sum_basicmodel <- summary(robustml_m) 
sum_basicmodel
# AIC = 31.3216, Q=87.7453, p-val < .0001
# estimate      se    zval    pval   ci.lb   ci.ub 
#  0.2743    0.1101  2.4904  0.0128  0.0584  0.4901   
# p=value = 0.0128 < 0.05 so there is a significant effect of CWD Volume on fungal diversity.

# HETEROGENEITY BASIC MODEL

# Check heterogeneity of the basic model
# Calculate the I^2 = (sigma2.1 + sigma2.2)/(sigma2.1 + sigma2.2 + mean_v)
# estimate v as the mean of the vi = average mean variance
# Save sigma2.1 and 2.2 for heterogeneity calculation

sigmasBasic <- sum_basicmodel$sigma2
sigmaBasic2.1 <- sigmasBasic[1]
sigmaBasic2.2 <- sigmasBasic[2]
mean_SEV <- mean(dat3$SEV)
I2Basic <- (sigmaBasic2.1 + sigmaBasic2.2)/(sigmaBasic2.1 + sigmaBasic2.2 + mean_SEV)
I2Basic
#I^2 is 0.5235751 which is 52.36% 
# This means that 52.36% of the variation in effect sizes is due to the between-study variance
# This is Medium Heterogeneity 

# Moderators 
#A meta regression on selected moderators will be run the result will explain the amount of heterogeneity accounted for by the moderators. a p-value will be given for each moderator indicating significance level
#Add the moderators to the Model 

#Moderator Biome
robustml_Biome <- rma.mv(yi = Zr, V = SEV, W = wi, random = list(~1 | Paper, ~1 | ID), mod = ~Biome, method = "REML", data = dat3)
sum_Biome <- summary(robustml_Biome)
sum_Biome
#AIC 32.2184
#Test of Moderators (Coefficient 2)
#QM(df = 1) = 0.1474, p-val = 0.7011  
#p-value of test of moderators is not significant, > 0.05, so Biome does not significantly account for any variation in the data


# HETEROGENEITY Biome
# Calculate the I^2 = (sigma2.1 + sigma2.2)/(sigma2.1 + sigma2.2 + mean_v)
# estimate v as the mean of the vi = average mean variance
# Save sigma2.1 and 2.2 for heterogeneity calculation
sigmasBiome <- sum_Biome$sigma2
sigmaBiome2.1 <- sigmasBiome[1]
sigmaBiome2.2 <- sigmasBiome[2]
mean_SEV <- mean(dat3$SEV)
I2Biome <- (sigmaBiome2.1 + sigmaBiome2.2)/(sigmaBiome2.1 + sigmaBiome2.2 + mean_SEV)
I2Biome
# I^2 is 0.536737 = 53.67% this is Medium Heterogeneity 

# VISUALISATION Biome
# Visualize biome using boxplot (PlotCI)
# Save the estimates and CI for biome
xbiome <- 1:2
biome_estimate <- sum_Biome$beta
biome_ci.lb <- sum_Biome$ci.lb
biome_ci.ub <- sum_Biome$ci.ub
biome_Names <- c("Boreal","Temperate")

require(plotrix)
plotCI(xbiome,biome_estimate, ui=biome_ci.ub, li=biome_ci.lb,
       ylab="Estimate with confidence intervals", xlab="Biome", 
       main="Estimate with confidence intervals for biome", xaxt = 'n', cex.axis=1.3,cex.lab=1.3,cex.main=1.5)
axis(1, at=1:2,labels=biome_Names, cex.lab=1.2)
# Add a horizontal line at y=0 in color blue
abline(h=0, col='blue')


# Moderator For_Type
robustml_MA1_For_Type <- rma.mv(yi = Zr, V = SEV, W = wi, random = list(~1 | Paper, ~1 | ID), mod = ~For_Type , method = "REML", data = dat3)
sum_For_Type  <- summary(robustml_MA1_For_Type )
sum_For_Type
#AIC = 36.6193
#Test of Moderators (coefficient 2):
#QM(df = 2) =0.0823, p-val = 0.9597
#p-value of test of moderators is not significant, > 0.05, so For_type does not significantly account for variation in the data

#                  estimate      se     zval    pval    ci.lb   ci.ub 
#intrcpt              0.3150  0.2016   1.5627  0.1181  -0.0801  0.7101 
#For_TypeDeciduous   -0.0485  0.2656  -0.1828  0.8550  -0.5691  0.4720 
#For_TypeMixed       -0.1056  0.3857  -0.2737  0.7843  -0.8615  0.6504 

# HETEROGENEITY Type
# Calculate the I^2 = (sigma2.1 + sigma2.2)/(sigma2.1 + sigma2.2 + mean_v)
# estimate v as the mean of the vi = average mean variance
# Save sigma2.1 and 2.2 for heterogeneity calculation
sigmasFor_Type <- sum_For_Type$sigma2
sigmaFor_Type2.1 <- sigmasFor_Type[1]
sigmaFor_Type2.2 <- sigmasFor_Type[2]
mean_SEV <- mean(dat3$SEV)
I2For_Type <- (sigmaFor_Type2.1 + sigmaFor_Type2.2)/(sigmaFor_Type2.1 + sigmaFor_Type2.2 + mean_SEV)
I2For_Type
#I^2 is 0.5998545 = 59.99% this is Medium Heterogeneity


# VISUALISATION For_Type
# Visualize biome using boxplot (PlotCI)
# Save the estimates and CI for Type
xFor_Type <- 1:3
For_Type_estimate <- sum_For_Type$beta
For_Type_ci.lb <- sum_For_Type$ci.lb
For_Type_ci.ub <- sum_For_Type$ci.ub
For_Type_Names <- c("Coniferious","Deciduous","Mixed")

require(plotrix)
plotCI(xFor_Type,For_Type_estimate, ui=For_Type_ci.ub, li=For_Type_ci.lb,
       ylab="Estimate with confidence intervals", xlab="For_Type", 
       main="Estimate with confidence intervals for Forest Type", xaxt = 'n', cex.axis=1.3,cex.lab=1.3,cex.main=1.5)
axis(1, at=1:3,labels=For_Type_Names, cex.lab=1.2)
# Add a horizontal line at y=0 in color blue
abline(h=0, col='blue')

# Moderator CWD_spp
robustml_MA1_CWD_spp <- rma.mv(yi = Zr, V = SEV, W = wi, random = list(~1 | Paper, ~1 | ID), mod = ~CWD_spp , method = "REML", data = dat3)
sum_CWD_spp  <- summary(robustml_MA1_CWD_spp )
sum_CWD_spp

# HETEROGENEITY CWD_spp
# Calculate the I^2 = (sigma2.1 + sigma2.2)/(sigma2.1 + sigma2.2 + mean_v)
# estimate v as the mean of the vi = average mean variance
# Save sigma2.1 and 2.2 for heterogeneity calculation
sigmasCWD_spp <- sum_CWD_spp$sigma2
sigmaCWD_spp2.1 <- sigmasCWD_spp[1]
sigmaCWD_spp2.2 <- sigmasCWD_spp[2]
mean_SEV <- mean(dat3$SEV)
I2CWD_spp <- (sigmaCWD_spp2.1 + sigmaCWD_spp2.2)/(sigmaCWD_spp2.1 + sigmaCWD_spp2.2 + mean_SEV)
I2CWD_spp

# VISUALISATION CWD_spp
# Visualize biome using boxplot (PlotCI)
# Save the estimates and CI for Type
xCWD_spp <- 1:5
CWD_spp_estimate <- sum_CWD_spp$beta
CWD_spp_ci.lb <- sum_CWD_spp$ci.lb
CWD_spp_ci.ub <- sum_CWD_spp$ci.ub
CWD_spp_Names <- c("Beech","Mix","Oak","Pine","Spruce")

require(plotrix)
plotCI(xCWD_spp,CWD_spp_estimate, ui=CWD_spp_ci.ub, li=CWD_spp_ci.lb,
       ylab="Estimate with confidence intervals", xlab="CWD_spp", 
       main="Estimate with confidence intervals for CWD_spp", xaxt = 'n', cex.axis=1.3,cex.lab=1.3,cex.main=1.5)
axis(1, at=1:5,labels=CWD_spp_Names, cex.lab=1.2)
# Add a horizontal line at y=0 in color blue
abline(h=0, col='blue')


##################################################################################################
# Moderator FO Type
robustml_MA1_Type <- rma.mv(yi = Zr, V = SEV, W = wi, random = list(~1 | Paper, ~1 | ID), mod = ~Type , method = "REML", data = dat3)
sum_Type  <- summary(robustml_MA1_Type )
sum_Type
#AIC = 39.0493
#Test of Moderators (coefficient 2):
#QM(df = 3) = 5.2958, p-val = 0.1514
#p-value of test of moderators is not significant, > 0.05, so Fo Type does not significantly account for variation in the data

# HETEROGENEITY FO Type
# Calculate the I^2 = (sigma2.1 + sigma2.2)/(sigma2.1 + sigma2.2 + mean_v)
# estimate v as the mean of the vi = average mean variance
# Save sigma2.1 and 2.2 for heterogeneity calculation
sigmasType <- sum_Type$sigma2
sigmaType2.1 <- sigmasType[1]
sigmaType2.2 <- sigmasType[2]
mean_SEV <- mean(dat3$SEV)
I2Type <- (sigmaType2.1 + sigmaType2.2)/(sigmaType2.1 + sigmaType2.2 + mean_SEV)
I2Type
#I^2 is 0.5710074= 57.10% this is Medium Heterogeneity


# VISUALISATION Type
# Visualize biome using boxplot (PlotCI)
# Save the estimates and CI for Type
xType <- 1:4
Type_estimate <- sum_Type$beta
Type_ci.lb <- sum_Type$ci.lb
Type_ci.ub <- sum_Type$ci.ub
Type_Names <- c("Boreal_Coniferious","Boreal_Deciduous","Temperate_Deciduous","Boreal_Mixed")

require(plotrix)
plotCI(xType,Type_estimate, ui=Type_ci.ub, li=Type_ci.lb,
       ylab="Estimate with confidence intervals", xlab="Type", 
       main="Estimate with confidence intervals for Forest Type", xaxt = 'n', cex.axis=1.3,cex.lab=1.3,cex.main=1.5)
axis(1, at=1:4,labels=Type_Names, cex.lab=1.2)
# Add a horizontal line at y=0 in color blue
abline(h=0, col='blue')


##########################################################
# PUBLICATION BIAS MA1

# Check for publication bias using the funnel()
# Use to detect a funnel asymmetry (a sign of publication bias),
funnel(robustml_m, main="Funnel plot - Publication bias MA1")

##############################################################
#Time Lag Bias MA1
# Check for time-lag bias
# following http://www.metafor-project.org/doku.php/plots:cumulative_forest_plot

# Time-lag bias test
### Fit random-effects models for time-lag bias purposes only ()
res <- rma(Zr, SEV, data=dat3, slab=paste(Paper, Year, sep=", "))
### Cumulative meta-analysis (in the order of publication year)
tmp <- cumul(res, order=dat3$Year)
### cumulative forest plot
forest(tmp, xlim=c(-4,2), at=log(c(0.25, 0.5, 1, 2, 2.5, 3)),
       atransf=exp, digits=c(2,3), cex=0.75, header="Author(s) and Year MA1")


###########################################################################################
###########################################################################################
###########################################################################################

#Load MA2 Data 
library(readxl)
dat4 <- MA2_CWD_Type <- read_excel(paste0(getwd(),"/Data/MA2_CWD_Type.xlsx"), 
                           col_types = c("numeric", "text", "numeric", 
                                         "text", "text", "text", "text", "text", 
                                         "text", "numeric", "numeric", "text", 
                                         "numeric", "numeric", "numeric", 
                                         "text", "text", "text", "text", "numeric"))
View(dat4)

############################################################################################
# ANALYSIS
# VISUALISATION OF EFFECT SIZES PER STUDY
# Visualize data$yi effect sizes of the different studies
# following Q7 and fig.5d of https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-017-0357-7#Fig5
# Visualisation forest plot
# "Forest plot visualizes point estimates (effect sizes) and their 95% confidence intervals (CI) 
# (based on sampling error variance) by using the forest function"

# Create a forest plot
forest(dat4$Zr, dat4$SEV, main="Standarised Effect Size MA2", slab=paste(dat4$In_Ref))

# Create a basic robust model

# Adding weights to account for publication bias
# vi = variance of the effect sizes per study
# We make weights, which is 1/vi and stick that into the argument, weights
dat4[, "wi"] <- 1/dat4$SEV
View(dat4)

# This model uses the weights (make it robust) and a multilevel model (paper bias)
# Model without moderators
robustml_MA2 <- rma.mv(yi = Zr, V = SEV, W = wi, random = list(~1 | Paper, ~1 | ID), method = "REML", data = dat4)
sum_basicmodelma2 <- summary(robustml_MA2)
sum_basicmodelma2
#AIC = 46.2692  Q(df=12) = 52.7649, p-val < .0001
#pvalue= 0.0527  >0.05 so there is no significant effect of CWD Volume and fungal diversity. this is due to the smaller sample size. 
#estimate = 0.2797         se = 0.1443   zval = 1.9377   pval = 0.0527    ci.lb = -0.0032   ci.ub = 0.5626 

# HETEROGENEITY BASIC MODEL

# Check heterogeneity of the basic model
# Calculate the I^2 = (sigma2.1 + sigma2.2)/(sigma2.1 + sigma2.2 + mean_v)
# estimate v as the mean of the vi = average mean variance
# Save sigma2.1 and 2.2 for heterogeneity calculation

sigmasBasicma2 <- sum_basicmodelma2$sigma2
sigmaBasicma22.1 <- sigmasBasicma2[1]
sigmaBasicma22.2 <- sigmasBasicma2[2]
mean_SEVma2 <- mean(dat4$SEV)
I2Basicma2 <- (sigmaBasicma22.1 + sigmaBasicma22.2)/(sigmaBasicma22.1 + sigmaBasicma22.2 + mean_SEVma2)
I2Basicma2
#I^2 = 0.2219574 = 22.2% of the variation in effect sizes is due to the between-study variance
# This is Low Heterogeneity

###########################################################

# Moderators 
#A meta regression on selected moderators will be run the result will explain the amount of heterogeneity accounted for by the moderators. a p-value will be given for each moderator indicating significance level
#Add the moderators to the Model 


# Moderator CWD_Type
robustml_MA2_CWD_Type <- rma.mv(yi = Zr, V = SEV, W = wi, random = list(~1 | Paper, ~1 | ID), mod = ~CWD_Type , method = "REML", data = dat4)
sum_CWD_Type  <- summary(robustml_MA2_CWD_Type)
sum_CWD_Type 
#AIC = 29.9257
#Test of Moderators (coefficient 2:3):
#QM(df = 2) =  5.3256, p-val = 0.0698
#p-value of test of moderators is not significant, > 0.05, so CWD_type does not significantly account for variation in the data

#        estimate      se     zval    pval    ci.lb   ci.ub 
#intrcpt          0.4068  0.2862   1.4214  0.1552  -0.1541  0.9678    
#CWD_TypeSnag    -0.1331  0.3604  -0.3693  0.7119  -0.8394  0.5732    
#CWD_TypeStump   -0.6544  0.3574  -1.8311  0.0671  -1.3548  0.0461  . 


# HETEROGENEITY CWD_Type
# Calculate the I^2 = (sigma2.1 + sigma2.2)/(sigma2.1 + sigma2.2 + mean_v)
# estimate v as the mean of the vi = average mean variance
# Save sigma2.1 and 2.2 for heterogeneity calculation
sigmasCWD_Type <- sum_CWD_Type$sigma2
sigmaCWD_Type2.1 <- sigmasCWD_Type[1]
sigmaCWD_Type2.2 <- sigmasCWD_Type[2]
mean_SEVma2 <- mean(dat4$SEV)
I2CWD_Type <- (sigmaCWD_Type2.1 + sigmaCWD_Type2.2)/(sigmaCWD_Type2.1 + sigmaCWD_Type2.2 + mean_SEVma2)
I2CWD_Type
#I^2 is 0.4599882 = 46% this is Medium Heterogeneity

# VISUALISATION CWD_Type
# Visualize CWD_Type using boxplot (PlotCI)
# Save the estimates and CI for CWD_Type
xCWD_Type <- 1:3
CWD_Type_estimate <- sum_CWD_Type$beta
CWD_Type_ci.lb <- sum_CWD_Type$ci.lb
CWD_Type_ci.ub <- sum_CWD_Type$ci.ub
CWD_Type_Names <- c("Fallen","Snag","Stump")

require(plotrix)
plotCI(xCWD_Type, CWD_Type_estimate, ui=CWD_Type_ci.ub, li=CWD_Type_ci.lb,
       ylab="Estimate with confidence intervals", xlab="CWD Type", 
       main="Estimate with confidence intervals for CWD Type", xaxt = 'n', cex.axis=1.3,cex.lab=1.3,cex.main=1.5)
axis(1, at=1:3,labels=CWD_Type_Names, cex.lab=1.2)
# Add a horizontal line at y=0 in color blue
abline(h=0, col='blue')

##################################################################################################
# PUBLICATION BIAS

# Check for publication bias using the funnel()
# Use to detect a funnel asymmetry (a sign of publication bias),
funnel(robustml_MA2, main="Funnel plot - Publication bias MA2")

##############################################################################################

# Time-lag bias test
### Fit random-effects models for time-lag bias purposes only ()
res2 <- rma(Zr, SEV, data=dat4, slab=paste(Paper, Year, sep=", "))
### Cumulative meta-analysis (in the order of publication year)
tmp2 <- cumul(res2, order=dat4$Year)
### cumulative forest plot
forest(tmp2, xlim=c(-4,2), at=log(c(0.25, 0.5, 1, 2, 2.5, 3 )),
       atransf=exp, digits=c(2,3), cex=0.75, header="Author(s) and Year MA2")

