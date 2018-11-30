#Required libraries
library(R.matlab)
library(abind)
library(psych)
library(MASS)
library(ggplot2)
library(reshape2)

#Set working directory to bring in the COBRE functions
setwd("c:/Users/jenseale/Dropbox/MS_Thesis_Work/R_Code")
source("COBRE_Functions.R")

#######################################
######### Full COBRE Dataset ##########
#######################################

#Bringing in COBRE parameters dataset
setwd('C:/Users/jenseale/Dropbox/MS_Thesis_Work/')
cobre<-read.csv('COBRE_Phenotypic_Parameters.csv')
head(cobre)

#Dropping the Dx_Code and Dx_Detail variables as they only pertain to one group
cobre_full<-cobre[c(-6,-7)]
head(cobre_full)

#Dropping the following patients
##40045: unable to preprocess rs-fMRI dataset fully
##40070: disenrolled
##40083: disenrolled
cobre_full<-cobre_full[!(cobre_full$ID=="40045" | cobre_full$ID=="40070" | cobre_full$ID=="40083"),]

#Creating dichotomous / numeric coding for categorical variables
cobre_full$Group<-ifelse(cobre_full$Group == "Control", 0, 
                       ifelse(cobre_full$Group == "Patient", 1, NA))
cobre_full$Gender<-ifelse(cobre_full$Gender == "Male", 0, 
                        ifelse(cobre_full$Gender == "Female", 1, NA))

cobre_full$Age_Years<-as.numeric(levels(cobre_full$Age_Years))[cobre_full$Age_Years]


#Bringing in the correlation matrix files (.mat file format)
setwd("C:/Users/jenseale/Dropbox/MS_Thesis_Work/For_import_R/")
full_files<-list.files(pattern="*.mat")
full_mydata<-lapply(full_files,function(X) readMat(X)$Z)
full_corr<-array(dim=c(length(full_files),2,132,132))

##Adding the data to the empty arrays and adding empty degree cells
for (i in 1:nrow(full_corr)) {
  p<-nrow(full_corr[1,1,,])
  full_corr[i,1,,]<-full_mydata[[i]][,-(133:167)]
  full_corr[i,2,,]<-matrix(0,p,p)
}

#######################################
## Zeroing out negative correlations ##
#######################################

##Applying functions to correlation array to get RPD matrix
full_corr_zero<-Matlab_fMRI(full_corr)
full_corr_zero<-degree(full_corr_zero)
full_eff_res_zero<-resist(full_corr_zero)
full_RPD_matrix_zero<-distance(full_eff_res_zero)

bounds_zero<-up_low(full_RPD_matrix_zero)
U_zero<-max(bounds_zero)*100
L_zero<-min(bounds_zero[bounds_zero>0])*0.1

#####  Semiparametric score test  #####
score(cobre_full,full_RPD_matrix_zero,L_zero,U_zero,5000) #p=0.5344345

#####  Nonparametric score test  #####
score_null(cobre_full,full_RPD_matrix_zero,L_zero,U_zero,5000) #p=0.0182281

########################################
## Keeping negative correlations as is##
########################################

##Applying functions to correlation array to get RPD matrix
full_corr_neg<-Matlab_fMRI_neg(full_corr)
full_corr_neg<-degree(full_corr_neg)
full_eff_res_neg<-resist(full_corr_neg)
full_RPD_matrix_neg<-distance(full_eff_res_neg)

bounds_neg<-up_low(full_RPD_matrix_neg)
U_neg<-max(bounds_neg)*100
L_neg<-min(bounds_neg[bounds_neg>0])*0.1

#####  Semiparametric score test  #####
score(cobre_full,full_RPD_matrix_neg,L_neg,U_neg,5000) #p=0.6146681

#####  Nonparametric score test  #####
score_null(cobre_full,full_RPD_matrix_neg,L_neg,U_neg,5000) #p=0.4833895
