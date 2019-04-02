#Required libraries
library(R.matlab)
library(abind)
library(psych)
library(MASS)
library(ggplot2)
library(reshape2)

#Setting working directory (optional) and bringing in functions
source("COBRE_Functions.R")

################################################################
######### COBRE Dataset: Paranoid Schizophrenia only ##########
################################################################

#Bringing in COBRE parameters dataset
#NOTE: You may have to change the working directory and file name for your 
#      file of demographic variables from the COBRE dataset
cobre<-read.csv('COBRE_ParanoidOnly.csv')
head(cobre)

#Dropping the file name, Dx_Code, and Dx_Detail variables as they only pertain to one group
cobre_par<-cobre[c(-1,-7,-8)]
head(cobre_par)

#Creating dichotomous / numeric coding for categorical variables
cobre_par$Group<-ifelse(cobre_par$Group == "Control", 0, 
                       ifelse(cobre_par$Group == "Patient", 1, NA))
cobre_par$Gender<-ifelse(cobre_par$Gender == "Male", 0, 
                        ifelse(cobre_par$Gender == "Female", 1, NA))

#Bringing in the correlation matrix files (.mat file format)
#NOTE: You may have to change the working directory for the location of your .mat file
#      adjacency matrices. This is assuming pre-processing occured in the CONN toolbox of MatLab.
files_par<-list.files(pattern="*.mat")
mydata_par<-lapply(files_par,function(X) readMat(X)$Z)
corr_par<-array(dim=c(length(files_par),2,132,132))

##Adding the data to the empty arrays and adding empty degree cells
for (i in 1:nrow(corr_par)) {
  p<-nrow(corr_par[1,1,,])
  corr_par[i,1,,]<-mydata_par[[i]][,-(133:167)]
  corr_par[i,2,,]<-matrix(0,p,p)
}

#######################################
## Zeroing out negative correlations ##
#######################################

##Applying functions to correlation array to get RPD matrix
par_corr_zero<-Matlab_fMRI(corr_par)
par_corr_zero<-degree(par_corr_zero)
par_eff_res_zero<-resist(par_corr_zero)
par_RPD_matrix_zero<-distance(par_eff_res_zero)

bounds_par_zero<-up_low(par_RPD_matrix_zero)
U_pzero<-max(bounds_par_zero)*100
L_pzero<-min(bounds_par_zero[bounds_par_zero>0])*0.1

#####  Semiparametric score test  #####
score(cobre_par,par_RPD_matrix_zero,L_pzero,U_pzero,5000) #p=0.4909631

#####  Nonparametric score test  #####
score_null(cobre_par,par_RPD_matrix_zero,L_pzero,U_pzero,5000) #p=0.1496169

########################################
## Keeping negative correlations as is##
########################################

##Applying functions to correlation array to get RPD matrix
par_corr_neg<-Matlab_fMRI_neg(corr_par)
par_corr_neg<-degree(par_corr_neg)
par_eff_res_neg<-resist(par_corr_neg)
par_RPD_matrix_neg<-distance(par_eff_res_neg)

bounds_pneg<-up_low(par_RPD_matrix_neg)
U_pneg<-max(bounds_pneg)*100
L_pneg<-min(bounds_pneg[bounds_pneg>0])*0.1

#####  Semiparametric score test  #####
score(cobre_par,par_RPD_matrix_neg,L_pneg,U_pneg,5000) #p=0.5855368

#####  Nonparametric score test  #####
score_null(cobre_par,par_RPD_matrix_neg,L_pneg,U_pneg,5000) #p=0.4601876
