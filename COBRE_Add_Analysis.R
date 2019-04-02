#Required libraries
library(R.matlab)
library(abind)
library(psych)
library(MASS)
library(ggplot2)
library(reshape2)
library(igraph)
library(brainGraph)

#Setting working directory (optional) and bringing in functions
source("COBRE_Functions.R")

#######################################
######### Full COBRE Dataset ##########
#######################################

#Bringing in COBRE parameters dataset
#NOTE: You may have to change the working directory and file name for your 
#      file of demographic variables from the COBRE dataset
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
#NOTE: You may have to change the working directory for the location of your .mat file
#      adjacency matrices. This is assuming pre-processing occured in the CONN toolbox of MatLab.
full_files<-list.files(pattern="*.mat")
full_mydata<-lapply(full_files,function(X) readMat(X)$Z)
full_corr<-array(dim=c(length(full_files),2,132,132))

##Adding the data to the empty arrays and adding empty degree cells
for (i in 1:nrow(full_corr)) {
  p<-nrow(full_corr[1,1,,])
  full_corr[i,1,,]<-full_mydata[[i]][,-(133:167)]
  full_corr[i,2,,]<-matrix(0,p,p)
}

full_corr_zero<-Matlab_fMRI(full_corr)

#Calculating global efficiency and rich club coefficient for each adjacency matrix,
#then adding to the cobre_full dataset
cobre_full$Global_Eff<-NA
cobre_full$Rich_Club<-NA

for (i in 1:nrow(full_corr_zero)){
  igraph_adjmat<-graph_from_adjacency_matrix(full_corr_zero[i,1,,],mode="undirected",weighted=TRUE,diag=FALSE)
  cobre_full$Global_Eff[i]<-efficiency(igraph_adjmat,type="global",weights=NULL)
  cobre_full$Rich_Club[i]<-rich_club_coeff(igraph_adjmat,k=1,weighted=FALSE)$`phi`
}

#Logistic regression - global efficiency
globeff_mod<-glm(Group~Age_Years+Gender+as.factor(Handedness)+Global_Eff,data=cobre_full,family="binomial")
summary(globeff_mod)

#Logistic regression - rich club coefficient
richclub_mod<-glm(Group~Age_Years+Gender+as.factor(Handedness)+Rich_Club,data=cobre_full,family="binomial")
summary(richclub_mod)

##########################################
######### Paranoid Only Dataset ##########
##########################################

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

par_corr_zero<-Matlab_fMRI(corr_par)

#Calculating global efficiency and rich club coefficient for each adjacency matrix,
#then adding to the cobre_full dataset
cobre_par$Global_Eff<-NA
cobre_par$Rich_Club<-NA

for (i in 1:nrow(par_corr_zero)){
  igraph_adjmat2<-graph_from_adjacency_matrix(par_corr_zero[i,1,,],mode="undirected",weighted=TRUE,diag=FALSE)
  cobre_par$Global_Eff[i]<-efficiency(igraph_adjmat2,type="global",weights=NULL)
  cobre_par$Rich_Club[i]<-rich_club_coeff(igraph_adjmat2,k=1,weighted=FALSE)$`phi`
}

#Logistic regression - global efficiency - can't have handedness in model because nearly everyone is left-handed
globeff_mod2<-glm(Group~Age_Years+Gender+Global_Eff,data=cobre_par,family="binomial")
summary(globeff_mod2)

#Logistic regression - rich club coefficient - can't have handedness in model because nearly everyone is left-handed
richclub_mod2<-glm(Group~Age_Years+Gender+Rich_Club,data=cobre_par,family="binomial")
summary(richclub_mod2)
