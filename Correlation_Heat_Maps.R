#Required libraries
library(R.matlab)
library(abind)
library(psych)
library(MASS)
library(ggplot2)
library(reshape2)
library(gridExtra)

#Setting working directory (optional) and bringing in functions
source("COBRE_Functions.R")

########################################
######### Full COBRE Dataset ##########
########################################

#Bringing in the correlation matrix files for controls (.mat file format)
#NOTE: You may have to change the working directory for the location of your .mat file
#      adjacency matrices. This is assuming pre-processing occured in the CONN toolbox of MatLab.
control_files<-list.files(pattern="*.mat")
control_data<-lapply(control_files,function(X) readMat(X)$Z)
control_mat<-array(dim=c(length(control_files),1,132,132))

#Bringing in the correlation matrix files for patients (.mat file format)
#NOTE: You may have to change the working directory for the location of your .mat file
#      adjacency matrices. This is assuming pre-processing occured in the CONN toolbox of MatLab.
case_files<-list.files(pattern="*.mat")
case_data<-lapply(case_files,function(X) readMat(X)$Z)
case_mat<-array(dim=c(length(case_files),1,132,132))

#Adding the data to the empty arrays
for (i in 1:nrow(control_mat)) {
  p<-nrow(control_mat[1,1,,])
  control_mat[i,1,,]<-control_data[[i]][,-(133:167)]
}

for (i in 1:nrow(case_mat)) {
  p<-nrow(case_mat[1,1,,])
  case_mat[i,1,,]<-case_data[[i]][,-(133:167)]
}

#####################################
#Zeroing out negative correlations#
#####################################

#Creating a correlation heat map of average case and control correlations
control_mat_zero<-Matlab_fMRI(control_mat)

control_list_zero<-list();
for (i in 1:length(control_files)){
  control_list_zero[[i]]<-control_mat_zero[i,1,,]
}
control_corr_zero<-Reduce("+",control_list_zero)/length(control_list_zero)
colnames(control_corr_zero)<-c(1:132)

melted_contpt_zero<-melt(control_corr_zero)

full_cont_zero_plot<-ggplot(data=melted_contpt_zero,aes(x=Var1,y=Var2,fill=value))+
       geom_tile()+scale_fill_gradient("Corr",low="yellow",high="red",
                                       limits=c(0.0,0.20),breaks=seq(0.0,0.20,by=0.04))+
      labs(x="Regions of Interest", y="Regions of Interest",
           title="Control Subjects Average Functional Connectivity \nCorrelation Heat Map - Zeroing out Negative Correlations")+
      theme(plot.title=element_text(hjust=0.5))

case_mat_zero<-Matlab_fMRI(case_mat)
case_list_zero<-list();
for (i in 1:length(case_files)){
  case_list_zero[[i]]<-case_mat_zero[i,1,,]
}
case_corr_zero<-Reduce("+",case_list_zero)/length(case_list_zero)
colnames(case_corr_zero)<-c(1:132)

melted_casept_zero<-melt(case_corr_zero)
full_case_zero_plot<-ggplot(data=melted_casept_zero,aes(x=Var1,y=Var2,fill=value))+
       geom_tile()+scale_fill_gradient("Corr",low="yellow",high="red",
                                       limits=c(0.0,0.20),breaks=seq(0.0,0.20,by=0.04))+
       labs(x="Regions of Interest", y="Regions of Interest",
            title="Case Subjects Average Functional Connectivity \nCorrelation Heat Map - Zeroing out Negative Correlations")+
       theme(plot.title=element_text(hjust=0.5))

######################################
#Keeping negative correlations as is#
######################################

#Creating a correlation heat map of average case and control correlations
control_mat_neg<-Matlab_fMRI_neg(control_mat)

control_list_neg<-list();
for (i in 1:length(control_files)){
  control_list_neg[[i]]<-control_mat_neg[i,1,,]
}
control_corr_neg<-Reduce("+",control_list_neg)/length(control_list_neg)
colnames(control_corr_neg)<-c(1:132)

melted_contpt_neg<-melt(control_corr_neg)
full_cont_neg_plot<-ggplot(data=melted_contpt_neg,aes(x=Var1,y=Var2,fill=value))+
       geom_tile()+scale_fill_gradient("Corr",low="yellow",high="red",
                                      limits=c(-0.20,0.20),breaks=seq(-0.20,0.20,by=0.05))+
       labs(x="Regions of Interest", y="Regions of Interest",
            title="Control Subjects Average Functional Connectivity \nCorrelation Heat Map - No Change to Negative Correlations")+
       theme(plot.title=element_text(hjust=0.5))

case_mat_neg<-Matlab_fMRI_neg(case_mat)

case_list_neg<-list();
for (i in 1:length(case_files)){
  case_list_neg[[i]]<-case_mat_neg[i,1,,]
}
case_corr_neg<-Reduce("+",case_list_neg)/length(case_list_neg)
colnames(case_corr_neg)<-c(1:132)

melted_casept_neg<-melt(case_corr_neg)
full_case_neg_plot<-ggplot(data=melted_casept_neg,aes(x=Var1,y=Var2,fill=value))+
       geom_tile()+scale_fill_gradient("Corr",low="yellow",high="red",
                                       limits=c(-0.20,0.20),breaks=seq(-0.20,0.20,by=0.05))+
       labs(x="Regions of Interest", y="Regions of Interest",
            title="Case Subjects Average Functional Connectivity \nCorrelation Heat Map - No Change to Negative Correlations")+
       theme(plot.title=element_text(hjust=0.5))

grid.arrange(full_case_zero_plot,full_cont_zero_plot,full_case_neg_plot,full_cont_neg_plot)

################################################################
######### COBRE Dataset: Paranoid Schizophrenia only ##########
################################################################

#Bringing in the correlation matrix files for controls (.mat file format)
#NOTE: You may have to change the working directory for the location of your .mat file
#      adjacency matrices. This is assuming pre-processing occured in the CONN toolbox of MatLab.
control_pfiles<-list.files(pattern="*.mat")
control_pdata<-lapply(control_pfiles,function(X) readMat(X)$Z)
control_pmat<-array(dim=c(length(control_pfiles),1,132,132))

#Bringing in the correlation matrix files for patients (.mat file format)
#NOTE: You may have to change the working directory for the location of your .mat file
#      adjacency matrices. This is assuming pre-processing occured in the CONN toolbox of MatLab.
case_pfiles<-list.files(pattern="*.mat")
case_pdata<-lapply(case_pfiles,function(X) readMat(X)$Z)
case_pmat<-array(dim=c(length(case_pfiles),1,132,132))

#Adding the data to the empty arrays
for (i in 1:nrow(control_pmat)) {
  p<-nrow(control_pmat[1,1,,])
  control_pmat[i,1,,]<-control_pdata[[i]][,-(133:167)]
}

for (i in 1:nrow(case_pmat)) {
  p<-nrow(case_pmat[1,1,,])
  case_pmat[i,1,,]<-case_pdata[[i]][,-(133:167)]
}

#####################################
#Zeroing out negative correlations#
#####################################

#Creating a correlation heat map of average case and control correlations
control_mat_pzero<-Matlab_fMRI(control_pmat)

control_list_pzero<-list();
for (i in 1:length(control_pfiles)){
  control_list_pzero[[i]]<-control_mat_pzero[i,1,,]
}
control_corr_pzero<-Reduce("+",control_list_pzero)/length(control_list_pzero)
colnames(control_corr_pzero)<-c(1:132)

melted_contpt_pzero<-melt(control_corr_pzero)

para_cont_zero_plot<-ggplot(data=melted_contpt_pzero,aes(x=Var1,y=Var2,fill=value))+
  geom_tile()+scale_fill_gradient("Corr",low="yellow",high="red",
                                  limits=c(0.0,0.20),breaks=seq(0.0,0.20,by=0.04))+
  labs(x="Regions of Interest", y="Regions of Interest",
       title="Paranoid Schizophrenia Control Subjects Average Functional Connectivity \nCorrelation Heat Map - Zeroing out Negative Correlations")+
  theme(plot.title=element_text(hjust=0.5))

case_mat_pzero<-Matlab_fMRI(case_pmat)
case_list_pzero<-list();
for (i in 1:length(case_pfiles)){
  case_list_pzero[[i]]<-case_mat_pzero[i,1,,]
}
case_corr_pzero<-Reduce("+",case_list_pzero)/length(case_list_pzero)
colnames(case_corr_pzero)<-c(1:132)

melted_casept_pzero<-melt(case_corr_pzero)
para_case_zero_plot<-ggplot(data=melted_casept_pzero,aes(x=Var1,y=Var2,fill=value))+
  geom_tile()+scale_fill_gradient("Corr",low="yellow",high="red",
                                  limits=c(0.0,0.20),breaks=seq(0.0,0.20,by=0.04))+
  labs(x="Regions of Interest", y="Regions of Interest",
       title="Paranoid Schizophrenia Case Subjects Average Functional Connectivity \nCorrelation Heat Map - Zeroing out Negative Correlations")+
  theme(plot.title=element_text(hjust=0.5))

######################################
#Keeping negative correlations as is#
######################################

#Creating a correlation heat map of average case and control correlations
control_mat_pneg<-Matlab_fMRI_neg(control_pmat)

control_list_pneg<-list();
for (i in 1:length(control_pfiles)){
  control_list_pneg[[i]]<-control_mat_pneg[i,1,,]
}
control_corr_pneg<-Reduce("+",control_list_pneg)/length(control_list_pneg)
colnames(control_corr_pneg)<-c(1:132)

melted_contpt_pneg<-melt(control_corr_pneg)
para_cont_neg_plot<-ggplot(data=melted_contpt_pneg,aes(x=Var1,y=Var2,fill=value))+
  geom_tile()+scale_fill_gradient("Corr",low="yellow",high="red",
                                  limits=c(-0.20,0.20),breaks=seq(-0.20,0.20,by=0.05))+
  labs(x="Regions of Interest", y="Regions of Interest",
       title="Paranoid Schizophrenia Control Subjects Average Functional Connectivity \nCorrelation Heat Map - No Change to Negative Correlations")+
  theme(plot.title=element_text(hjust=0.5))

case_mat_pneg<-Matlab_fMRI_neg(case_pmat)

case_list_pneg<-list();
for (i in 1:length(case_pfiles)){
  case_list_pneg[[i]]<-case_mat_pneg[i,1,,]
}
case_corr_pneg<-Reduce("+",case_list_pneg)/length(case_list_pneg)
colnames(case_corr_pneg)<-c(1:132)

melted_casept_pneg<-melt(case_corr_pneg)
para_case_neg_plot<-ggplot(data=melted_casept_pneg,aes(x=Var1,y=Var2,fill=value))+
  geom_tile()+scale_fill_gradient("Corr",low="yellow",high="red",
                                  limits=c(-0.20,0.20),breaks=seq(-0.20,0.20,by=0.05))+
  labs(x="Regions of Interest", y="Regions of Interest",
       title="Paranoid Schizophrenia Case Subjects Average Functional Connectivity \nCorrelation Heat Map - No Change to Negative Correlations")+
  theme(plot.title=element_text(hjust=0.5))

grid.arrange(para_case_zero_plot,para_cont_zero_plot,para_case_neg_plot,para_cont_neg_plot)