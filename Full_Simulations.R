#-------Libraries to be used in simulations and setting the seed-------#
library(MNS)
library(abind)
library(psych)
library(gridExtra)
library(MASS)
library(ggplot2)
library(reshape2)
library(tidyr)

set.seed(80045)

#Setting working directory and bringing in functions
#setwd("C:/Users/jenseale/Dropbox/MS_Thesis_Work/R_Code/")
setwd("/Users/alexandriajensen/Dropbox/MS_Thesis_Work/R_Code/")
source("Simulation_Functions.R")

#----------Example plot of MNS package---------#
set.seed(80045)
N<-3
Nets<-gen.Network(method="cohort",p=15,Nsub=N,sparsity=0.2,                    
                  REsize=10,REprob=0.5,REnoise=1)

tiff(filename="C:/Users/jenseale/Dropbox/MS_Thesis_Work/PLOS_One_Manuscript/Figures/MNS_Example.tiff",width=9,height=6,units="in",res=300)
plot(Nets,view="sub")
dev.off()

#-------Koutra et al principle simulations-------#
##Principle 1 (Edge Importance): changes that create disconnected components should
##be penalized more than changes that maintain the connectivity properties of the graph
disconnected<-function(data){
  q<-nrow(data)
  data2<-data
  for (i in 1:(q/2)){
    for (j in (q/2):q){
      data2[i,j]<-0
      data2[j,i]<-0
    }
  }
  return(data2)
}

edge_importance<-function(nperm){
  
  RPD_offdiag_discon<-c()
  RPD_offdiag_random<-c()
  
  for (a in 1:nperm){
    Net_P1<-gen.Network(method="cohort",p=10,Nsub=1,sparsity=0.75,REsize=2,REprob=0.40,REnoise=2)
    Matrix_Orig<-Net_P1$Networks[[1]]
    Matrix_Orig<-fMRI_mimic2(Matrix_Orig)
    
    Matrix_Targeted<-disconnected(Matrix_Orig)
    Matrix_Diff<-Matrix_Orig-Matrix_Targeted
    
    num_edges<-sum(Matrix_Diff!=0)
    r_edges<-sample(1:100,num_edges)
    
    Vectorize<-c(Matrix_Orig)
    
    for (b in 1:length(r_edges)){
      Vectorize[r_edges[b]]<-0
    }
    
    Matrix_Random<-matrix(Vectorize,nrow=10,ncol=10)
    
    sim_array<-array(dim=c(3,2,10,10))
    sim_array[1,1,,]<-Matrix_Orig
    sim_array[2,1,,]<-Matrix_Targeted
    sim_array[3,1,,]<-Matrix_Random
    
    for (b in 1:nrow(sim_array)) {
      p<-nrow(sim_array[1,1,,])
      sim_array[b,2,,]<-matrix(0,p,p)
    }
    sim_array<-degree(sim_array)
    eff_res<-resist(sim_array)
    RPD_matrix<-distance(eff_res)
    RPD_offdiag_discon[a]<-RPD_matrix[1,2]
    RPD_offdiag_random[a]<-RPD_matrix[1,3]
  }
  
  RPD_diag<-cbind(RPD_offdiag_discon,RPD_offdiag_random)
  return(RPD_diag)
}

RPD_edgeimport<-edge_importance(1000)
RPD_edgeimport<-cbind(1:nrow(RPD_edgeimport),RPD_edgeimport)
RPD_edgeimport<-as.data.frame(RPD_edgeimport)

summary(RPD_edgeimport[,2])
summary(RPD_edgeimport[,3])

RPD_edge_long<-gather(RPD_edgeimport,paradigm,RPD,RPD_offdiag_discon:RPD_offdiag_random,factor_key=TRUE)
RPD_edge_long$paradigm<-ifelse(RPD_edge_long$paradigm=="RPD_offdiag_discon","Disconnected Components","Random Removal")
edge_plot<-ggplot(RPD_edge_long,aes(x=paradigm,y=log(RPD),fill=paradigm))+geom_boxplot()+
           stat_summary(fun.y=mean,geom="point",shape=23,size=4)+scale_color_brewer(palette="Dark2")+
           labs(title="Edge Awareness 1000 Iteration Simulations",
                x="Paradigm",y="Natural Logarithm of RPD")+guides(fill=FALSE)+
           theme(plot.title=element_text(hjust=0.5))

# edge_discon_plot<-ggplot(RPD_edgeimport,aes(x=RPD_edgeimport[,1],y=RPD_edgeimport[,2]),size=12)+
#                   geom_point(color="#106D0F")+ggtitle("Edge Awareness Simulations: \nEdges Removed Resulting in Disconnected Components")+
#                   labs(x="Iteration",y="Resistance Perturbation Distance")+theme(plot.title=element_text(hjust=0.5))+
#                   geom_hline(yintercept=summary(RPD_edgeimport[,2])[[2]],linetype="dashed",color="red",size=1.25)+
#                   geom_hline(yintercept=summary(RPD_edgeimport[,2])[[5]],linetype="dashed",color="red",size=1.25)
# 
# edge_random_plot<-ggplot(RPD_edgeimport,aes(x=RPD_edgeimport[,1],y=RPD_edgeimport[,3]),size=12)+
#                   geom_point(color="#323472")+ggtitle("Edge Awareness Simulations: \nRandom Removal of Edges")+
#                   labs(x="Iteration",y="Resistance Perturbation Distance")+theme(plot.title=element_text(hjust=0.5))+
#                   geom_hline(yintercept=summary(RPD_edgeimport[,3])[[2]],linetype="dashed",color="red",size=1.25)+
#                   geom_hline(yintercept=summary(RPD_edgeimport[,3])[[5]],linetype="dashed",color="red",size=1.25)

tiff(filename="/Users/alexandriajensen/Dropbox/MS_Thesis_Work/PLOS_One_Manuscript/Figures/Edge_Awareness_Simulplot_25Nov2018.tiff",width=13,height=10,units="in",res=300)
grid.arrange(edge_plot)
dev.off()

##Principle 2 (Weight Awareness): in weighted graphs, the larger the weight of the removed edge,
##the greater the impact on the distance should be
weight_aware<-function(nperm){
  
  RPD_offdiag_min<-c()
  RPD_offdiag_max<-c()
  
  for (a in 1:nperm){
    Net_P2<-gen.Network(method="cohort",p=10,Nsub=1,sparsity=0.55,REsize=2,REprob=0.40,REnoise=2)
    Matrix_Orig<-Net_P2$Networks[[1]]
    Matrix_Orig<-fMRI_mimic2(Matrix_Orig)
    
    Vectorize<-c(Matrix_Orig)
    max<-which(Vectorize==max(Vectorize))
    min<-which(Vectorize==min(Vectorize[Vectorize>0]))
    Vectorize_min<-Vectorize
    Vectorize_min[min]<-0
    Vectorize_max<-Vectorize
    Vectorize_max[max]<-0
    
    Matrix_Max<-matrix(Vectorize_max,nrow=10,ncol=10)
    Matrix_Min<-matrix(Vectorize_min,nrow=10,ncol=10)
    
    sim_array<-array(dim=c(3,2,10,10))
    sim_array[1,1,,]<-Matrix_Orig
    sim_array[2,1,,]<-Matrix_Max
    sim_array[3,1,,]<-Matrix_Min
    
    for (b in 1:nrow(sim_array)) {
      p<-nrow(sim_array[1,1,,])
      sim_array[b,2,,]<-matrix(0,p,p)
    }
    sim_array<-degree(sim_array)
    eff_res<-resist(sim_array)
    RPD_matrix<-distance(eff_res)
    RPD_offdiag_max[a]<-RPD_matrix[1,2]
    RPD_offdiag_min[a]<-RPD_matrix[1,3]
  }
  
  RPD_diag<-cbind(RPD_offdiag_max,RPD_offdiag_min)
  return(RPD_diag)
}

RPD_weight<-weight_aware(1000)
RPD_weight<-cbind(1:nrow(RPD_weight),RPD_weight)
RPD_weight<-as.data.frame(RPD_weight)

RPD_weight_long<-gather(RPD_weight,paradigm,RPD,RPD_offdiag_max:RPD_offdiag_min,factor_key=TRUE)
RPD_weight_long$paradigm<-ifelse(RPD_weight_long$paradigm=="RPD_offdiag_max","Maximum Non-Zero Correlation","Minimum Non-Zero Correlation")
weight_plot<-ggplot(RPD_weight_long,aes(x=paradigm,y=log(RPD),fill=paradigm))+geom_boxplot()+
  stat_summary(fun.y=mean,geom="point",shape=23,size=4)+scale_color_brewer(palette="Dark2")+
  labs(title="Weight Awareness 1000 Iteration Simulations",
       x="Paradigm",y="Natural Logarithm of RPD")+guides(fill=FALSE)+
  theme(plot.title=element_text(hjust=0.5))

summary(RPD_weight[,2])
summary(RPD_weight[,3])
# weight_min_plot<-ggplot(RPD_weight,aes(x=RPD_weight[,1],y=RPD_weight[,3]),size=12)+
#   geom_point(color="#106D0F")+ggtitle("Weight Awareness Simulations: Minimum Correlation")+
#   labs(x="Iteration",y="Resistance Perturbation Distance")+theme(plot.title=element_text(hjust=0.5))+
#   geom_hline(yintercept=summary(RPD_weight[,2])[[2]],linetype="dashed",color="red",size=1.25)+
#   geom_hline(yintercept=summary(RPD_weight[,2])[[5]],linetype="dashed",color="red",size=1.25)

# weight_max_plot<-ggplot(RPD_weight,aes(x=RPD_weight[,1],y=RPD_weight[,2]),size=12)+
#   geom_point(color="#323472")+ggtitle("Weight Awareness Simulations: Maximum Correlation")+
#   labs(x="Iteration",y="Resistance Perturbation Distance")+theme(plot.title=element_text(hjust=0.5))+
#   geom_hline(yintercept=summary(RPD_weight[,3])[[2]],linetype="dashed",color="red",size=1.25)+
#   geom_hline(yintercept=summary(RPD_weight[,3])[[5]],linetype="dashed",color="red",size=1.25)

tiff(filename="/Users/alexandriajensen/Dropbox/MS_Thesis_Work/PLOS_One_Manuscript/Figures/Weight_Aware_Simulplot_25Nov2018.tiff",width=13,height=10,units="in",res=300)
grid.arrange(weight_plot)
dev.off()


##Principle 3 (Edge-"Submodularity"): a specific change is more important in a graph with few edges
##than in a much denser, but equally sized graph
submod<-function(nperm,sparse){
  
  RPD_offdiag<-c()
  
  for (a in 1:nperm){
    Net_P3<-gen.Network(method="cohort",p=10,Nsub=1,sparsity=sparse,REsize=2,REprob=0.40,REnoise=2)
    Matrix_Orig<-Net_P3$Networks[[1]]
    Matrix_Orig<-fMRI_mimic2(Matrix_Orig)
    
    Vectorize<-c(Matrix_Orig)
    max<-which(Vectorize==max(Vectorize))
    Vectorize[max]<-0
    
    Matrix_Mod<-matrix(Vectorize,nrow=10,ncol=10)
    
    sim_array<-array(dim=c(2,2,10,10))
    sim_array[1,1,,]<-Matrix_Orig
    sim_array[2,1,,]<-Matrix_Mod
    
    for (b in 1:nrow(sim_array)) {
      p<-nrow(sim_array[1,1,,])
      sim_array[b,2,,]<-matrix(0,p,p)
    }
    sim_array<-degree(sim_array)
    eff_res<-resist(sim_array)
    RPD_matrix<-distance(eff_res)
    RPD_offdiag[a]<-RPD_matrix[1,2]
  }
  return(RPD_offdiag)
}

RPD_sparse<-as.data.frame(submod(1000,0.45))
RPD_sparse$Iteration<-1:nrow(RPD_sparse)
names(RPD_sparse)<-c("RPD_sparse","Iteration")

RPD_dense<-as.data.frame(submod(1000,0.95))
RPD_dense$Iteration<-1:nrow(RPD_dense)
names(RPD_dense)<-c("RPD_dense","Iteration")

RPD_density<-merge(x=RPD_dense,y=RPD_sparse,by="Iteration",all.x=TRUE)
RPD_density_long<-gather(RPD_density,paradigm,RPD,RPD_dense:RPD_sparse,factor_key=TRUE)
RPD_density_long$paradigm<-ifelse(RPD_density_long$paradigm=="RPD_dense","Dense Network","Sparse Network")
density_plot<-ggplot(RPD_density_long,aes(x=paradigm,y=log(RPD),fill=paradigm))+geom_boxplot()+
  stat_summary(fun.y=mean,geom="point",shape=23,size=4)+scale_color_brewer(palette="Dark2")+
  labs(title="Edge-'Submodularity' 1000 Iteration Simulations",
       x="Paradigm",y="Natural Logarithm of RPD")+guides(fill=FALSE)+
  theme(plot.title=element_text(hjust=0.5))

summary(RPD_sparse$RPD_sparse)
summary(RPD_dense$RPD_dense)

# edgemod_sparse_plot<-ggplot(RPD_sparse,aes(x=Iteration,y=RPD_sparse),size=12)+
#   geom_point(color="#106D0F")+ggtitle("Edge-'Submodularity' Simulations: Sparse Network")+
#   labs(x="Iteration",y="Resistance Perturbation Distance")+theme(plot.title=element_text(hjust=0.5))+
#   geom_hline(yintercept=summary(RPD_sparse[,1])[[2]],linetype="dashed",color="red",size=1.25)+
#   geom_hline(yintercept=summary(RPD_sparse[,1])[[5]],linetype="dashed",color="red",size=1.25)
# 
# edgemod_dense_plot<-ggplot(RPD_dense,aes(x=Iteration,y=RPD_dense),size=12)+
#   geom_point(color="#323472")+ggtitle("Edge-'Submodularity' Simulations: Dense Network")+
#   labs(x="Iteration",y="Resistance Perturbation Distance")+theme(plot.title=element_text(hjust=0.5))+
#   geom_hline(yintercept=summary(RPD_dense[,1])[[2]],linetype="dashed",color="red",size=1.25)+
#   geom_hline(yintercept=summary(RPD_dense[,1])[[5]],linetype="dashed",color="red",size=1.25)
# 

tiff(filename="/Users/alexandriajensen/Dropbox/MS_Thesis_Work/PLOS_One_Manuscript/Figures/Edge_Submod_Simulplot_25Nov2018.tiff",width=13,height=10,units="in",res=300)
grid.arrange(density_plot)
dev.off()

##Principle 4 (Focus Awareness): random changes in graphs are less important than targeted
##changes of the same extent

###NOTE: KOUTRA ET AL TESTED FOR THIS BY EITHER REMOVING ALL EDGES CONNECTED TO A VERTEX (TARGETED)
###      OR RANDOMLY REMOVING THE SAME NUMBER OF EDGES FROM THE WHOLE GRAPH (RANDOM)
targeted_deletion<-function(data,x){
  q<-nrow(data)
  data2<-data
  for (i in 1:q){
    data2[x,i]<-0
    data2[i,x]<-0
  }
  return(data2)
}

focus<-function(nperm){
  
  RPD_offdiag_target<-c()
  RPD_offdiag_random<-c()
  
  for (a in 1:nperm){
    Net_P2<-gen.Network(method="cohort",p=10,Nsub=1,sparsity=0.75,REsize=2,REprob=0.40,REnoise=2)
    Matrix_Orig<-Net_P2$Networks[[1]]
    Matrix_Orig<-fMRI_mimic2(Matrix_Orig)
    
    r_num<-sample(1:10,1)
    Matrix_Targeted<-targeted_deletion(Matrix_Orig,r_num)
    
    num_edges<-colSums(Matrix_Orig!=0)[r_num]
    r_edges<-sample(1:100,num_edges)
    
    Vectorize<-c(Matrix_Orig)
    
    for (b in 1:length(r_edges)){
      Vectorize[r_edges[b]]<-0
    }

    Matrix_Random<-matrix(Vectorize,nrow=10,ncol=10)

    sim_array<-array(dim=c(3,2,10,10))
    sim_array[1,1,,]<-Matrix_Orig
    sim_array[2,1,,]<-Matrix_Targeted
    sim_array[3,1,,]<-Matrix_Random
    
    for (b in 1:nrow(sim_array)) {
      p<-nrow(sim_array[1,1,,])
      sim_array[b,2,,]<-matrix(0,p,p)
    }
    sim_array<-degree(sim_array)
    eff_res<-resist(sim_array)
    RPD_matrix<-distance(eff_res)
    RPD_offdiag_target[a]<-RPD_matrix[1,2]
    RPD_offdiag_random[a]<-RPD_matrix[1,3]
  }
  
  RPD_diag<-cbind(RPD_offdiag_target,RPD_offdiag_random)
  return(RPD_diag)
}

RPD_focus<-focus(1000)
RPD_focus<-cbind(1:nrow(RPD_focus),RPD_focus)
RPD_focus<-as.data.frame(RPD_focus)

RPD_focus_long<-gather(RPD_focus,paradigm,RPD,RPD_offdiag_target:RPD_offdiag_random,factor_key=TRUE)
RPD_focus_long$paradigm<-ifelse(RPD_focus_long$paradigm=="RPD_offdiag_target","Targeted Deletions","Random Deletions")
focus_plot<-ggplot(RPD_focus_long,aes(x=paradigm,y=log(RPD),fill=paradigm))+geom_boxplot()+
  stat_summary(fun.y=mean,geom="point",shape=23,size=4)+scale_color_brewer(palette="Dark2")+
  labs(title="Focus Awareness 1000 Iteration Simulations",
       x="Paradigm",y="Natural Logarithm of RPD")+guides(fill=FALSE)+
  theme(plot.title=element_text(hjust=0.5))

summary(RPD_focus[,2])
summary(RPD_focus[,3])

# focus_discon_plot<-ggplot(RPD_focus,aes(x=RPD_focus[,1],y=RPD_focus[,2]),size=12)+
#   geom_point(color="#106D0F")+ggtitle("Focus Awareness Simulations: Targeted Changes")+
#   labs(x="Iteration",y="Resistance Perturbation Distance")+theme(plot.title=element_text(hjust=0.5))+
#   geom_hline(yintercept=summary(RPD_focus[,2])[[2]],linetype="dashed",color="red",size=1.25)+
#   geom_hline(yintercept=summary(RPD_focus[,2])[[5]],linetype="dashed",color="red",size=1.25)

# focus_random_plot<-ggplot(RPD_edgeimport,aes(x=RPD_edgeimport[,1],y=RPD_edgeimport[,3]),size=12)+
#   geom_point(color="#323472")+ggtitle("Focus Awareness Simulations: Random Changes")+
#   labs(x="Iteration",y="Resistance Perturbation Distance")+theme(plot.title=element_text(hjust=0.5))+
#   geom_hline(yintercept=summary(RPD_focus[,3])[[2]],linetype="dashed",color="red",size=1.25)+
#   geom_hline(yintercept=summary(RPD_focus[,3])[[5]],linetype="dashed",color="red",size=1.25)

tiff(filename="/Users/alexandriajensen/Dropbox/MS_Thesis_Work/PLOS_One_Manuscript/Figures/Focus_Awareness_Simulplot_25Nov2018.tiff",width=13,height=10,units="in",res=300)
grid.arrange(focus_plot)
dev.off()


#-------Two generation process full simulation: Power-------#
power_sim<-function(N1,N2,N3,x,z,nperm){
  
  y<-c(rep(0,(N1+x)),rep(1,(z+N3)))
  upper_tails<-c()
  tail_prob<-0
  
  for (a in 1:nperm){
    Net1<-gen.Network(method="cohort", p=90, Nsub=N1, sparsity=0.75, REsize=10, REprob=0.65, REnoise=3)
    Net2<-gen.Network(method="cohort", p=90, Nsub=N2, sparsity=0.75, REsize=10, REprob=0.65, REnoise=3)
    Net3<-gen.Network(method="cohort", p=90, Nsub=N3, sparsity=0.75, REsize=10, REprob=0.65, REnoise=3)
    
    sim_array1a<-array(dim=c(N1,2,90,90))
    sim_array1b<-array(dim=c(x,2,90,90)) 
    sim_array2a<-array(dim=c(z,2,90,90)) 
    sim_array2b<-array(dim=c(N3,2,90,90))
    
    for (b in 1:nrow(sim_array1a)) {
      p<-nrow(sim_array1a[1,1,,])
      sim_array1a[b,1,,]<-Net1$Networks[[b]]
      sim_array1a[b,2,,]<-matrix(0,p,p)
    }
    for (c in 1:nrow(sim_array1b)) {
      p<-nrow(sim_array1b[1,1,,])
      sim_array1b[c,1,,]<-Net2$Networks[[c]]
      sim_array1b[c,2,,]<-matrix(0,p,p)
    }
    for (d in 1:nrow(sim_array2a)) {
      p<-nrow(sim_array2a[1,1,,])
      sim_array2a[d,1,,]<-Net2$Networks[[x+d]]
      sim_array2a[d,2,,]<-matrix(0,p,p)
    }
    for (e in 1:nrow(sim_array2b)) {
      p<-nrow(sim_array2b[1,1,,])
      sim_array2b[e,1,,]<-Net3$Networks[[e]]
      sim_array2b[e,2,,]<-matrix(0,p,p)
    }
   sim_array<-abind(sim_array1a,sim_array1b,sim_array2a,sim_array2b,along=1) 

    sim_array<-fMRI_mimic(sim_array)
    sim_array<-degree(sim_array)
    eff_res<-resist(sim_array)
    RPD_matrix<-distance(eff_res)
    
    bounds<-up_low(RPD_matrix)
    U<-max(bounds)*100
    L<-min(bounds[bounds>0])*0.1
    upper_tails[a]<-sim_score(RPD_matrix,y,L,U,5000)
  }
  for (j in 1:length(upper_tails)){
    if(upper_tails[j]>0.05)
      tail_prob<-tail_prob+1
  }
  power<-1-(tail_prob/length(upper_tails))
  return(power)
}

power_sim(N1=10,N2=80,N3=10,x=35,z=45,nperm=100)

#-------One generation process full simulation: Type I error-------#
type1err_sim<-function(N,prob,nperm){
  
  z<-rbinom(N,1,prob)
  upper_tails<-c()
  tail_prob<-0
  
  for (a in 1:nperm){
    Net<-gen.Network(method="cohort", p=90, Nsub=N, sparsity=0.75, REsize=10, REprob=0.65, REnoise=3)
    
    sim_array<-array(dim=c(N,2,90,90))

    for (b in 1:nrow(sim_array)) {
      p<-nrow(sim_array[1,1,,])
      sim_array[b,1,,]<-Net$Networks[[b]]
      sim_array[b,2,,]<-matrix(0,p,p)
    }
    sim_array<-fMRI_mimic(sim_array)
    sim_array<-degree(sim_array)
    eff_res<-resist(sim_array)
    RPD_matrix<-distance(eff_res)
    
    bounds<-up_low(RPD_matrix)
    U<-max(bounds)*100
    L<-min(bounds[bounds>0])*0.1
    upper_tails[a]<-sim_score(RPD_matrix,z,L,U,5000)
  }
  for (j in 1:length(upper_tails)){
    if(upper_tails[j]<0.05)
      tail_prob<-tail_prob+1
  }
  type_one<-tail_prob/length(upper_tails)
  return(type_one)
}

type1err_sim(N=100,prob=0.45,nperm=10000)

#-------One generation process full simulation: Control/Patient Splits-------#
prob_sim<-function(N,lower,upper,inc){
  
  grid<-seq(lower,upper,by=inc)
  upper_tails<-c()
  
  for (k in 1:length(grid)){
    Net<-gen.Network(method="cohort", p=90, Nsub=N, sparsity=0.75, REsize=10, REprob=0.65, REnoise=3)
    
    sim_array<-array(dim=c(N,2,90,90))
    
    for (j in 1:nrow(sim_array)) {
      p<-nrow(sim_array[1,1,,])
      sim_array[j,1,,]<-Net$Networks[[j]]
      sim_array[j,2,,]<-matrix(0,p,p)
    }
    sim_array<-fMRI_mimic(sim_array)
    sim_array<-degree(sim_array)
    eff_res<-resist(sim_array)
    RPD_matrix<-distance(eff_res)
    
    bounds<-up_low(RPD_matrix)
    U<-max(bounds)*100
    L<-min(bounds[bounds>0])*0.1
    
    z<-rbinom(N,1,grid[k])
    upper_tails[k]<-sim_score(RPD_matrix,z,L,U,5000)
  }
  return(upper_tails)
}

prob_sim<-prob_sim(N=100,lower=0.05,upper=0.95,inc=0.01)

tiff(filename="C:/Users/jenseale/Dropbox/MS_Thesis_Work/PLOS_One_Manuscript/Figures/Control_Case_OneGen_Simulplot.tiff",width=13,height=10,units="in",res=300)
plot_col<-rgb(89,145,142,names=NULL,maxColorValue=255)
plot(x=seq(0.05,0.95,by=0.01),y=prob_sim,xlab="Control/Case Split",ylab="P-Value",main="P-Values of Varying Control/Case Splits under \nOne Generation Simulation Process for N=100",type="p",col=plot_col,pch=16,axes=FALSE)
axis(1,at=seq(0.0,1.0,by=0.10),las=2)
axis(2,at=seq(0.0,0.8,by=0.10),las=1)
abline(h=0.05,col="firebrick4",lwd=3,lty=2)
dev.off()

#-------Two generation process full simulation: Power Noise Splits-------#
noise_sim<-function(N1,N2,N3,lower,upper,inc,nperm){
  
  grid<-rep(seq(lower,upper,inc),nperm)
  upper_tails<-rep(NA,length(grid))
  
  for (k in 1:length(grid)){
    y<-c(rep(0,(N1+(grid[k]*N2))),rep(1,N3+((1-grid[k])*N2)))
    x<-grid[k]*N2
      
    Net1<-gen.Network(method="cohort", p=90, Nsub=N1, sparsity=0.75, REsize=10, REprob=0.65, REnoise=3)
    Net2<-gen.Network(method="cohort", p=90, Nsub=N2, sparsity=0.75, REsize=10, REprob=0.65, REnoise=3)
    Net3<-gen.Network(method="cohort", p=90, Nsub=N3, sparsity=0.75, REsize=10, REprob=0.65, REnoise=3)
      
    sim_array1a<-array(dim=c(N1,2,90,90))
    sim_array1b<-array(dim=c((grid[k]*N2),2,90,90)) 
    sim_array2a<-array(dim=c(((1-grid[k])*N2),2,90,90)) 
    sim_array2b<-array(dim=c(N3,2,90,90))
      
    for (b in 1:nrow(sim_array1a)) {
      p<-nrow(sim_array1a[1,1,,])
      sim_array1a[b,1,,]<-Net1$Networks[[b]]
      sim_array1a[b,2,,]<-matrix(0,p,p)
    }
    for (c in 1:nrow(sim_array1b)) {
      p<-nrow(sim_array1b[1,1,,])
      sim_array1b[c,1,,]<-Net2$Networks[[c]]
      sim_array1b[c,2,,]<-matrix(0,p,p)
    }
    for (d in 1:nrow(sim_array2a)) {
      p<-nrow(sim_array2a[1,1,,])
      sim_array2a[d,1,,]<-Net2$Networks[[x+d]]
      sim_array2a[d,2,,]<-matrix(0,p,p)
    }
    for (e in 1:nrow(sim_array2b)) {
      p<-nrow(sim_array2b[1,1,,])
      sim_array2b[e,1,,]<-Net3$Networks[[e]]
      sim_array2b[e,2,,]<-matrix(0,p,p)
    }
    sim_array<-abind(sim_array1a,sim_array1b,sim_array2a,sim_array2b,along=1) 
      
    sim_array<-fMRI_mimic(sim_array)
    sim_array<-degree(sim_array)
    eff_res<-resist(sim_array)
    RPD_matrix<-distance(eff_res)
      
    bounds<-up_low(RPD_matrix)
    U<-max(bounds)*100
    L<-min(bounds[bounds>0])*0.1
    upper_tails[k]<-sim_score(RPD_matrix,y,L,U,5000)
  }
  return(upper_tails)
}

noise<-noise_sim(N1=10,N2=80,N3=10,lower=0.05,upper=0.95,inc=0.05,nperm=100)

tiff(filename="C:/Users/jenseale/Dropbox/MS_Thesis_Work/PLOS_One_Manuscript/Figures/Control_TwoGen_OneGen_Simulplot.tiff",width=13,height=10,units="in",res=300)
plot_col<-rgb(89,145,142,names=NULL,maxColorValue=255)
plot(x=rep(seq(0.05,0.95,by=0.05),100),y=noise,xlab="Noise Split",ylab="P-Value",main="P-Values of Varying Noise Splits under \nTwo Generation Simulation Process for N=100",type="p",col=plot_col,pch=16,axes=FALSE)
axis(1,at=seq(0.0,1.0,by=0.05),las=2)
axis(2,at=seq(0.0,0.5,by=0.05),las=1)
abline(h=0.05,col="firebrick4",lwd=3,lty=2)
dev.off()

x<-rep(seq(0.05,0.95,by=0.05),100)
noise_data<-data.frame(x,noise)
noise_data$insig<-ifelse(noise_data$noise>=0.05,1,0)
tapply(noise_data$insig,noise_data$x,sum)

write.table(noise,"C:/Users/jenseale/Dropbox/MS_Thesis_Work/noise_splits.txt",sep="\t")
#0.05  0.1 0.15  0.2 0.25  0.3 0.35  0.4 0.45  0.5 0.55  0.6 0.65  0.7 0.75  0.8 0.85  0.9 0.95 
#   0    0    0    0    0    2    4   12    4   14   11    4    6    0    0    0    0    0    0 


#-------One generation process full simulation: Bounds-------#
bounds_one_sim<-function(N,low_bnd,up_bnd,len,prob,nperm){
  
  low_grid<-rep(low_bnd*10^(0:len),nperm)
  low_grid<-low_grid[order(low_grid)]
  
  up_grid<-rep(up_bnd*10^(0:len),nperm)
  up_grid<-up_grid[order(up_grid)]
  
  z<-rbinom(N,1,prob)
  tails<-rep(NA,(len+1)*nperm)
  
  for (g in 1:length(low_grid)){
    Net<-gen.Network(method="cohort", p=90, Nsub=N, sparsity=0.75, REsize=10, REprob=0.65, REnoise=3)
      
    sim_array<-array(dim=c(N,2,90,90))
      
     for (h in 1:nrow(sim_array)) {
       p<-nrow(sim_array[1,1,,])
       sim_array[h,1,,]<-Net$Networks[[h]]
       sim_array[h,2,,]<-matrix(0,p,p)
     }
     sim_array<-fMRI_mimic(sim_array)
     sim_array<-degree(sim_array)
     eff_res<-resist(sim_array)
     RPD_matrix<-distance(eff_res)
      
     bounds<-up_low(RPD_matrix)
     U<-max(bounds)*up_grid[g]
     L<-min(bounds[bounds>0])*low_grid[g]
     tails[g]<-sim_score(RPD_matrix,z,L,U,5000)
    }
    
  return(tails)
}

onesim_bounds_pvalues<-bounds_one_sim(N=100,low_bnd=0.00001,up_bnd=0.01,len=10,prob=0.45,nperm=10)

tiff(filename="C:/Users/jenseale/Dropbox/MS_Thesis_Work/PLOS_One_Manuscript/Figures/Bounds_OneGen_Simulplot.tiff",width=13,height=10,units="in",res=300)
plot_col<-rgb(89,145,89,names=NULL,maxColorValue=255)
plot(x=rep(seq(0,10,by=1),10),y=onesim_bounds_pvalues,xlab="Iteration (10^i) for i=0,...,10",ylab="P-Value",main="P-Values of Varying Upper and Lower Bounds under \nOne Generation Simulation Process for N=100",type="p",col=plot_col,pch=16)
abline(h=0.05,col="firebrick4",lwd=3,lty=2)
dev.off()


#-------Two generation process full simulation: Bounds-------#
bounds_two_sim<-function(N1,N2,low_bnd,up_bnd,len,nperm){
  
  low_grid<-rep(low_bnd*10^(0:len),nperm)
  low_grid<-low_grid[order(low_grid)]
  
  up_grid<-rep(up_bnd*10^(0:len),nperm)
  up_grid<-up_grid[order(up_grid)]
  
  y<-c(rep(0,N1),rep(1,N2))
  tails<-rep(NA,(len+1)*nperm)
  
  for (g in 1:length(low_grid)){
    Net1<-gen.Network(method="cohort", p=90, Nsub=N1, sparsity=0.75, REsize=10, REprob=0.65, REnoise=3)
    Net2<-gen.Network(method="cohort",p=90, Nsub=N2, sparsity=0.75, REsize=10, REprob=0.65, REnoise=3)
    
    sim_array1<-array(dim=c(N1,2,90,90))
    sim_array2<-array(dim=c(N2,2,90,90)) 
    
    for (b in 1:nrow(sim_array1)) {
      p<-nrow(sim_array1[1,1,,])
      sim_array1[b,1,,]<-Net1$Networks[[b]]
      sim_array1[b,2,,]<-matrix(0,p,p)
    }
    for (c in 1:nrow(sim_array2)) {
      p<-nrow(sim_array2[1,1,,])
      sim_array2[c,1,,]<-Net2$Networks[[c]]
      sim_array2[c,2,,]<-matrix(0,p,p)
    }
    
    sim_array<-abind(sim_array1,sim_array2,along=1)
    
    sim_array<-fMRI_mimic(sim_array)
    sim_array<-degree(sim_array)
    eff_res<-resist(sim_array)
    RPD_matrix<-distance(eff_res)
    
    bounds<-up_low(RPD_matrix)
    U<-max(bounds)*up_grid[g]
    L<-min(bounds[bounds>0])*low_grid[g]
    tails[g]<-sim_score(RPD_matrix,y,L,U,5000)
  }
  
  return(tails)
}

twosim_bounds_pvalues<-bounds_two_sim(N1=45,N2=55,low_bnd=0.00001,up_bnd=0.01,len=10,nperm=10)

tiff(filename="C:/Users/jenseale/Dropbox/MS_Thesis_Work/PLOS_One_Manuscript/Figures/Bounds_TwoGen_Simulplot.tiff",width=13,height=10,units="in",res=300)
plot_col<-rgb(89,145,89,names=NULL,maxColorValue=255)
plot(x=rep(seq(0,10,by=1),10),y=twosim_bounds_pvalues,xlab="Iteration (10^i) for i=0,...,10",
     ylab="P-Value",main="P-Values of Varying Upper and Lower Bounds under \nTwo Generation Simulation Process for N=100",
     type="p",col=plot_col,pch=16,ylim=c(0.0,0.10))
abline(h=0.05,col="firebrick4",lwd=3,lty=2)
dev.off()


tiff(filename="C:/Users/jenseale/Dropbox/MS_Thesis_Work/PLOS_One_Manuscript/Figures/Bounds_Simulplot.tiff",width=13,height=10,units="in",res=300)
par(mfrow=c(2,1))
plot_col<-rgb(89,145,89,names=NULL,maxColorValue=255)
plot(x=rep(seq(0,10,by=1),10),y=onesim_bounds_pvalues,xlab="Iteration (10^i) for i=0,...,10",ylab="P-Value",main="P-Values of Varying Upper and Lower Bounds under \nOne Generation Simulation Process for N=100",type="p",col=plot_col,pch=16)
abline(h=0.05,col="firebrick4",lwd=3,lty=2)

plot_col<-rgb(89,145,89,names=NULL,maxColorValue=255)
plot(x=rep(seq(0,10,by=1),10),y=twosim_bounds_pvalues,xlab="Iteration (10^i) for i=0,...,10",
     ylab="P-Value",main="P-Values of Varying Upper and Lower Bounds under \nTwo Generation Simulation Process for N=100",
     type="p",col=plot_col,pch=16,ylim=c(0.0,0.10))
abline(h=0.05,col="firebrick4",lwd=3,lty=2)
dev.off()

