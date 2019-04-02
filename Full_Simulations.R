#-------Libraries to be used in simulations and setting the seed-------#
library(MNS)
library(abind)
library(psych)
library(gridExtra)
library(MASS)
library(ggplot2)
library(reshape2)
library(tidyverse)

set.seed(80045)

#Setting working directory (optional) and bringing in functions
source("Simulation_Functions.R")

#----------Example plot of MNS package---------#
set.seed(80045)
N<-3
Nets<-gen.Network(method="cohort",p=15,Nsub=N,sparsity=0.2,                    
                  REsize=10,REprob=0.5,REnoise=1)

plot(Nets,view="sub")

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

edge_importance<-function(nperm,nodes,spars,size,prob,noise){
  
  RPD_offdiag_discon<-c()
  RPD_offdiag_random<-c()
  
  for (a in 1:nperm){
    Net_P1<-gen.Network(method="cohort",p=nodes,Nsub=1,sparsity=spars,REsize=size,REprob=prob,REnoise=noise)
    Matrix_Orig<-Net_P1$Networks[[1]]
    Matrix_Orig<-fMRI_mimic2(Matrix_Orig)
    
    Matrix_Targeted<-disconnected(Matrix_Orig)
    Matrix_Diff<-Matrix_Orig-Matrix_Targeted
    
    num_edges<-sum(Matrix_Diff!=0)
    r_edges<-sample(1:(nodes^2),num_edges)
    
    Vectorize<-c(Matrix_Orig)
    
    for (b in 1:length(r_edges)){
      Vectorize[r_edges[b]]<-0
    }
    
    Matrix_Random<-matrix(Vectorize,nrow=nodes,ncol=nodes)
    
    sim_array<-array(dim=c(3,2,nodes,nodes))
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

RPD_edgeimport<-edge_importance(nperm=1000,nodes=10,sparse=0.75,size=2,prob=0.40,noise=2)
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

grid.arrange(edge_plot)

##Principle 2 (Weight Awareness): in weighted graphs, the larger the weight of the removed edge,
##the greater the impact on the distance should be
weight_aware<-function(nperm,nodes,spars,size,prob,noise){
  
  RPD_offdiag_min<-c()
  RPD_offdiag_max<-c()
  
  for (a in 1:nperm){
    Net_P2<-gen.Network(method="cohort",p=nodes,Nsub=1,sparsity=spars,REsize=size,REprob=prob,REnoise=noise)
    Matrix_Orig<-Net_P2$Networks[[1]]
    Matrix_Orig<-fMRI_mimic2(Matrix_Orig)
    
    Vectorize<-c(Matrix_Orig)
    max<-which(Vectorize==max(Vectorize))
    min<-which(Vectorize==min(Vectorize[Vectorize>0]))
    Vectorize_min<-Vectorize
    Vectorize_min[min]<-0
    Vectorize_max<-Vectorize
    Vectorize_max[max]<-0
    
    Matrix_Max<-matrix(Vectorize_max,nrow=nodes,ncol=nodes)
    Matrix_Min<-matrix(Vectorize_min,nrow=nodes,ncol=nodes)
    
    sim_array<-array(dim=c(3,2,nodes,nodes))
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

RPD_weight<-weight_aware(nperm=1000,nodes=10,sparse=0.55,size=2,prob=0.40,noise=2)
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

grid.arrange(weight_plot)

##Principle 3 (Edge-"Submodularity"): a specific change is more important in a graph with few edges
##than in a much denser, but equally sized graph
submod<-function(nperm,nodes,sparse,size,prob,noise){
  
  RPD_offdiag<-c()
  
  for (a in 1:nperm){
    Net_P3<-gen.Network(method="cohort",p=nodes,Nsub=1,sparsity=sparse,REsize=size,REprob=prob,REnoise=noise)
    Matrix_Orig<-Net_P3$Networks[[1]]
    Matrix_Orig<-fMRI_mimic2(Matrix_Orig)
    
    Vectorize<-c(Matrix_Orig)
    max<-which(Vectorize==max(Vectorize))
    Vectorize[max]<-0
    
    Matrix_Mod<-matrix(Vectorize,nrow=nodes,ncol=nodes)
    
    sim_array<-array(dim=c(2,2,nodes,nodes))
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

RPD_sparse<-as.data.frame(submod(nperm=1000,nodes=10,sparse=0.45,size=2,prob=0.40,noise=2))
RPD_sparse$Iteration<-1:nrow(RPD_sparse)
names(RPD_sparse)<-c("RPD_sparse","Iteration")

RPD_dense<-as.data.frame(submod(nperm=1000,nodes=10,sparse=0.95,size=2,prob=0.40,noise=2))
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

grid.arrange(density_plot)

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

focus<-function(nperm,nodes,sparse,size,prob,noise){
  
  RPD_offdiag_target<-c()
  RPD_offdiag_random<-c()
  
  for (a in 1:nperm){
    Net_P2<-gen.Network(method="cohort",p=nodes,Nsub=1,sparsity=sparse,REsize=size,REprob=prob,REnoise=noise)
    Matrix_Orig<-Net_P2$Networks[[1]]
    Matrix_Orig<-fMRI_mimic2(Matrix_Orig)
    
    r_num<-sample(1:nodes,1)
    Matrix_Targeted<-targeted_deletion(Matrix_Orig,r_num)
    
    num_edges<-colSums(Matrix_Orig!=0)[r_num]
    r_edges<-sample(1:(nodes^2),num_edges)
    
    Vectorize<-c(Matrix_Orig)
    
    for (b in 1:length(r_edges)){
      Vectorize[r_edges[b]]<-0
    }

    Matrix_Random<-matrix(Vectorize,nrow=nodes,ncol=nodes)

    sim_array<-array(dim=c(3,2,nodes,nodes))
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

RPD_focus<-focus(nperm=1000,nodes=10,sparse=0.75,size=2,prob=0.40,noise=2)
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

grid.arrange(focus_plot)


#-------Two generation process full simulation: Power-------#
power_sim<-function(N1,N2,N3,x,z,nperm,nodes,sparse,size,prob,noise){
  
  y<-c(rep(0,(N1+x)),rep(1,(z+N3)))
  upper_tails<-c()
  tail_prob<-0
  
  for (a in 1:nperm){
    Net1<-gen.Network(method="cohort", p=nodes, Nsub=N1, sparsity=sparse, REsize=size, REprob=prob, REnoise=noise)
    Net2<-gen.Network(method="cohort", p=nodes, Nsub=N2, sparsity=sparse, REsize=size, REprob=prob, REnoise=noise)
    Net3<-gen.Network(method="cohort", p=nodes, Nsub=N3, sparsity=sparse, REsize=size, REprob=prob, REnoise=noise)
    
    sim_array1a<-array(dim=c(N1,2,nodes,nodes))
    sim_array1b<-array(dim=c(x,2,nodes,nodes)) 
    sim_array2a<-array(dim=c(z,2,nodes,nodes)) 
    sim_array2b<-array(dim=c(nodes))
    
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

power_sim(N1=10,N2=80,N3=10,x=35,z=45,nperm=100,nodes=90,sparse=0.75,size=10,prob=0.65,noise=3)

#-------One generation process full simulation: Type I error-------#
type1err_sim<-function(N,bin_prob,nperm,nodes,sparse,size,prob,noise){
  
  z<-rbinom(N,1,bin_prob)
  upper_tails<-c()
  tail_prob<-0
  
  for (a in 1:nperm){
    Net<-gen.Network(method="cohort",p=nodes,Nsub=N,sparsity=sparse,REsize=size,REprob=prob,REnoise=noise)
    
    sim_array<-array(dim=c(N,2,nodes,nodes))

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

type1err_sim(N=100,bin_prob=0.45,nperm=10000,nodes=90,sparse=0.75,size=10,prob=0.65,noise=3)

#-------One generation process full simulation: Control/Patient Splits-------#
prob_sim<-function(N,lower,upper,inc,nodes,sparse,size,prob,noise){
  
  grid<-seq(lower,upper,by=inc)
  upper_tails<-c()
  
  for (k in 1:length(grid)){
    Net<-gen.Network(method="cohort",p=nodes,Nsub=N,sparsity=sparse,REsize=size,REprob=prob,REnoise=noise)
    
    sim_array<-array(dim=c(N,2,nodes,nodes))
    
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

prob_sim<-prob_sim(N=100,lower=0.05,upper=0.95,inc=0.01,nodes=90,sparse=0.75,size=10,prob=0.65,noise=3)

plot_col<-rgb(89,145,142,names=NULL,maxColorValue=255)
plot(x=seq(0.05,0.95,by=0.01),y=prob_sim,xlab="Control/Case Split",ylab="P-Value",main="P-Values of Varying Control/Case Splits under \nOne Generation Simulation Process for N=100",type="p",col=plot_col,pch=16,axes=FALSE)
axis(1,at=seq(0.0,1.0,by=0.10),las=2)
axis(2,at=seq(0.0,0.8,by=0.10),las=1)
abline(h=0.05,col="firebrick4",lwd=3,lty=2)

#-------Two generation process full simulation: Power Noise Splits-------#
noise_sim<-function(N1,N2,N3,lower,upper,inc,nperm,nodes,sparse,size,prob,noise){
  
  grid<-rep(seq(lower,upper,inc),nperm)
  upper_tails<-rep(NA,length(grid))
  
  for (k in 1:length(grid)){
    y<-c(rep(0,(N1+(grid[k]*N2))),rep(1,N3+((1-grid[k])*N2)))
    x<-grid[k]*N2
      
    Net1<-gen.Network(method="cohort",p=nodes,Nsub=N1,sparsity=sparse,REsize=size,REprob=prob,REnoise=noise)
    Net2<-gen.Network(method="cohort",p=nodes,Nsub=N2,sparsity=sparse,REsize=size,REprob=prob,REnoise=noise)
    Net3<-gen.Network(method="cohort",p=nodes,Nsub=N3,sparsity=sparse,REsize=size,REprob=prob,REnoise=noise)
      
    sim_array1a<-array(dim=c(N1,2,nodes,nodes))
    sim_array1b<-array(dim=c((grid[k]*N2),2,nodes,nodes)) 
    sim_array2a<-array(dim=c(((1-grid[k])*N2),2,nodes,nodes)) 
    sim_array2b<-array(dim=c(nodes))
      
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

noise<-noise_sim(N1=10,N2=80,N3=10,lower=0.05,upper=0.95,inc=0.05,nperm=100,
                 nodes=90,sparse=0.75,size=10,prob=0.65,noise=3)

plot_col<-rgb(89,145,142,names=NULL,maxColorValue=255)
plot(x=rep(seq(0.05,0.95,by=0.05),100),y=noise,xlab="Noise Split",ylab="P-Value",main="P-Values of Varying Noise Splits under \nTwo Generation Simulation Process for N=100",type="p",col=plot_col,pch=16,axes=FALSE)
axis(1,at=seq(0.0,1.0,by=0.05),las=2)
axis(2,at=seq(0.0,0.5,by=0.05),las=1)
abline(h=0.05,col="firebrick4",lwd=3,lty=2)

x<-rep(seq(0.05,0.95,by=0.05),100)
noise_data<-data.frame(x,noise)
noise_data$insig<-ifelse(noise_data$noise>=0.05,1,0)
tapply(noise_data$insig,noise_data$x,sum)


#-------One generation process full simulation: Bounds-------#
bounds_one_sim<-function(N,low_bnd,up_bnd,len,prob,nperm,nodes,sparse,size,prob,noise){
  
  low_grid<-rep(low_bnd*10^(0:len),nperm)
  low_grid<-low_grid[order(low_grid)]
  
  up_grid<-rep(up_bnd*10^(0:len),nperm)
  up_grid<-up_grid[order(up_grid)]
  
  z<-rbinom(N,1,prob)
  tails<-rep(NA,(len+1)*nperm)
  
  for (g in 1:length(low_grid)){
    Net<-gen.Network(method="cohort",p=nodes,Nsub=N,sparsity=sparse,REsize=size,REprob=prob,REnoise=noise)
      
    sim_array<-array(dim=c(N,2,nodes,nodes))
      
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

onesim_bounds_pvalues<-bounds_one_sim(N=100,low_bnd=0.00001,up_bnd=0.01,len=10,prob=0.45,nperm=10,
                                      nodes=90,sparse=0.75,size=10,prob=0.65,noise=3)

plot_col<-rgb(89,145,89,names=NULL,maxColorValue=255)
plot(x=rep(seq(0,10,by=1),10),y=onesim_bounds_pvalues,xlab="Iteration (10^i) for i=0,...,10",ylab="P-Value",main="P-Values of Varying Upper and Lower Bounds under \nOne Generation Simulation Process for N=100",type="p",col=plot_col,pch=16)
abline(h=0.05,col="firebrick4",lwd=3,lty=2)


#-------Two generation process full simulation: Bounds-------#
bounds_two_sim<-function(N1,N2,low_bnd,up_bnd,len,nperm,nodes,sparse,size,prob,noise){
  
  low_grid<-rep(low_bnd*10^(0:len),nperm)
  low_grid<-low_grid[order(low_grid)]
  
  up_grid<-rep(up_bnd*10^(0:len),nperm)
  up_grid<-up_grid[order(up_grid)]
  
  y<-c(rep(0,N1),rep(1,N2))
  tails<-rep(NA,(len+1)*nperm)
  
  for (g in 1:length(low_grid)){
    Net1<-gen.Network(method="cohort",p=nodes,Nsub=N1,sparsity=sparse,REsize=size,REprob=prob,REnoise=noise)
    Net2<-gen.Network(method="cohort",p=nodes,Nsub=N2,sparsity=sparse,REsize=size,REprob=prob,REnoise=noise)
    
    sim_array1<-array(dim=c(N1,2,nodes,nodes))
    sim_array2<-array(dim=c(N2,2,nodes,nodes)) 
    
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

twosim_bounds_pvalues<-bounds_two_sim(N1=45,N2=55,low_bnd=0.00001,up_bnd=0.01,len=10,nperm=10,
                                      nodes=90,sparse=0.75,size=10,prob=0.65,noise=3)

plot_col<-rgb(89,145,89,names=NULL,maxColorValue=255)
plot(x=rep(seq(0,10,by=1),10),y=twosim_bounds_pvalues,xlab="Iteration (10^i) for i=0,...,10",
     ylab="P-Value",main="P-Values of Varying Upper and Lower Bounds under \nTwo Generation Simulation Process for N=100",
     type="p",col=plot_col,pch=16,ylim=c(0.0,0.10))
abline(h=0.05,col="firebrick4",lwd=3,lty=2)

par(mfrow=c(2,1))
plot_col<-rgb(89,145,89,names=NULL,maxColorValue=255)
plot(x=rep(seq(0,10,by=1),10),y=onesim_bounds_pvalues,xlab="Iteration (10^i) for i=0,...,10",ylab="P-Value",main="P-Values of Varying Upper and Lower Bounds under \nOne Generation Simulation Process for N=100",type="p",col=plot_col,pch=16)
abline(h=0.05,col="firebrick4",lwd=3,lty=2)

plot_col<-rgb(89,145,89,names=NULL,maxColorValue=255)
plot(x=rep(seq(0,10,by=1),10),y=twosim_bounds_pvalues,xlab="Iteration (10^i) for i=0,...,10",
     ylab="P-Value",main="P-Values of Varying Upper and Lower Bounds under \nTwo Generation Simulation Process for N=100",
     type="p",col=plot_col,pch=16,ylim=c(0.0,0.10))
abline(h=0.05,col="firebrick4",lwd=3,lty=2)