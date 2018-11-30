#-------Functions to be used in full simulations-------#
#Zeroing out the diagonal and zeroing out any negative correlation
fMRI_mimic<-function(data){
  q<-nrow(data)
  p<-nrow(data[1,1,,])
  data2<-data
  for (k in 1:q){
    for (i in 1:p){
      data2[k,1,i,i]<-0
    }
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        ifelse(data2[k,1,i,j]<0,data2[k,1,i,j]<-0,data2[k,1,i,j])
        ifelse(data2[k,1,j,i]<0,data2[k,1,j,i]<-0,data2[k,1,j,i])
      }
    }
  }
  return(data2)
}

#Zeroing out the diagonal and zeroing out any negative correlation - only when dealing with a matrix!
fMRI_mimic2<-function(data){
  q<-nrow(data)
  data2<-data
  for (k in 1:q){
    data2[k,k]<-0
  }
  for (i in 1:(q-1)){
    for (j in (i+1):q){
      ifelse(data2[i,j]<0,data2[i,j]<-0,data2[i,j])
      ifelse(data2[j,i]<0,data2[j,i]<-0,data2[j,i])
    }
  }
  return(data2)
}

#Calculating the degree matrix for each simulation
degree<-function(data){
  q<-nrow(data)
  p<-nrow(data[1,1,,])
  data_sub<-data
  for (j in 1:q){
    for (i in 1:p){
      data_sub[j,2,i,i]<-colSums(data[j,1,,])[i]
    }
  }
  return(data_sub)
}

#Calculating the effective resistance matrix for each simulation
resist<-function(data) {
  p<-nrow(data)
  q<-nrow(data[1,1,,])
  eff_res<-array(dim=c(p,1,q,q))
  ones<-matrix(rep(1,len=q),nrow=1,ncol=q)
  for (i in 1:p) {
    eff_res[i,1,,]<-diag(ginv(data[i,2,,]-data[i,1,,]))%*%ones+t(ones)%*%t(diag(ginv(data[i,2,,]-data[i,1,,])))-2*ginv(data[i,2,,]-data[i,1,,])
  }
  return(eff_res)
}

#Calculating the resistance pertubation distance matrix
distance<-function(eff_res) {
  p<-nrow(eff_res)
  RPD_matrix<-matrix(NA,p,p)
  for(i in 1:p){
    RPD_matrix[i,i]=0
  }
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      RPD_matrix[i,j]=RPD_matrix[j,i]=sqrt(sum(abs(eff_res[i,1,,]-eff_res[j,1,,])^2))
    }
  }
  return(RPD_matrix)
}

#Kernel function
kernel<-function(x,rho){   
  K=exp(-x^2/rho)
  return(K)
}

#Inverse logit function
invlogit<-function(x){
  return(exp(x)/(exp(x)+1))
}

#Score function
sim_score<-function(x,y,lower,upper,m){  
  n<-ncol(x)
  fit<-glm(y~1,family=binomial)  ## fits the null model
  betahat<-fit$coef  ##pulls the coefficient from the null model (intercept)
  mu_0<-invlogit(betahat)   ##inverse logit - gets back the proportion of ones in 0/1 y variable
  D_0<-diag(mu_0*(1-mu_0),n)  ##diagonal matrix of size nxn where diagonal entries are p(1-p)
  X<-rep(1,n)  ##row vector of ones of length n
  P_0<-D_0-D_0%*%X%*%solve(t(X)%*%D_0%*%X)%*%t(X)%*%D_0
  rho<-seq(lower,upper,length=m)  ##grid of rho values for kernel
  S<-rep(0,m)  ##row vector of zeros of grid length
  for(l in 1:m){  ##for each of the grid values, calculates the kernel
    
    k<-kernel(x,rho[l])
    
    Q<-t(y-mu_0)%*%k%*%(y-mu_0)  ##test statistic
    mu<-tr(P_0%*%k)    
    sigma<-sqrt(2*tr(P_0%*%k%*%P_0%*%k))
    S[l]<-(Q-mu)/sigma  ##standardized version of test statistic for each value of kernel
  }
  M<-max(S)  ##finds the max value of the standardized test statistic
  W=0
  for(i in 1:(m-1)){
    W<-W+abs(S[i+1]-S[i])  ##total variation of S in grid
  }
  pvalue<-pnorm(-M)+W*exp(-M^2/2)/sqrt(8*pi) ##(12) in Liu et al. (2008)
  return(pvalue)
}

#Function to find upper and lower bounds of grid search
up_low<-function(x){
  p<-nrow(x)
  q<-ncol(x)
  sq_diff<-diag(p)
  for (i in 1:p){
    for (j in 1:p){
      sum<-0
      for (k in 1:q){
        sum=sum+((x[i,k]-x[j,k])**2)
      }
      sq_diff[i,j]<-sum
    }
  }
  return(sq_diff)
}