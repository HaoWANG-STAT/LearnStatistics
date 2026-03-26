
#R code used for the co-primary EPs - 
#  Used to adjust the original sample size calculation to ensure an Overall Power >80% for both co-primary tests.


library(stats)
library(MASS)
##Co_Primary_Endpoints

CV_AUC<-0.398
CV_Cth<-0.568 # original CV assumption: 0.55

t_Cth_SC<-log(0.95)#underlying true GMR
t_AUC_SC<-log(0.95)

Corr<-0.787 #correlation between the two endpoints

n_SC<-180 #original PK evaluable 261: 174
n_IV<- 90 #original PK evaluable 261: 87

split_alpha<-c(0.025,0.025)


var_AUC<-log(CV_AUC^2 + 1)
var_Cth<-log(CV_Cth^2+1)
Cov<-Corr*sqrt(var_AUC)*sqrt(var_Cth)

sim_trial<-function(){
  rand<-mvrnorm(n = n_SC+n_IV,mu=c(0,0),Sigma=cbind(c(var_Cth,Cov),c(Cov,var_AUC)))
  SC<-rbind(c(t_Cth_SC,t_AUC_SC))[rep(1,n_SC), ]+rand[1:n_SC,]
  IV<-rand[(n_SC+1):(n_SC+n_IV),]
  GMRobs<-colMeans(SC)-colMeans(IV)
  
  t_critic<-qt(split_alpha,n_SC+n_IV-2)
  lowCI<-GMRobs+sqrt( 2*(colSums((SC-colMeans(SC))^2)+colSums((IV-colMeans(IV))^2))/(n_SC+n_IV-2) /(2/(1/n_SC+1/n_IV)))*t_critic
  
  t_critic_recycl<-qt(c(0.05,0.05),n_SC+n_IV-2)
  lowCI_recycl<-GMRobs+sqrt( 2*(colSums((SC-colMeans(SC))^2)+colSums((IV-colMeans(IV))^2))/(n_SC+n_IV-2) /(2/(1/n_SC+1/n_IV)))*t_critic_recycl
  
  return(rbind(exp(GMRobs),exp(lowCI),exp(lowCI_recycl)))
}


n_sim<-100000
success_AUC<-0
success_Cthrough<-0
success_one<-0
success_both<-0
for (i in 1:n_sim){
  t<-sim_trial()
  if(t[3,2]>0.8 & t[3,1]>0.8){
    success_AUC<-success_AUC+1
    success_one<-success_one+1
    success_Cthrough<-success_Cthrough+1
    success_both<-success_both+1
  }else{
    if(t[2,1]>0.8){
      success_Cthrough<-success_Cthrough+1
      success_one<-success_one+1
    }
    if(t[2,2]>0.8){
      success_AUC<-success_AUC+1
      success_one<-success_one+1
    }
  }
}

print(success_AUC/n_sim)
print(success_Cthrough/n_sim)
print(success_one/n_sim)
print(success_both/n_sim)
