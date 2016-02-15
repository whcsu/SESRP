rm(list=ls())
#delete temp files

unlink("spmat*")

library(rpart)
library(MASS)
library(gbm)
library(Matrix)
library(survival)
#library(Hmisc)
library(survcomp)
library(survivalROC)
#library(ff)

#for High-dimensional 
library(affy)
#for dnet TCGA_mutations data
#library(foreach)
#library(doSNOW)

#cl <- makeCluster(2, type="SOCK")
#registerDoSNOW(cl)



set.seed(123)

datasetname="transbig"

library(breastCancerTRANSBIG)
data("transbig")
transpdata=phenoData(transbig)[,-c(1:5,10:12,16,17,20,21)]
#change the covariates order
transpdata=transpdata[ ,colnames(transpdata)[c(6,5,1:4,7:9)]]
colnames(transpdata)[1] <- "time"
colnames(transpdata)[2] <- "status"

# log 2 transform
tm=ExpressionSet(assayData=(exprs(transbig)),phenoData=transpdata)
#delete all missing values



tm=tm[,rowSums(is.na(pData(tm)))==0]


#delere necessary variables to save memorys
rm(transbig)


tm=tm[,pData(tm)$time>0]



#save(tm,file="mainz.RData")
#sample covariates 3096 samples * 7 covariates 
#pd=pData(TCGA_mutations)
#pd=pData(tm)
#gene expression data
#time=pd[,1]
#status=pd[,2]
#pd=pd[,-c(1,2)]

#save(tm,file="tm.RData")



#write.csv(tm,"tm.csv")


n=dim(tm)[2]


Rn=10
#testresult <-data.frame(row.names=1:Rn)

#AUC

#AUC
ci_srp<-c(rep(0,Rn))

#ci_rsf<-c(rep(0,Rn))
#ci_cox<-c(rep(0,Rn))
#ci_gbm<-c(rep(0,Rn))
source("funs-srp.R")

for (i in seq(1, Rn, by=2))
{
  print(i)
  
  L=sample(1:n,ceiling(n*0.5))
  trset<-tm[,L]
  teset<-tm[,-L]
  
  filen1=paste0(datasetname,"trset",i)
  filen2=paste0(datasetname,"teset",i)
  save(trset,file=filen1)
  save(teset,file=filen2)
  
  print("SESRP")
  #delere necessary variables to save memory 
  tesetx= ExpressionSet(assayData=exprs(teset),phenoData=phenoData(teset)[,-c(1,2)])
  
  print(Sys.time())
  srp.fit=srp2(Surv( time, status)~.,data=trset,trlength=500,blocksize=10)
  
  srp_pre=srp2.predict(srp.fit,tesetx,trlength=500)
  
  ci=concordance.index(srp_pre,pData(teset)$time,pData(teset)$status)
  print(ci[1])
  ci_srp[i]=unlist(ci[1])
  rm(srp.fit,srp_pre)
  gc()
  
  
  
  print(i+1)
  
  print("Fold 2: sesrp")
  
  tesetx= ExpressionSet(assayData=exprs(trset),phenoData=phenoData(trset)[,-c(1,2)])
  
  print(Sys.time())
  srp.fit=srp2(Surv( time, status)~.,data=teset,trlength=500,blocksize=10)
  
  srp_pre=srp2.predict(srp.fit,tesetx,trlength=500)
  
  ci=concordance.index(srp_pre,pData(trset)$time,pData(trset)$status)
  print(ci[1])
  ci_srp[i+1]=unlist(ci[1])
  rm(srp.fit,srp_pre)
  gc()
  print(Sys.time())
  
  
}



#testresult=cbind(testresult,data.frame(unlist(ci_rotsf),unlist(ci_rsf),unlist(ci_cox),unlist(ci_gbm)))

csvfn=paste0("srp-",datasetname,"-n5.csv")
write.csv(ci_srp,csvfn)

