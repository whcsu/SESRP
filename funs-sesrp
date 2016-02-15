srp2<-function (formula, data,trlength,blocksize=100,k=2,control = control, na.action =  na.omit,tra=FALSE) 
{

  tim <- format(Sys.time(), format = "%Y%j%H%M%S") 
  filen=paste0("spmat",tim)
  pd=pData(data)
  time=pd$time
  status=pd$status
  
 # pd[,5]=as.factor(pd[,5])
  #exp=t(exprs(data))
  #tm=data.frame(pd,exp)
  
  
  if (any(time <= 0)) stop("Observation time must be > 0")
  
  if (all(status == 0)) stop("No deaths in data set")
  
  
  if (!missing(control)) 
    controls[names(control)] <- control
  
  p1=dim(data)[1]
  p2=dim(pd)[2]
 #excluding time+status variable
  
  
  
  rotms1<- vector(mode = "list", length = blocksize)
 
  rotms2<- vector(mode = "list", length = blocksize)
 
  
  trees <- vector(mode = "list", length = trlength)
  
  oobci<-NULL
  varimp<-NULL 
  
  for (i in 1:trlength)
  {
    
    
    
    trainindex=sample(ncol(data),replace=T)
    
    #create a bagging version for rotation
    trainb=data[,trainindex]      
    
    # oobset  
    train_posp<-1:ncol(data) %in% trainindex
    
    oobset=data[,!train_posp]      
    
    #expdata
    bexp=exprs(trainb)
  
    #covariates
    bpd=pData(trainb)
    bpdy=bpd[,c(1,2)]
    bpdx=bpd[,-c(1,2)]
    
    
    bpdx[,2]=as.factor(bpdx[,2])
    bpdx[,3]=as.factor(bpdx[,3])
    bpdx=data.matrix( bpdx)
    mode(bpdx)="numeric"
    
    # rotation for gr i bagging version  
    rp=p1
     
    
    #random projections
    reduced=floor(sqrt(rp))
    
    # Projection matrix using The Achlioptas projection
    RPA <- floor(runif(rp*reduced,1,7)); # generate a vector 1 to 6 valued
    sqr3 <- sqrt(3);
    RPA[RPA==1] <- sqr3;
    RPA[RPA==6] <- -sqr3;
    RPA[RPA==2 | RPA==3 | RPA==4 | RPA==5] <- 0;
    rotationm1 <- matrix(RPA, ncol=reduced,nrow=rp);
    
    
    rp=p2-2
      
    
    #random projections
    reduced=floor(sqrt(rp))
    
    # Projection matrix using The Achlioptas projection
    RPA <- floor(runif(rp*reduced,1,7)); # generate a vector 1 to 6 valued
    sqr3 <- sqrt(3);
    RPA[RPA==1] <- sqr3;
    RPA[RPA==6] <- -sqr3;
    RPA[RPA==2 | RPA==3 | RPA==4 | RPA==5] <- 0;
    rotationm2 <- matrix(RPA, ncol=reduced,nrow=rp);
    
    
    
    
    # the final the rotation matrix 
    
    
    if (i%%blocksize==0)
    {
      rotms1[[blocksize]]=rotationm1
      rotms2[[blocksize]]=rotationm2
      save(rotms1,rotms2,file=paste(filen,floor(i/blocksize),sep=""))
    }
    else
      {
        rotms1[[i%%blocksize]]=rotationm1
        rotms2[[i%%blocksize]]=rotationm2
      }
   
    
    
    # the final bagging training set 
  
    
    rmatrix1=as.matrix(t(bexp)) %*% (rotationm1) 
    
    rmatrix2=bpdx%*%(rotationm2) 
    rmatrixnew1=as.data.frame(rmatrix1)
    rmatrixnew2=as.data.frame(rmatrix2)
    names(rmatrixnew2)=names(bpdx[,1:reduced])
    rmatrixnew=cbind(bpdy,rmatrixnew1,rmatrixnew2)  
    #print(names(rmatrixnew))
    #    names(rmatrixnew)=names(trainb)	       
    trees[[i]]=rpart(Surv(time,status)~.,data=rmatrixnew,control = control)
    
    if (tra)   print(i) 
      
    
    
  }
  if (i%%blocksize!=0)  
    save(rotms1,rotms2,,file=paste(filen,floor(i/blocksize)+1,sep=""))    
  
  rm(rotms1)
  rm(rotms2)
  
  fit=trees
  vari=NULL
  #make the difference apparent
  
  
  #print(oobci)
  ooberror=NULL
  
  class(fit) <- "rotsfhd"  
  
  
  return(list(trees=trees,blocksize=blocksize,filen=filen,ooberror=ooberror))
  
}



srp2.predict<-function(srpfit,newdata,trlength=100,del=TRUE){
  
  trees=srpfit$trees
  blocksize=srpfit$blocksize
  filen=srpfit$filen
  
  
  if (trlength>length(trees)) stop("Number of Trees for prediction should not be more than Number of Trees Fitted")
  
  # classify the test data   
  testpre<-NULL
  
  load(file=paste(filen,1,sep=""))   
  
  #test expdata
  texp=exprs(newdata)
  
  #covariates
  tpd=pData(newdata)  
  
  tpd[,2]=as.factor(tpd[,2])
  tpd[,3]=as.factor(tpd[,3])
  tpd=data.matrix( tpd)
  mode(tpd)="numeric"
  
  
  
  
  
  
  for (j in 1:trlength) {
    #if (oobacc[i]<=avroobacc)
    #print(j)
    
    # preparing for testing 
    if (j%%blocksize==0) {
      
      testmatrix1=as.matrix(t(texp)) %*% (rotms1[[(j-1)%%blocksize+1]])
      
      testmatrix2=tpd%*%(rotms2[[(j-1)%%blocksize+1]]) 
      testmatrix1=as.data.frame(testmatrix1)
      testmatrix2=as.data.frame(testmatrix2)
      names(testmatrix2)=names(tpd[,1:dim(testmatrix2)[2]])
      testmatrixnew=cbind(testmatrix1,testmatrix2) 
      
  
      
      if(j!=trlength) 
        load(file=paste(filen,floor(j/blocksize)+1,sep=""))
    }
    else
    {
      testmatrix1=as.matrix(t(texp)) %*% (rotms1[[j%%blocksize]])
      
      testmatrix2=tpd%*%(rotms2[[j%%blocksize+1]]) 
      testmatrix1=as.data.frame(testmatrix1)
      testmatrix2=as.data.frame(testmatrix2)
    
      names(testmatrix2)=names(tpd[,1:dim(testmatrix2)[2]])
      testmatrixnew=cbind(testmatrix1,testmatrix2)       
     
    }
    testmatrixnew=as.data.frame(testmatrixnew)
    
    # names(rmatrixnew2)=names(newdata)  
    
    predicts<-predict(trees[[j]],testmatrixnew)
    
    testpre<-cbind(predicts,testpre)
    
    
    
    
  }
  
  ensemble_predictions<-rowMeans(testpre)
  
  #delete temp files
  
  if (del)
    unlink(paste0(filen,"*"))
  
  return(ensemble_predictions)
  
}
