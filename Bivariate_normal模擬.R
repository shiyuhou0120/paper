library(mvtnorm)
library(survival)
library(survAUC)
library(risksetROC)
library(survivalROC)
library(parallel)
library(pbapply)
library(timeROC)
library(pbapply)
library(parallel)

########################### n = 500,m = 1000
#Cox
##20%,n=500,m=1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 1.19;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores-2)
clusterExport(cl, c("mean","sigma","meanlog","sdlog","pt","MyrisksetAUC.wai.new","MyrisksetAUC.wai.old",
                    "EST.cindex.std","CoxWeights","CoxWeights.new","rmvnorm"))
clusterEvalQ(cl, c(library(aftgee),library(survival),library(risksetROC),library(timereg)))

# 無信賴區間估計 ####
clusterSetRNGStream(cl, 0120)     
PH.bn.20.500 = pblapply(X = lapply(1:1000,function(x) c(500)),FUN = function(n){
  AUC.bn.PH.new = Cindex.bn.PH.new = c()
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  beta.ph <- coxph(Surv(Stime,status)~M)$coefficients
  fit.ph.new = MyrisksetAUC.wai.new( Stime , entry = NULL, status, beta.ph*M, "Cox"
                                     , tmax = max(pt) , AFT.method.1 = "1",AFT.method.2 = "1"
                                     , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                     , weight = "rescale", plot=T, type="l",xlab="Time", ylab="AUC")
  for( r in 1:9){
    AUC.bn.PH.new[r] = fit.ph.new$AUC[which.min(abs(fit.ph.new$utimes-exp(-2.5+0.5*r)))]
  }
  Cindex.bn.PH.new = fit.ph.new$Cindex.new
  A = data.frame(cox = c(AUC.bn.PH.new,Cindex.bn.PH.new))
  return(A)
},cl = cl)

c = matrix(unlist(PH.bn.20.500,use.names = F),nrow = 1000,byrow = T)
round(apply(c,2,mean),3)
round(apply(c,2,sd),3)

#40%
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 0.36;sdlog = 1

cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores-2)
clusterExport(cl, c("mean","sigma","meanlog","sdlog","pt","MyrisksetAUC.wai.new","MyrisksetAUC.wai.old",
                    "EST.cindex.std","CoxWeights","CoxWeights.new","rmvnorm"))
clusterEvalQ(cl, c(library(aftgee),library(survival),library(risksetROC),library(timereg)))

clusterSetRNGStream(cl, 0120)     
PH.bn.40.500 = pblapply(X = lapply(1:1000,function(x) c(500)),FUN = function(n){
  AUC.bn.PH.new = Cindex.bn.PH.new = c()
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  beta.ph <- coxph(Surv(Stime,status)~M)$coefficients
  fit.ph.new = MyrisksetAUC.wai.new( Stime , entry = NULL, status, beta.ph*M, "Cox"
                                     , tmax = max(pt) , AFT.method.1 = "1",AFT.method.2 = "1"
                                     , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                     , weight = "rescale", plot=T, type="l",xlab="Time", ylab="AUC")
  for( r in 1:9){
    AUC.bn.PH.new[r] = fit.ph.new$AUC[which.min(abs(fit.ph.new$utimes-exp(-2.5+0.5*r)))]
  }
  Cindex.bn.PH.new = fit.ph.new$Cindex.new
  A = data.frame(cox = c(AUC.bn.PH.new,Cindex.bn.PH.new))
  return(A)
},cl = cl)

Cox.40 = matrix(unlist(PH.bn.40.500,use.names = F),nrow = 1000,byrow = T)
round(apply(Cox.40,2,mean),3)
round(apply(Cox.40,2,sd),3)
old

#residual smooth
##20%,n=500,m=1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 1.19;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores-1)
clusterExport(cl, c("mean","sigma","meanlog","sdlog","pt","rmvnorm","weightedKM","IntegrateAUC","risksetROC.new",
                    "CoxWeights.new"))
clusterEvalQ(cl, c(library(risksetROC),library(survival)))
clusterSetRNGStream(cl, 0120)
res = pblapply(X = lapply(1:1000,function(x) c(500)),FUN = function(n){
  res20.500 = AUC.bn.res.new = c()
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  for (i in 1:length(utimes)) {
    res.fit = risksetROC.new(Stime = Stime,status = status,marker = M,predict.time = utimes[i],
                             method = "Schoenfeld",type = "b",span = n^(-0.2))
    res20.500[i] = res.fit$AUC
  }
  for( r in 1:9){
    AUC.bn.res.new[r] = res20.500[which.min(abs(utimes-exp(-2.5+0.5*r)))]
  }
  km.out = weightedKM(Stime = Stime , status = status)
  Cindex = IntegrateAUC(res20.500 , utimes , km.out$survival , max(utimes), weight = "rescale") #concordance估計
  A = data.frame(res.20 = c(AUC.bn.res.new,Cindex))
  return(A)
},cl = cl)

r = matrix(unlist(res,use.names = F),nrow = 1000,byrow = T)
round(apply(r,2,mean),3)
res.old
round(apply(r,2,sd),3)


##40%,n=500,m=1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 0.36;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores <- detectCores() 
cl <- makeCluster(cpu.cores-1)
clusterExport(cl, c("mean","sigma","meanlog","sdlog","pt","rmvnorm","weightedKM",
                    "IntegrateAUC","risksetROC.new","CoxWeights.new"))
clusterEvalQ(cl, c(library(risksetROC),library(survival)))
clusterSetRNGStream(cl, 0120)
res.40 = pblapply(X = lapply(1:1000,function(x) c(500)),FUN = function(n){
  res40.500 = AUC.bn.res.new = c()
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  for (i in 1:length(utimes)) {
    res.fit = risksetROC.new(Stime = Stime,status = status,marker = M,predict.time = utimes[i],
                             method = "Schoenfeld",type = "b",span = n^(-0.2))
    res40.500[i] = res.fit$AUC
  }
  for( r in 1:9){
    AUC.bn.res.new[r] = res40.500[which.min(abs(utimes-exp(-2.5+0.5*r)))]
  }
  km.out = weightedKM(Stime = Stime , status = status)
  Cindex = IntegrateAUC(res40.500 , utimes , km.out$survival , max(utimes), weight = "rescale") #concordance估計
  A = data.frame(res.40 = c(AUC.bn.res.new,Cindex))
  return(A)
},cl = cl)

r.40= matrix(unlist(res.40,use.names = F),nrow = 1000,byrow = T)
round(apply(r.40,2,mean),3)
round(apply(r.40,2,sd),3)


#PO
##20%,n=500,m=1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 1.19;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-1)
clusterExport(cl, c("mean","weightedKM","EST.cindex.std","sigma","meanlog","sdlog","MyrisksetAUC.wai.new",
                    "CoxWeights","MyrisksetAUC.wai.old","IntegrateAUC","POWeights","POWeights.new",
                    "AFTWeights.1","prop.odds","Event","rmvnorm","pt"))
clusterEvalQ(cl, c(library(muhaz),library(aftgee),library(survival),library(risksetROC),library(timereg)))

clusterSetRNGStream(cl, 0120)     
EST.PO.bn.20.500 = pblapply(X = lapply(1:1000,function(x) c(500)),FUN = function(n){
  AUC.bn.PO.new = c()
  s = rmvnorm(n,mean,sigma)
  t = exp(s[,2])
  M = s[,1]
  c = rlnorm(n, meanlog, sdlog)
  status <- ifelse(t<c,1,0)
  Stime = ifelse(t<c,t,c)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  data = data.frame(time=rep(0,length(Stime)),time2=Stime,status = status,eta = M)
  fit = prop.odds(Event(time2,status)~M,max.time=max(Stime),data = data,n.sim=100)
  gamma.t = as.numeric(fit$gamma)
  fit.po.new = MyrisksetAUC.wai.new(Stime , entry = NULL, status, gamma.t*M, "PO" 
                                    , max(Stime) , AFT.method.1 = "1",AFT.method.2 = "1"
                                    , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                    , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  for( r in 1:9){
    AUC.bn.PO.new[r] = fit.po.new$AUC[which.min(abs(fit.po.new$utimes-exp(-2.5+0.5*r)))]
  }
  Cindex.bn.PO.np = fit.po.new$Cindex.new
  A = data.frame(EST = c(AUC.bn.PO.new,Cindex.bn.PO.np))
  return(A)
},cl = cl)

PO20.500 = matrix(unlist(EST.PO.bn.20.500,use.names = F),nrow = 1000,byrow = T)

round(apply(PO20.500,2,mean),3)
bn.PO.20.500
round(apply(PO20.500,2,sd),3)

##40%,n=500,m=1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 0.36;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-2)
clusterExport(cl, c("mean","weightedKM","EST.cindex.std","sigma","meanlog","sdlog","MyrisksetAUC.wai.new",
                    "CoxWeights","MyrisksetAUC.wai.old","IntegrateAUC","POWeights","POWeights.new",
                    "AFTWeights.1","prop.odds","Event","rmvnorm","pt"))
clusterEvalQ(cl, c(library(muhaz),library(aftgee),library(survival),library(risksetROC),library(timereg)))

clusterSetRNGStream(cl, 0120)     
EST.PO.bn.40.500 = pblapply(X = lapply(1:1000,function(x) c(500)),FUN = function(n){
  AUC.bn.PO.new = c()
  s = rmvnorm(n,mean,sigma)
  t = exp(s[,2])
  M = s[,1]
  c = rlnorm(n, meanlog, sdlog)
  pt = exp(seq(-2,2,0.5))
  status <- ifelse(t<c,1,0)
  Stime = ifelse(t<c,t,c)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  data = data.frame(time=rep(0,length(Stime)),time2=Stime,status = status,eta = M)
  fit = prop.odds(Event(time2,status)~M,max.time=max(Stime),data = data,n.sim=100)
  gamma.t = as.numeric(fit$gamma)
  fit.po.new = MyrisksetAUC.wai.new(Stime , entry = NULL, status, gamma.t*M, "PO" 
                                    , max(Stime) , AFT.method.1 = "1",AFT.method.2 = "1"
                                    , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                    , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  for( r in 1:9){
    AUC.bn.PO.new[r] = fit.po.new$AUC[which.min(abs(fit.po.new$utimes-exp(-2.5+0.5*r)))]
  }
  Cindex.bn.PO.np = fit.po.new$Cindex.new
  A = data.frame(EST = c(AUC.bn.PO.new,Cindex.bn.PO.np))
  return(A)
},cl = cl)

PO40.500 = matrix(unlist(EST.PO.bn.40.500,use.names = F),nrow = 1000,byrow = T)

round(apply(PO40.500,2,mean),3)
bn.PO.40.500
round(apply(PO40.500,2,sd),3)

#AFT
##20%,n=500,m=1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 1.19;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-2)
clusterExport(cl, c("MyrisksetAUC.wai.new","AFTWeights.1","AFTWeights.2","EST.cindex.std",
                    "mean","sigma","meanlog","sdlog","rmvnorm","pt","AFTWeights.1.new","AFTWeights.2.new"))
clusterEvalQ(cl, c(library(muhaz),library(aftgee),library(survival),library(risksetROC)))

clusterSetRNGStream(cl, 0120)     
ANS.AFT.bn.20.500 = pblapply(X = lapply(1:1000,function(x) c(500)),FUN = function(n){
  AUC.bn.AFT.new = c()
  s = rmvnorm(n,mean,sigma)
  t = exp(s[,2])
  M = s[,1]
  c = rlnorm(n, meanlog, sdlog)
  Stime = ifelse(t<c,t,c)
  status = ifelse(t<c,1,0)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  beta.aft <- -aftsrr(Surv(Stime,status)~M)$beta
  fit.aft.new = MyrisksetAUC.wai.new(Stime , entry = NULL, status, beta.aft*M, "AFT" 
                                     , max(Stime) , AFT.method.1 = "1", AFT.method.2 = "1"
                                     , bw.1="1/2",bw.2 = "1/2", bw.mul.1=10, bw.mul.2=0.5, CI=F, confi=0.95
                                     , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  for( r in 1:9){
    AUC.bn.AFT.new[r] = fit.aft.new$AUC[which.min(abs(fit.aft.new$utimes-exp(-2.5+0.5*r)))]
  }
  Cindex.bn.AFT = fit.aft.new$Cindex.new
  A = data.frame(EST = c(AUC.bn.AFT.new,Cindex.bn.AFT))
  return(A)
},cl = cl)

AFT20.500 = matrix(unlist(ANS.AFT.bn.20.500,use.names = F),nrow = 1000,byrow = T)
which(AFT20.500<0.5)
which(is.na(AFT20.500))

round(apply(AFT20.500,2,mean),3)
bn.AFT.20.500
round(apply(AFT20.500,2,sd),3)

##40%,n=500,m=1000 
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 0.36;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-1)
clusterExport(cl, c("MyrisksetAUC.wai.new","AFTWeights.1","AFTWeights.2","EST.cindex.std",
                    "mean","sigma","meanlog","sdlog","rmvnorm","pt","AFTWeights.1.new","AFTWeights.2.new"))
clusterEvalQ(cl, c(library(muhaz),library(aftgee),library(survival),library(risksetROC)))

clusterSetRNGStream(cl, 0120)     
ANS.AFT.bn.40.500 = pblapply(X = lapply(1:1000,function(x) c(500)),FUN = function(n){
  AUC.bn.AFT.new = c()
  s = rmvnorm(n,mean,sigma)
  t = exp(s[,2])
  M = s[,1]
  c = rlnorm(n, meanlog, sdlog)
  pt = exp(seq(-2,2,0.5))
  Stime = ifelse(t<c,t,c)
  status = ifelse(t<c,1,0)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  beta.aft <- -aftsrr(Surv(Stime,status)~M)$beta
  fit.aft.new = MyrisksetAUC.wai.new(Stime , entry = NULL, status, beta.aft*M, "AFT" 
                                     , max(Stime) , AFT.method.1 = "1", AFT.method.2 = "1"
                                     , bw.1="1/2",bw.2 = "1/2", bw.mul.1=10, bw.mul.2=0.5, CI=F, confi=0.95
                                     , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  for( r in 1:9){
    AUC.bn.AFT.new[r] = fit.aft.new$AUC[which.min(abs(fit.aft.new$utimes-exp(-2.5+0.5*r)))]
  }
  Cindex.bn.AFT = fit.aft.new$Cindex.new
  A = data.frame(EST = c(AUC.bn.AFT.new,Cindex.bn.AFT))
  return(A)
},cl = cl)

AFT40.500 = matrix(unlist(ANS.AFT.bn.40.500,use.names = F),nrow = 1000,byrow = T)
which(AFT40.500<0.5)
which(is.na(AFT40.500))

round(apply(AFT40.500,2,mean),3)
round(apply(AFT40.500,2,sd),3)
bn.AFT.40.500
bn.AFT.40.500.sd

##method = KM,CD1 
##20%
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 1.19;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-2)
clusterExport(cl, c("mean","sigma","meanlog","sdlog","rmvnorm","pt"))
clusterEvalQ(cl, c(library(survival),library(risksetROC),library(survivalROC)))
clusterSetRNGStream(cl, 0120)
CD1.20 = pblapply(X = lapply(1:1000,function(x) c(500)),FUN = function(n){
  CD1.20.500 = c()
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  for (i in 1:9) {
    CD1.20.500[i] = survivalROC(Stime = Stime,status = status,marker = M, method = "KM",
                                predict.time = utimes[which.min(abs(utimes-exp(-2.5+0.5*i)))])$AUC
  }
  A = data.frame(CD1 = CD1.20.500)
  return(A)
},cl = cl)

CD1.20.500 = matrix(unlist(CD1.20,use.names = F),nrow = 1000,byrow = T)
which(CD1.20.500<0.5)
which(is.na(CD1.20.500))

round(apply(na.omit(CD1.20.500),2,mean),3)
bn.CD1.20.500
round(apply(na.omit(CD1.20.500),2,sd),3)

##40%
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 0.36;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-2)
clusterExport(cl, c("mean","sigma","meanlog","sdlog","rmvnorm"))
clusterEvalQ(cl, c(library(survival),library(risksetROC),library(survivalROC)))
clusterSetRNGStream(cl, 0120)
CD1.40 = pblapply(X = lapply(1:1000,function(x) c(500)),FUN = function(n){
  CD1.40.500 = c()
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  for (i in 1:9) {
    CD1.40.500[i] = survivalROC(Stime = Stime,status = status,marker = M, method = "KM",
                                predict.time = utimes[which.min(abs(utimes-exp(-2.5+0.5*i)))])$AUC
  }
  A = data.frame(CD1 = CD1.40.500)
  return(A)
},cl = cl)

CD1.40.500 = matrix(unlist(CD1.40,use.names = F),nrow = 1000,byrow = T)
which(CD1.40.500<0.5)
which(is.na(CD1.40.500))

round(apply(na.omit(CD1.40.500),2,mean),3)
bn.CD1.40.500
round(apply(na.omit(CD1.40.500),2,sd),3)

#method = NNE,CD2
#20%
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 1.19;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-2)
clusterExport(cl, c("mean","sigma","meanlog","sdlog","rmvnorm"))
clusterEvalQ(cl, c(library(survival),library(risksetROC),library(survivalROC)))
clusterSetRNGStream(cl, 0120)
CD2.20 = pblapply(X = lapply(1:1000,function(x) c(500)),FUN = function(n){
  CD2.20.500 = c()
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  for (i in 1:9) {
    CD2.20.500[i] = survivalROC(Stime = Stime,status = status,marker = M, 
                                predict.time = utimes[which.min(abs(utimes-exp(-2.5+0.5*i)))],
                                method = "NNE",span = n^(-0.2))$AUC
  }
  A = data.frame(CD2 = CD2.20.500)
  return(A)
},cl = cl)

CD2.20.500 = matrix(unlist(CD2.20,use.names = F),nrow = 1000,byrow = T)
which(CD2.20.500<0.5)
which(is.na(CD2.20.500))

round(apply(na.omit(CD2.20.500),2,mean),3)
bn.CD2.20.500
round(apply(na.omit(CD2.20.500),2,sd),3)

##40%
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 0.36;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-2)
clusterExport(cl, c("mean","sigma","meanlog","sdlog","rmvnorm"))
clusterEvalQ(cl, c(library(survival),library(risksetROC),library(survivalROC)))
clusterSetRNGStream(cl, 0120)
CD2.40 = pblapply(X = lapply(1:1000,function(x) c(500)),FUN = function(n){
  CD2.40.500 = c()
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  for (i in 1:9) {
    CD2.40.500[i] = survivalROC(Stime = Stime,status = status,marker = M, 
                                predict.time = utimes[which.min(abs(utimes-exp(-2.5+0.5*i)))],
                                method = "NNE",span = n^(-0.2))$AUC
  }
  A = data.frame(CD2 = CD2.40.500)
  return(A)
},cl = cl)

CD2.40.500 = matrix(unlist(CD2.40,use.names = F),nrow = 1000,byrow = T)
which(CD2.40.500<0.5)
which(is.na(CD2.40.500))

round(apply(na.omit(CD2.40.500),2,mean),3)
round(apply(na.omit(CD2.40.500),2,sd),3)

##檢查為何有NA
s = rmvnorm(n,mean,sigma)
T = exp(s[,2])
M = s[,1]
C = rlnorm(n, meanlog, sdlog)
status = ifelse(T<C,1,0)
Stime = ifelse(T<C,T,C)
abc = survivalROC(Stime = Stime,
                     status = status,
                     marker = M, 
                     predict.time = pt[9],
                     method = "NNE",
                     span = n^(-0.2)) ##Span for the NNE(or lambda)
abc$AUC

--------------------------------------------------------##沒用到
##CD5,with a Cox model for computing the weights.
##20%
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 1.19;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-2)
clusterExport(cl, c("mean","sigma","meanlog","sdlog","rmvnorm","pt"))
clusterEvalQ(cl, c(library(survival),library(risksetROC),library(survivalROC),library(timeROC)))
clusterSetRNGStream(cl, 0120)
CD5.20 = pblapply(X = lapply(1:1000,function(x) c(500)),FUN = function(n){
  CD5.20.500 = c()
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  
  CD5.20.500 = timeROC(T = Stime,delta = status,marker = M,cause = 1,
                       weighting = "cox",times = pt)$AUC
  A = data.frame(CD5 = CD5.20.500)
  return(A)
},cl = cl)

CD5.20.500 = matrix(unlist(CD5.20,use.names = F),nrow = 1000,byrow = T)
which(CD5.20.500<0.5)
which(is.na(CD5.20.500))

round(apply(na.omit(CD5.20.500),2,mean),3)
round(apply(na.omit(CD5.20.500),2,sd),3)

m=1
CD5.20.500 = matrix(rep(0,length(pt)*m),nrow = m,byrow = T)
#CD6.20 = matrix(rep(0,length(pt)*m),nrow = m,byrow = T)
cen = c()
for (j in 1:m) {
  p.time = c()
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  for (i in 1:9) {
    p.time[i] = utimes[which.min(abs(utimes-exp(-2.5+0.5*i)))]
  }
  fit_C5 = timeROC(T = Stime,delta = status,marker = M,cause = 1,
                   weighting = "cox",times = p.time)
  #fit_C6 = timeROC(T = Stime,delta = status,marker = M,cause = 1,
  #                 weighting = "marginal",times = pt,ROC = F,iid = TRUE)
  CD5.20.500[j,] = fit_C5$AUC
  #CD6.20[j,] = fit_C6$AUC
}
mean(cen)
CD5.20.500
which(CD5.20.500<0.5)
which(is.na(CD5.20.500))
#CD6.20

##將AUC=NA的改成0.5
CD5.20.500[which(is.na(CD5.20.500))] = 0.5
##將AUC<0.5的改成0.5
CD5.20.500[which(CD5.20.500<0.5)] = 0.5

bn.CD5.20.nci.500 = data.frame(EST = round(apply(CD5.20.500,2,mean),3))
#bn.CD6.20.nci.200 = data.frame(EST = round(apply(C6,2,mean),3))

bn.CD5.SD.20.500 = sapply(1:9,function(i)sd(sapply(1:length(CD5.20.500[,1]),function(k) CD5.20.500[k,i])))
bn.CD5.SD.20.500 = data.frame(SD = round(bn.CD5.SD.20.500,3))
bn.CD5.SD.20.500
cbind(bn.CD5.20.nci.500,bn.CD5.SD.20.500)

##40%
set.seed(2021)
n = 500;m = 1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 0.36;sdlog = 1
pt = exp(seq(-2,2,0.5))
CD5.40.500 = matrix(rep(0,length(pt)*m),nrow = m,byrow = T)
##CD6.40 = matrix(rep(0,length(pt)*m),nrow = m,byrow = T)
cen = c()
for (j in 1:m) {
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  cen[j] = table(status)[1]/sum(table(status))
  fit_C5 = timeROC(T = Stime,delta = status,marker = M,cause = 1,
                   weighting = "cox",times = pt,ROC = F)
  ##fit_C6 = timeROC(T = Stime,delta = status,marker = M,cause = 1,
  ##                 weighting = "marginal",times = pt,ROC = F,iid = TRUE)
  CD5.40.500[j,] = fit_C5$AUC
  ##CD6.40[j,] = fit_C6$AUC
}
mean(cen)
CD5.40.500
which(CD5.40.500<0.5)
which(is.na(CD5.40.500))
##CD6.40

##將AUC=NA的改成0.5
CD5.40.500[which(is.na(CD5.40.500))] = 0.5
##將AUC<0.5的改成0.5
CD5.40.500[which(CD5.40.500<0.5)] = 0.5

bn.CD5.40.nci.500 = data.frame(EST = round(apply(CD5.40.500,2,mean),3))
#bn.CD6.40.nci.200 = data.frame(EST = round(apply(CD6.40,2,mean),3))

bn.CD5.SD.40.500 = sapply(1:9,function(i)sd(sapply(1:length(CD5.40.500[,1]),function(k) CD5.40.500[k,i])))
bn.CD5.SD.40.500 = data.frame(SD = round(bn.CD5.SD.40.500,3))
bn.CD5.SD.40.500

cbind(bn.CD5.40.nci.500,bn.CD5.SD.40.500)

n.500 = data.frame(Logtime = seq(-2,2,0.5),
                   Cox = bn.cox.20.nci.500,
                   Residual = bn.residual.20.nci.500,
                   PO = bn.PO.20.nci.500[-10,],
                   AFT = bn.AFT.20.nci.500[-10,],
                   CD1 = bn.CD1.20.nci.500,
                   CD2 = bn.CD2.20.nci.500,
                   CD5 = bn.CD5.20.nci.500)
write.csv(n.500,"C:\\Users\\user\\Desktop\\模擬paper\\模擬.csv")
-------------------------------------------------------------------------------------

########################### n = 200,m = 1000
#Cox
##20%,n=200,m=1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 1.19;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores-2)
clusterExport(cl, c("mean","sigma","meanlog","sdlog","pt","MyrisksetAUC.wai.new","MyrisksetAUC.wai.old",
                    "EST.cindex.std","CoxWeights","CoxWeights.new","rmvnorm"))
clusterEvalQ(cl, c(library(aftgee),library(survival),library(risksetROC),library(timereg)))

# 無信賴區間估計 ####
clusterSetRNGStream(cl, 0120)     
PH.bn.20.200 = pblapply(X = lapply(1:1000,function(x) c(200)),FUN = function(n){
  AUC.bn.PH.new = Cindex.bn.PH.new = c()
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  beta.ph <- coxph(Surv(Stime,status)~M)$coefficients
  fit.ph.new = MyrisksetAUC.wai.new( Stime , entry = NULL, status, beta.ph*M, "Cox"
                                     , tmax = max(pt) , AFT.method.1 = "1",AFT.method.2 = "1"
                                     , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                     , weight = "rescale", plot=T, type="l",xlab="Time", ylab="AUC")
  for( r in 1:9){
    AUC.bn.PH.new[r] = fit.ph.new$AUC[which.min(abs(fit.ph.new$utimes-exp(-2.5+0.5*r)))]
  }
  Cindex.bn.PH.new = fit.ph.new$Cindex.new
  A = data.frame(cox = c(AUC.bn.PH.new,Cindex.bn.PH.new))
  return(A)
},cl = cl)

c.200 = matrix(unlist(PH.bn.20.200,use.names = F),nrow = 1000,byrow = T)
round(apply(c.200,2,mean),3)
round(apply(c.200,2,sd),3)

#40%
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 0.36;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores-2)
clusterExport(cl, c("mean","sigma","meanlog","sdlog","pt","MyrisksetAUC.wai.new","MyrisksetAUC.wai.old",
                    "EST.cindex.std","CoxWeights","CoxWeights.new","rmvnorm"))
clusterEvalQ(cl, c(library(aftgee),library(survival),library(risksetROC),library(timereg)))

# 無信賴區間估計 ####
clusterSetRNGStream(cl, 0120)     
PH.bn.40.200 = pblapply(X = lapply(1:1000,function(x) c(200)),FUN = function(n){
  AUC.bn.PH.new = Cindex.bn.PH.new = c()
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  beta.ph <- coxph(Surv(Stime,status)~M)$coefficients
  fit.ph.new = MyrisksetAUC.wai.new( Stime , entry = NULL, status, beta.ph*M, "Cox"
                                     , tmax = max(pt) , AFT.method.1 = "1",AFT.method.2 = "1"
                                     , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                     , weight = "rescale", plot=T, type="l",xlab="Time", ylab="AUC")
  for( r in 1:9){
    AUC.bn.PH.new[r] = fit.ph.new$AUC[which.min(abs(fit.ph.new$utimes-exp(-2.5+0.5*r)))]
  }
  Cindex.bn.PH.new = fit.ph.new$Cindex.new
  A = data.frame(cox = c(AUC.bn.PH.new,Cindex.bn.PH.new))
  return(A)
},cl = cl)

c.200.40 = matrix(unlist(PH.bn.40.200,use.names = F),nrow = 1000,byrow = T)
round(apply(c.200.40,2,mean),3)
round(apply(c.200.40,2,sd),3)


#residual smooth
##20%,n=200,m=1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 1.19;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores-1)
clusterExport(cl, c("mean","sigma","meanlog","sdlog","pt","rmvnorm","weightedKM","IntegrateAUC","risksetROC.new",
                    "CoxWeights.new"))
clusterEvalQ(cl, c(library(risksetROC),library(survival)))
clusterSetRNGStream(cl, 0120)
res.20.200 = pblapply(X = lapply(1:1000,function(x) c(200)),FUN = function(n){
  res20.200 = AUC.bn.res.new = c()
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  for (i in 1:length(utimes)) {
    res.fit = risksetROC.new(Stime = Stime,status = status,marker = M,predict.time = utimes[i],
                             method = "Schoenfeld",type = "b",span = n^(-0.2))
    res20.200[i] = res.fit$AUC
  }
  for( r in 1:9){
    AUC.bn.res.new[r] = res20.200[which.min(abs(utimes-exp(-2.5+0.5*r)))]
  }
  km.out = weightedKM(Stime = Stime , status = status)
  Cindex = IntegrateAUC(res20.200 , utimes , km.out$survival , max(utimes), weight = "rescale") #concordance估計
  A = data.frame(res.200 = c(AUC.bn.res.new,Cindex))
  return(A)
},cl = cl)

r.20.200 = matrix(unlist(res.20.200,use.names = F),nrow = 1000,byrow = T)
round(apply(r.20.200,2,mean),3)
round(apply(r.20.200,2,sd),3)


##40%,n=200,m=1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 0.36;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores <- detectCores() 
cl <- makeCluster(cpu.cores-1)
clusterExport(cl, c("mean","sigma","meanlog","sdlog","pt","rmvnorm","weightedKM",
                    "IntegrateAUC","risksetROC.new","CoxWeights.new"))
clusterEvalQ(cl, c(library(risksetROC),library(survival)))
clusterSetRNGStream(cl, 0120)
res.40.200 = pblapply(X = lapply(1:1000,function(x) c(200)),FUN = function(n){
  res40.200 = AUC.bn.res.new = c()
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  for (i in 1:length(utimes)) {
    res.fit = risksetROC.new(Stime = Stime,status = status,marker = M,predict.time = utimes[i],
                             method = "Schoenfeld",type = "b",span = n^(-0.2))
    res40.200[i] = res.fit$AUC
  }
  for( r in 1:9){
    AUC.bn.res.new[r] = res40.200[which.min(abs(utimes-exp(-2.5+0.5*r)))]
  }
  km.out = weightedKM(Stime = Stime , status = status)
  Cindex = IntegrateAUC(res40.200 , utimes , km.out$survival , max(utimes), weight = "rescale") #concordance估計
  A = data.frame(res.200 = c(AUC.bn.res.new,Cindex))
  return(A)
},cl = cl)

r.40.200 = matrix(unlist(res.40.200,use.names = F),nrow = 1000,byrow = T)
round(apply(r.40.200,2,mean),3)
round(apply(r.40.200,2,sd),3)


#PO
##20%,n=200,m=1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 1.19;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-1)
clusterExport(cl, c("mean","weightedKM","EST.cindex.std","sigma","meanlog","sdlog","MyrisksetAUC.wai.new",
                    "CoxWeights","MyrisksetAUC.wai.old","IntegrateAUC","POWeights","POWeights.new",
                    "AFTWeights.1","prop.odds","Event","rmvnorm","pt"))
clusterEvalQ(cl, c(library(muhaz),library(aftgee),library(survival),library(risksetROC),library(timereg)))

clusterSetRNGStream(cl, 0120)     
EST.PO.bn.20.200 = pblapply(X = lapply(1:1000,function(x) c(200)),FUN = function(n){
  AUC.bn.PO.new = c()
  s = rmvnorm(n,mean,sigma)
  t = exp(s[,2])
  M = s[,1]
  c = rlnorm(n, meanlog, sdlog)
  status <- ifelse(t<c,1,0)
  Stime = ifelse(t<c,t,c)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  data = data.frame(time=rep(0,length(Stime)),time2=Stime,status = status,eta = M)
  fit = prop.odds(Event(time2,status)~M,max.time=max(Stime),data = data,n.sim=100)
  gamma.t = as.numeric(fit$gamma)
  fit.po.new = MyrisksetAUC.wai.new(Stime , entry = NULL, status, gamma.t*M, "PO" 
                                    , max(Stime) , AFT.method.1 = "1",AFT.method.2 = "1"
                                    , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                    , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  for( r in 1:9){
    AUC.bn.PO.new[r] = fit.po.new$AUC[which.min(abs(fit.po.new$utimes-exp(-2.5+0.5*r)))]
  }
  Cindex.bn.PO.np = fit.po.new$Cindex.new
  A = data.frame(EST = c(AUC.bn.PO.new,Cindex.bn.PO.np))
  return(A)
},cl = cl)

PO20.200 = matrix(unlist(EST.PO.bn.20.200,use.names = F),nrow = 1000,byrow = T)

round(apply(PO20.200,2,mean),3)
round(apply(PO20.200,2,sd),3)

##40%,n=200,m=1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 0.36;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-2)
clusterExport(cl, c("mean","weightedKM","EST.cindex.std","sigma","meanlog","sdlog","MyrisksetAUC.wai.new",
                    "CoxWeights","MyrisksetAUC.wai.old","IntegrateAUC","POWeights","POWeights.new",
                    "AFTWeights.1","prop.odds","Event","rmvnorm","pt"))
clusterEvalQ(cl, c(library(muhaz),library(aftgee),library(survival),library(risksetROC),library(timereg)))

clusterSetRNGStream(cl, 0120)     
EST.PO.bn.40.200 = pblapply(X = lapply(1:1000,function(x) c(200)),FUN = function(n){
  AUC.bn.PO.new = c()
  s = rmvnorm(n,mean,sigma)
  t = exp(s[,2])
  M = s[,1]
  c = rlnorm(n, meanlog, sdlog)
  pt = exp(seq(-2,2,0.5))
  status <- ifelse(t<c,1,0)
  Stime = ifelse(t<c,t,c)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  data = data.frame(time=rep(0,length(Stime)),time2=Stime,status = status,eta = M)
  fit = prop.odds(Event(time2,status)~M,max.time=max(Stime),data = data,n.sim=100)
  gamma.t = as.numeric(fit$gamma)
  fit.po.new = MyrisksetAUC.wai.new(Stime , entry = NULL, status, gamma.t*M, "PO" 
                                    , max(Stime) , AFT.method.1 = "1",AFT.method.2 = "1"
                                    , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                    , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  for( r in 1:9){
    AUC.bn.PO.new[r] = fit.po.new$AUC[which.min(abs(fit.po.new$utimes-exp(-2.5+0.5*r)))]
  }
  Cindex.bn.PO.np = fit.po.new$Cindex.new
  A = data.frame(EST = c(AUC.bn.PO.new,Cindex.bn.PO.np))
  return(A)
},cl = cl)

PO40.200 = matrix(unlist(EST.PO.bn.40.200,use.names = F),nrow = 1000,byrow = T)

round(apply(PO40.200,2,mean),3)
round(apply(PO40.200,2,sd),3)

#AFT
##20%,n=500,m=1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 1.19;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-2)
clusterExport(cl, c("MyrisksetAUC.wai.new","AFTWeights.1","AFTWeights.2","EST.cindex.std",
                    "mean","sigma","meanlog","sdlog","rmvnorm","pt","AFTWeights.1.new","AFTWeights.2.new"))
clusterEvalQ(cl, c(library(muhaz),library(aftgee),library(survival),library(risksetROC)))

clusterSetRNGStream(cl, 0120)     
ANS.AFT.bn.20.200 = pblapply(X = lapply(1:1000,function(x) c(200)),FUN = function(n){
  AUC.bn.AFT.new = c()
  s = rmvnorm(n,mean,sigma)
  t = exp(s[,2])
  M = s[,1]
  c = rlnorm(n, meanlog, sdlog)
  Stime = ifelse(t<c,t,c)
  status = ifelse(t<c,1,0)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  beta.aft <- -aftsrr(Surv(Stime,status)~M)$beta
  fit.aft.new = MyrisksetAUC.wai.new(Stime , entry = NULL, status, beta.aft*M, "AFT" 
                                     , max(Stime) , AFT.method.1 = "1", AFT.method.2 = "1"
                                     , bw.1="1/2",bw.2 = "1/2", bw.mul.1=10, bw.mul.2=0.5, CI=F, confi=0.95
                                     , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  for( r in 1:9){
    AUC.bn.AFT.new[r] = fit.aft.new$AUC[which.min(abs(fit.aft.new$utimes-exp(-2.5+0.5*r)))]
  }
  Cindex.bn.AFT = fit.aft.new$Cindex.new
  A = data.frame(EST = c(AUC.bn.AFT.new,Cindex.bn.AFT))
  return(A)
},cl = cl)

AFT20.200 = matrix(unlist(ANS.AFT.bn.20.200,use.names = F),nrow = 1000,byrow = T)
which(AFT20.200<0.5)
which(is.na(AFT20.200))

round(apply(AFT20.200,2,mean),3)
bn.AFT.20.200
round(apply(AFT20.200,2,sd),3)

##40%,n=500,m=1000 
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 0.36;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-1)
clusterExport(cl, c("MyrisksetAUC.wai.new","AFTWeights.1","AFTWeights.2","EST.cindex.std",
                    "mean","sigma","meanlog","sdlog","rmvnorm","pt","AFTWeights.1.new","AFTWeights.2.new"))
clusterEvalQ(cl, c(library(muhaz),library(aftgee),library(survival),library(risksetROC)))

clusterSetRNGStream(cl, 0120)     
ANS.AFT.bn.40.200 = pblapply(X = lapply(1:1000,function(x) c(200)),FUN = function(n){
  AUC.bn.AFT.new = c()
  s = rmvnorm(n,mean,sigma)
  t = exp(s[,2])
  M = s[,1]
  c = rlnorm(n, meanlog, sdlog)
  pt = exp(seq(-2,2,0.5))
  Stime = ifelse(t<c,t,c)
  status = ifelse(t<c,1,0)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  beta.aft <- -aftsrr(Surv(Stime,status)~M)$beta
  fit.aft.new = MyrisksetAUC.wai.new(Stime , entry = NULL, status, beta.aft*M, "AFT" 
                                     , max(Stime) , AFT.method.1 = "1", AFT.method.2 = "1"
                                     , bw.1="1/2",bw.2 = "1/2", bw.mul.1=10, bw.mul.2=0.5, CI=F, confi=0.95
                                     , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  for( r in 1:9){
    AUC.bn.AFT.new[r] = fit.aft.new$AUC[which.min(abs(fit.aft.new$utimes-exp(-2.5+0.5*r)))]
  }
  Cindex.bn.AFT = fit.aft.new$Cindex.new
  A = data.frame(EST = c(AUC.bn.AFT.new,Cindex.bn.AFT))
  return(A)
},cl = cl)

AFT40.200 = matrix(unlist(ANS.AFT.bn.40.200,use.names = F),nrow = 1000,byrow = T)
which(is.na(AFT40.200 < 0.5))
which(is.na(AFT40.200))

round(apply(AFT40.200,2,mean),3)
round(apply(AFT40.200,2,sd),3)
bn.AFT.40.200
bn.AFT.40.200.sd


##method = KM,CD1 
##20%
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 1.19;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-2)
clusterExport(cl, c("mean","sigma","meanlog","sdlog","rmvnorm","pt"))
clusterEvalQ(cl, c(library(survival),library(risksetROC),library(survivalROC)))
clusterSetRNGStream(cl, 0120)
CD1.200.20 = pblapply(X = lapply(1:1000,function(x) c(200)),FUN = function(n){
  CD1.20.200 = c()
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  for (i in 1:9) {
    CD1.20.200[i] = survivalROC(Stime = Stime,status = status,marker = M, method = "KM",
                                predict.time = utimes[which.min(abs(utimes-exp(-2.5+0.5*i)))])$AUC
  }
  A = data.frame(CD1 = CD1.20.200)
  return(A)
},cl = cl)

CD1.20.200 = matrix(unlist(CD1.200.20,use.names = F),nrow = 1000,byrow = T)
which(CD1.20.200<0.5)
which(is.na(CD1.20.200))

round(apply(na.omit(CD1.20.200),2,mean),3)
bn.CD1.20.500
round(apply(na.omit(CD1.20.200),2,sd),3)

##40%
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 0.36;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-2)
clusterExport(cl, c("mean","sigma","meanlog","sdlog","rmvnorm"))
clusterEvalQ(cl, c(library(survival),library(risksetROC),library(survivalROC)))
clusterSetRNGStream(cl, 0120)
CD1.200.40 = pblapply(X = lapply(1:1000,function(x) c(200)),FUN = function(n){
  CD1.40.200 = c()
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  for (i in 1:9) {
    CD1.40.200[i] = survivalROC(Stime = Stime,status = status,marker = M, method = "KM",
                                predict.time = utimes[which.min(abs(utimes-exp(-2.5+0.5*i)))])$AUC
  }
  A = data.frame(CD1 = CD1.40.200)
  return(A)
},cl = cl)

CD1.40.200 = matrix(unlist(CD1.200.40,use.names = F),nrow = 1000,byrow = T)
which(CD1.40.200<0.5)
which(is.na(CD1.40.200))

round(apply(na.omit(CD1.40.200),2,mean),3)
bn.CD1.40.500
round(apply(na.omit(CD1.40.200),2,sd),3)

#method = NNE,CD2
#20%
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 1.19;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-2)
clusterExport(cl, c("mean","sigma","meanlog","sdlog","rmvnorm","survivalROC.fix"))
clusterEvalQ(cl, c(library(survival),library(risksetROC),library(survivalROC)))
clusterSetRNGStream(cl, 0120)
CD2.200.20 = pblapply(X = lapply(1:1000,function(x) c(200)),FUN = function(n){
  CD2.20.200 = c()
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  for (i in 1:9) {
    CD2.20.200[i] = survivalROC(Stime = Stime,status = status,marker = M, 
                                predict.time = utimes[which.min(abs(utimes-exp(-2.5+0.5*i)))],
                                method = "NNE",span = n^(-0.2))$AUC
  }
  A = data.frame(CD2 = CD2.20.200)
  return(A)
},cl = cl)

CD2.20.200 = matrix(unlist(CD2.200.20,use.names = F),nrow = 1000,byrow = T)
which(CD2.20.200<0.5)
which(is.na(CD2.20.200))

round(apply(na.omit(CD2.20.200),2,mean),3)
bn.CD2.20.200
round(apply(na.omit(CD2.20.200),2,sd),3)

##40%
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 0.36;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-2)
clusterExport(cl, c("mean","sigma","meanlog","sdlog","rmvnorm"))
clusterEvalQ(cl, c(library(survival),library(risksetROC),library(survivalROC)))
clusterSetRNGStream(cl, 0120)
CD2.200.40 = pblapply(X = lapply(1:1000,function(x) c(200)),FUN = function(n){
  CD2.40.200 = c()
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  for (i in 1:9) {
    CD2.40.200[i] = survivalROC(Stime = Stime,status = status,marker = M, 
                                predict.time = utimes[which.min(abs(utimes-exp(-2.5+0.5*i)))],
                                method = "NNE",span = n^(-0.2))$AUC
  }
  A = data.frame(CD2 = CD2.40.200)
  return(A)
},cl = cl)

CD2.40.200 = matrix(unlist(CD2.200.40,use.names = F),nrow = 1000,byrow = T)
which(CD2.40.200<0.5)
which(is.na(CD2.40.200))

round(apply(na.omit(CD2.40.200),2,mean),3)
round(apply(na.omit(CD2.40.200),2,sd),3)

##檢查為何有NA
s = rmvnorm(n,mean,sigma)
T = exp(s[,2])
M = s[,1]
C = rlnorm(n, meanlog, sdlog)
status = ifelse(T<C,1,0)
Stime = ifelse(T<C,T,C)
abc = survivalROC(Stime = Stime,
                  status = status,
                  marker = M, 
                  predict.time = pt[9],
                  method = "NNE",
                  span = n^(-0.2)) ##Span for the NNE(or lambda)
abc$AUC



-----------------------------------------------------------------------------------------------
########################### n = 200,m = 1000
  
#Cox & residual smooth
##20%,n=200,m=1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 1.19;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores-1)
clusterExport(cl, c("mean","sigma","meanlog","sdlog","pt","rmvnorm"))
clusterEvalQ(cl, c(library(risksetROC),library(survival)))
clusterSetRNGStream(cl, 0120)
Cox.res = pblapply(X = lapply(1:1000,function(x) c(200)),FUN = function(n){
  cox20.200 = res20.200 = c()
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  for (i in 1:length(pt)) {
    fit1 = risksetROC(Stime = Stime,status = status,marker = M,predict.time = pt[i],
                      method = "Cox",type = "b")
    fit2 = risksetROC(Stime = Stime,status = status,marker = M,predict.time = pt[i],
                      method = "Schoenfeld",type = "b",span = n^(-0.2))
    cox20.200[i] = fit1$AUC
    res20.200[i] = fit2$AUC
  }
  A = data.frame(Cox = cox20.200,
                 res = res20.200)
  return(A)
},cl = cl)

a = matrix(unlist(Cox.res,use.names = F),nrow = 1000,byrow = T)
cox20.200 = a[,1:9]
res20.200 = a[,10:18]
round(apply(cox20.200,2,mean),3)
round(apply(res20.200,2,mean),3)
round(apply(cox20.200,2,sd),3)
round(apply(res20.200,2,sd),3)

##40%,n=500,m=1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 0.36;sdlog = 1
pt = exp(seq(-2,2,0.5))

cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores-1)
clusterExport(cl, c("mean","sigma","meanlog","sdlog","pt","rmvnorm"))
clusterEvalQ(cl, c(library(risksetROC),library(survival)))
clusterSetRNGStream(cl, 0120)
Cox.res.40 = pblapply(X = lapply(1:1000,function(x) c(500)),FUN = function(n){
  cox40.500 = res40.500 = c()
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  for (i in 1:length(pt)) {
    fit1 = risksetROC(Stime = Stime,status = status,marker = M,predict.time = pt[i],
                      method = "Cox",type = "b")
    fit2 = risksetROC(Stime = Stime,status = status,marker = M,predict.time = pt[i],
                      method = "Schoenfeld",type = "b",span = n^(-0.2))
    cox40.500[i] = fit1$AUC
    res40.500[i] = fit2$AUC
  }
  A = data.frame(Cox = cox40.500,
                 res = res40.500)
  return(A)
},cl = cl)

a = matrix(unlist(Cox.res.40,use.names = F),nrow = 1000,byrow = T)
cox40.500 = a[,1:9]
res40.500 = a[,10:18]
round(apply(cox40.500,2,mean),3)
round(apply(res40.500,2,mean),3)
round(apply(cox40.500,2,sd),3)
round(apply(res40.500,2,sd),3)


#PO
##20%,n=200,m=1000
n = 200;m = 1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 1.19;sdlog = 1

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-2)
clusterExport(cl, c("mean","weightedKM","EST.cindex.std","sigma","meanlog","sdlog",
                    "CoxWeights","MyrisksetAUC.wai","IntegrateAUC","POWeights",
                    "AFTWeights.1","prop.odds","Event","rmvnorm"))
clusterEvalQ(cl, c(library(muhaz),library(aftgee),library(survival),library(risksetROC),library(timereg)))
# 無信賴區間估計 ####
set.seed(2021)
clusterSetRNGStream(cl, 1)     
EST.PO.bn.20.200 = pblapply(X = lapply(1:1000,function(x) c(200)),FUN = function(n){
  AUC.bn.PO.np = c()
  s = rmvnorm(n,mean,sigma)
  t = exp(s[,2])
  M = s[,1]
  c = rlnorm(n, meanlog, sdlog)
  status <- ifelse(t<c,1,0)
  Stime = ifelse(t<c,t,c)
  data = data.frame(time=rep(0,length(Stime)),time2=Stime,status = status,eta = M)
  fit = prop.odds(Event(time2,status)~M,max.time=max(Stime),data = data,n.sim=100)
  gamma.t = as.numeric(fit$gamma)
  fit.po = MyrisksetAUC.wai( Stime , entry = NULL, status, gamma.t*M, "PO" 
                             , max(Stime) , AFT.method.1 = "1",AFT.method.2 = "1"
                             , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                             , weight = "rescale", plot=TRUE, type="l",xlab="Time", ylab="AUC")
  for( r in 1:9){
    AUC.bn.PO.np[r] = fit.po$AUC[which.min(abs(fit.po$utimes-exp(-2.5+0.5*r)))]
  }
  Cindex.bn.PO.np = fit.po$Cindex
  A = data.frame(EST = c(AUC.bn.PO.np,Cindex.bn.PO.np))
  return(A)
},cl = cl)

##將AUC<0.5的改成0.5
PO20 = matrix(unlist(EST.PO.bn.20.200,use.names = F),nrow = 1000,byrow = T)
PO20[which(PO20 < 0.5)] = 0.5

bn.PO.20.nci.200 = data.frame(EST = round(apply(PO20,2,mean),3))

bn.PO.SD.20.200 = sapply(1:10,function(i)sd(sapply(1:length(PO20[,1]),function(k) PO20[k,i])))
bn.PO.SD.20.200 = data.frame(SD = round(bn.PO.SD.20.200,3))
bn.PO.SD.20.200

##40%,n=200,m=1000
n = 200;m = 1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 0.36;sdlog = 1

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-2)
clusterExport(cl, c("mean","weightedKM","EST.cindex.std","sigma","meanlog","sdlog",
                    "CoxWeights","MyrisksetAUC.wai","IntegrateAUC","POWeights",
                    "AFTWeights.1","prop.odds","Event","rmvnorm"))
clusterEvalQ(cl, c(library(muhaz),library(aftgee),library(survival),library(risksetROC),library(timereg)))
# 無信賴區間估計 ####
set.seed(2021)
clusterSetRNGStream(cl, 1)
ANS.PO.bn.40.200 = pblapply(X = lapply(1:1000,function(x) c(200)),FUN = function(n){
  AUC.bn.PO.np = c()
  s = rmvnorm(n,mean,sigma)
  t = exp(s[,2])
  M = s[,1]
  c = rlnorm(n, meanlog, sdlog)
  status <- ifelse(t<c,1,0)
  Stime = ifelse(t<c,t,c)
  data = data.frame(time=rep(0,length(Stime)),time2=Stime,status = status,eta = M)
  fit = prop.odds(Event(time2,status)~eta,max.time=max(Stime),data = data,n.sim=100)
  gamma.t = as.numeric(fit$gamma)
  fit.po = MyrisksetAUC.wai( Stime , entry = NULL, status, gamma.t*M, "PO" 
                             , max(Stime) , AFT.method.1 = "1",AFT.method.2 = "1"
                             , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                             , weight = "rescale", plot=TRUE, type="l",xlab="Time", ylab="AUC")
  for( r in 1:9){
    AUC.bn.PO.np[r] = fit.po$AUC[which.min(abs(fit.po$utimes-exp(-2.5+0.5*r)))]
  }
  Cindex.bn.PO.np = fit.po$Cindex
  A = data.frame(EST = c(AUC.bn.PO.np,Cindex.bn.PO.np))
  return(A)
},cl = cl)  

##將AUC<0.5的改成0.5
PO40 = matrix(unlist(ANS.PO.bn.40.200,use.names = F),nrow = 1000,byrow = T)
PO40[which(PO40 < 0.5)] = 0.5

bn.PO.40.nci.200 = data.frame(EST = round(apply(PO40,2,mean),3))

bn.PO.SD.40.200 = sapply(1:10,function(i)sd(sapply(1:length(PO40[,1]),function(k) PO40[k,i])))
bn.PO.SD.40.200 = data.frame(SD = round(bn.PO.SD.40.200,3))
bn.PO.SD.40.200

#AFT
##20%,n=200,m=1000
n = 200;m = 1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 1.19;sdlog = 1

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-2)
clusterExport(cl, c("MyrisksetAUC.wai","AFTWeights.1","AFTWeights.2","EST.cindex.std",
                    "mean","sigma","meanlog","sdlog","rmvnorm"))
clusterEvalQ(cl, c(library(muhaz),library(aftgee),library(survival),library(risksetROC)))
# 無信賴區間估計 ####
set.seed(2021)
clusterSetRNGStream(cl, 1)     
ANS.AFT.bn.20.200 = pblapply(X = lapply(1:1000,function(x) c(200)),FUN = function(n){
  AUC.bn.AFT.np = AUC.bn.AFT.np.se = AUC.bn.AFT.np.cp = c()
  s = rmvnorm(n,mean,sigma)
  t = exp(s[,2])
  M = s[,1]
  c = rlnorm(n, meanlog, sdlog)
  Stime = ifelse(t<c,t,c)
  status = ifelse(t<c,1,0)
  beta.aft <- -aftsrr(Surv(Stime,status)~M)$beta
  fit.aft.np = MyrisksetAUC.wai(Stime , entry = NULL, status, beta.aft*M, "AFT" 
                                , max(Stime) , AFT.method.1 = "1", AFT.method.2 = "1"
                                , bw.1="1/2",bw.2 = "1/2", bw.mul.1=10, bw.mul.2=0.5, CI=F, confi=0.95
                                , weight = "rescale", plot=TRUE, type="l",xlab="Time", ylab="AUC")
  for( r in 1:9){
    AUC.bn.AFT.np[r] = fit.aft.np$AUC[which.min(abs(fit.aft.np$utimes-exp(-2.5+0.5*r)))]
  }
  Cindex.bn.AFT.np = fit.aft.np$Cindex
  A = data.frame(EST = c(AUC.bn.AFT.np,Cindex.bn.AFT.np))
  return(A)
},cl = cl)

##將AUC<0.5的改成0.5
AFT20 = matrix(unlist(ANS.AFT.bn.20.200,use.names = F),nrow = 1000,byrow = T)
AFT20[which(AFT20 < 0.5)] = 0.5

bn.AFT.20.nci.200 = data.frame(EST = round(apply(AFT20,2,mean),3))

bn.AFT.SD.20.200 = sapply(1:10,function(i)sd(sapply(1:length(AFT20[,1]),function(k) AFT20[k,i])))
bn.AFT.SD.20.200 = data.frame(SD = round(bn.AFT.SD.20.200,3))
bn.AFT.SD.20.200
  
##40%,n=200,m=1000
n = 200;m = 1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 0.36;sdlog = 1

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-2)
clusterExport(cl, c("MyrisksetAUC.wai","AFTWeights.1","AFTWeights.2","EST.cindex.std",
                    "mean","sigma","meanlog","sdlog","rmvnorm"))
clusterEvalQ(cl, c(library(muhaz),library(aftgee),library(survival),library(risksetROC)))
# 無信賴區間估計 ####
set.seed(2021)
clusterSetRNGStream(cl, 1)     
ANS.AFT.bn.40.200 = pblapply(X = lapply(1:1000,function(x) c(200)),FUN = function(n){
  AUC.bn.AFT.np = AUC.bn.AFT.np.se = AUC.bn.AFT.np.cp = c()
  s = rmvnorm(n,mean,sigma)
  t = exp(s[,2])
  M = s[,1]
  c = rlnorm(n, meanlog, sdlog)
  Stime = ifelse(t<c,t,c)
  status = ifelse(t<c,1,0)
  beta.aft <- -aftsrr(Surv(Stime,status)~M)$beta
  fit.aft.np = MyrisksetAUC.wai(Stime , entry = NULL, status, beta.aft*M, "AFT" 
                                , max(Stime) , AFT.method.1 = "1", AFT.method.2 = "1"
                                , bw.1="1/2",bw.2 = "1/2", bw.mul.1=10, bw.mul.2=0.5, CI=F, confi=0.95
                                , weight = "rescale", plot=TRUE, type="l",xlab="Time", ylab="AUC")
  for( r in 1:9){
    AUC.bn.AFT.np[r] = fit.aft.np$AUC[which.min(abs(fit.aft.np$utimes-exp(-2.5+0.5*r)))]
  }
  Cindex.bn.AFT.np = fit.aft.np$Cindex
  A = data.frame(EST = c(AUC.bn.AFT.np,Cindex.bn.AFT.np))
  return(A)
},cl = cl)
  
#將AUC<0.5的改成0.5
AFT40 = matrix(unlist(ANS.AFT.bn.40.200,use.names = F),nrow = 1000,byrow = T)
AFT40[which(AFT40 < 0.5)] = 0.5

bn.AFT.40.nci.200 = data.frame(EST = round(apply(AFT40,2,mean),3))

bn.AFT.SD.40.200 = sapply(1:10,function(i)sd(sapply(1:length(AFT40[,1]),function(k) AFT40[k,i])))
bn.AFT.SD.40.200 = data.frame(SD = round(bn.AFT.SD.40.200,3))
bn.AFT.SD.40.200

##method = KM,CD1
##20%
set.seed(2021)
n = 200;m = 1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 1.19;sdlog = 1
pt = exp(seq(-2,2,0.5))
CD1.20 = matrix(rep(0,length(pt)*m),nrow = m,byrow = T)
cen = c()
for (j in 1:m) {
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  cen[j] = table(status)[1]/sum(table(status))
  for (i in 1:length(pt)) {
    fit_C1 = survivalROC(Stime = Stime,
                         status = status,
                         marker = M, 
                         predict.time = pt[i], ##time points of the ROC curve
                         method = "KM")
    CD1.20[j,i] = fit_C1$AUC
  }
}
mean(cen)
which(is.na(CD1.20))
which(CD1.20<0.5)

##將AUC=NA的改成0.5
CD1.20[which(is.na(CD1.20))] = 0.5
##將AUC<0.5的改成0.5
CD1.20[which(CD1.20<0.5)] = 0.5

bn.CD1.20.nci.200 = data.frame(EST = round(apply(CD1.20,2,mean),3))

bn.CD1.SD.20.200 = sapply(1:9,function(i)sd(sapply(1:length(CD1.20[,1]),function(k) CD1.20[k,i])))
bn.CD1.SD.20.200 = data.frame(SD = round(bn.CD1.SD.20.200,3))
bn.CD1.SD.20.200

##40%
set.seed(2021)
n = 200;m = 1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 0.36;sdlog = 1
pt = exp(seq(-2,2,0.5))
CD1.40 = matrix(rep(0,length(pt)*m),nrow = m,byrow = T)
cen = c()
for (j in 1:m) {
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  cen[j] = table(status)[1]/sum(table(status))
  for (i in 1:length(pt)) {
    fit_C1 = survivalROC(Stime = Stime,
                         status = status,
                         marker = M, 
                         predict.time = pt[i], ##time points of the ROC curve
                         method = "KM")
    CD1.40[j,i] = fit_C1$AUC
  }
}
mean(cen)
which(is.na(CD1.40))
which(CD1.40<0.5)

##將AUC=NA的改成0.5
CD1.40[which(is.na(CD1.40))] = 0.5
##將AUC<0.5的改成0.5
CD1.40[which(CD1.40<0.5)] = 0.5

bn.CD1.40.nci.200 = data.frame(EST = round(apply(CD1.40,2,mean),3))

bn.CD1.SD.40.200 = sapply(1:9,function(i)sd(sapply(1:length(CD1.40[,1]),function(k) CD1.40[k,i])))
bn.CD1.SD.40.200 = data.frame(SD = round(bn.CD1.SD.40.200,3))
bn.CD1.SD.40.200

#method = NNE,CD2
#20%
set.seed(2021)
n = 200;m = 1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 1.19;sdlog = 1
pt = exp(seq(-2,2,0.5))
CD2.20 = matrix(rep(0,length(pt)*m),nrow = m,byrow = T)

for (j in 1:m) {
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  for (i in 1:length(pt)) {
    fit_C2 = survivalROC.fix(Stime = Stime,
                             status = status,
                             marker = M, 
                             predict.time = pt[i],
                             method = "NNE",
                             span = n^(-0.2)) ##Span for the NNE(or lambda)
    CD2.20[j,i] = fit_C2$AUC
  }
}
which(is.na(CD2.20))
which(CD2.20<0.5)

#survival or censoring time > predict.time
survivalROC
##將AUC=NA的改成0.5
CD2.20[which(is.na(CD2.20))] = 0.5
##將AUC<0.5的改成0.5
CD2.20[which(CD2.20<0.5)] = 0.5

bn.CD2.20.nci.200 = data.frame(EST = round(apply(CD2.20,2,mean),3))

bn.CD2.SD.20.200 = sapply(1:9,function(i)sd(sapply(1:length(CD2.20[,1]),function(k) CD2.20[k,i])))
bn.CD2.SD.20.200 = data.frame(SD = round(bn.CD2.SD.20.200,3))
bn.CD2.SD.20.200


##40%
set.seed(2021)
n = 200;m = 1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 0.36;sdlog = 1
pt = exp(seq(-2,2,0.5))
CD2.40 = matrix(rep(0,length(pt)*m),nrow = m,byrow = T)
cen = c()
for (j in 1:m) {
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  cen[j] = table(status)[1]/sum(table(status))
  for (i in 1:length(pt)) {
    fit_C2 = survivalROC.fix(Stime = Stime,
                             status = status,
                             marker = M, 
                             predict.time = pt[i],
                             method = "NNE",
                             span = n^(-0.2)) ##Span for the NNE(or lambda)
    CD2.40[j,i] = fit_C2$AUC
  }
}
mean(cen)
which(is.na(CD2.40))
which(CD2.40<0.5)

##將AUC=NA的改成0.5
CD2.40[which(is.na(CD2.40))] = 0.5
##將AUC<0.5的改成0.5
CD2.40[which(CD2.40<0.5)] = 0.5

bn.CD2.40.nci.200 = data.frame(EST = round(apply(CD2.40,2,mean),3))

bn.CD2.SD.40.200 = sapply(1:9,function(i)sd(sapply(1:length(CD2.40[,1]),function(k) CD2.40[k,i])))
bn.CD2.SD.40.200 = data.frame(SD = round(bn.CD2.SD.40.200,3))
bn.CD2.SD.40.200

##檢查為何有NA
set.seed(0120)
s = rmvnorm(n,mean,sigma)
T = exp(s[,2])
M = s[,1]
C = rlnorm(n, meanlog, sdlog)
status = ifelse(T<C,1,0)
Stime = ifelse(T<C,T,C)
abc = survivalROC(Stime = Stime,
                  status = status,
                  marker = M, 
                  predict.time = pt[1],
                  method = "NNE",
                  span = n^(-0.2)) ##Span for the NNE(or lambda)
abc$AUC

##解決因沒有event發生導致fun.無法執行的問題
survivalROC.fix = function (Stime, status, marker, entry = NULL, predict.time, 
                            cut.values = NULL, method = "NNE", lambda = NULL, span = NULL, 
                            window = "symmetric") 
{
  times = Stime
  x <- marker
  if (is.null(entry)) 
    entry <- rep(0, length(times))
  bad <- is.na(times) | is.na(status) | is.na(x) | is.na(entry)
  entry <- entry[!bad]
  times <- times[!bad]
  status <- status[!bad]
  x <- x[!bad]
  if (sum(bad) > 0) 
    cat(paste("\n", sum(bad), "records with missing values dropped. \n"))
  if (is.null(cut.values)) 
    cut.values <- unique(x)
  cut.values <- cut.values[order(cut.values)]
  ncuts <- length(cut.values)
  ooo <- order(times)
  times <- times[ooo]
  status <- status[ooo]
  x <- x[ooo]
  s0 <- 1
  unique.t0 <- unique(times)
  unique.t0 <- unique.t0[order(unique.t0)]
  
  n.times <- sum(unique.t0 <= predict.time)
  for (j in 1:n.times) {
    n <- sum(entry <= unique.t0[j] & times >= unique.t0[j])
    d <- sum((entry <= unique.t0[j]) & (times == unique.t0[j]) & 
               (status == 1))
    if (n > 0) 
      s0 <- s0 * (1 - d/n)
  }
  s.pooled <- s0
  roc.matrix <- matrix(NA, ncuts, 2)
  roc.matrix[ncuts, 1] <- 0
  roc.matrix[ncuts, 2] <- 1
  if (method == "KM") {
    for (c in 1:(ncuts - 1)) {
      s0 <- 1
      subset <- as.logical(x > cut.values[c])
      e0 <- entry[subset]
      t0 <- times[subset]
      c0 <- status[subset]
      if (!is.null(t0)) {
        unique.t0 <- unique(t0)
        unique.t0 <- unique.t0[order(unique.t0)]
        n.times <- sum(unique.t0 <= predict.time)
        if (n.times > 0) {
          for (j in 1:n.times) {
            n <- sum(e0 <= unique.t0[j] & t0 >= unique.t0[j])
            d <- sum((e0 <= unique.t0[j]) & (t0 == unique.t0[j]) & 
                       (c0 == 1))
            if (n > 0) 
              s0 <- s0 * (1 - d/n)
          }
        }
      }
      p0 <- mean(subset)
      roc.matrix[c, 1] <- (1 - s0) * p0/(1 - s.pooled)
      roc.matrix[c, 2] <- 1 - s0 * p0/s.pooled
    }
  }
  if (method == "NNE") {
    if (is.null(lambda) & is.null(span)) {
      cat("method = NNE requires either lambda or span! \n")
      stop(0)
    }
    x.unique <- unique(x)
    x.unique <- x.unique[order(x.unique)]
    S.t.x <- rep(0, length(x.unique))
    t.evaluate <- unique(times[status == 1])
    t.evaluate <- t.evaluate[order(t.evaluate)]
    t.evaluate <- t.evaluate[t.evaluate <= predict.time]
    if(length(t.evaluate) == 0){
      predict.time = pt[2]
      t.evaluate <- unique(times[status == 1])
      t.evaluate <- t.evaluate[order(t.evaluate)]
      t.evaluate <- t.evaluate[t.evaluate <= predict.time]
    }
    for (j in 1:length(x.unique)) {
      if (!is.null(span)) {
        if (window == "symmetric") {
          ddd <- (x - x.unique[j])
          n <- length(x)
          ddd <- ddd[order(ddd)]
          index0 <- sum(ddd < 0) + 1
          index1 <- index0 + trunc(n * span + 0.5)
          if (index1 > n) 
            index1 <- n
          lambda <- ddd[index1]
          wt <- as.integer(((x - x.unique[j]) <= lambda) & 
                             ((x - x.unique[j]) >= 0))
          index0 <- sum(ddd <= 0)
          index2 <- index0 - trunc(n * span/2)
          if (index2 < 1) 
            index2 <- 1
          lambda <- abs(ddd[index1])
          set.index <- ((x - x.unique[j]) >= -lambda) & 
            ((x - x.unique[j]) <= 0)
          wt[set.index] <- 1
        }
        if (window == "asymmetric") {
          ddd <- (x - x.unique[j])
          n <- length(x)
          ddd <- ddd[order(ddd)]
          index0 <- sum(ddd < 0) + 1
          index <- index0 + trunc(n * span)
          if (index > n) 
            index <- n
          lambda <- ddd[index]
          wt <- as.integer(((x - x.unique[j]) <= lambda) & 
                             ((x - x.unique[j]) >= 0))
        }
      }
      else {
        wt <- exp(-(x - x.unique[j])^2/lambda^2)
      }
      s0 <- 1
      for (k in 1:length(t.evaluate)) {
        n <- sum(wt * (entry <= t.evaluate[k]) & (times >= 
                                                    t.evaluate[k]))
        d <- sum(wt * (entry <= t.evaluate[k]) & (times == 
                                                    t.evaluate[k]) * (status == 1))
        if (n > 0) 
          s0 <- s0 * (1 - d/n)
      }
      S.t.x[j] <- s0
    }
    S.all.x <- S.t.x[match(x, x.unique)]
    n <- length(times)
    S.marginal <- sum(S.all.x)/n
    for (c in 1:(ncuts - 1)) {
      p1 <- sum(x > cut.values[c])/n
      Sx <- sum(S.all.x[x > cut.values[c]])/n
      roc.matrix[c, 1] <- (p1 - Sx)/(1 - S.marginal)
      roc.matrix[c, 2] <- 1 - Sx/S.marginal
    }
  }
  sensitivity = roc.matrix[, 1]
  specificity = roc.matrix[, 2]
  x <- 1 - c(0, specificity)
  y <- c(1, sensitivity)
  n <- length(x)
  dx <- x[-n] - x[-1]
  mid.y <- (y[-n] + y[-1])/2
  area <- sum(dx * mid.y)
  list(cut.values = c(-Inf, cut.values), TP = y, FP = x, predict.time = predict.time, 
       Survival = s.pooled, AUC = area)
}


##CD5,with a Cox model for computing the weights.
##20%
set.seed(2021)
n = 200;m = 1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 1.19;sdlog = 1
pt = exp(seq(-2,2,0.5))
CD5.20 = matrix(rep(0,length(pt)*m),nrow = m,byrow = T)
#CD6.20 = matrix(rep(0,length(pt)*m),nrow = m,byrow = T)
cen = c()
for (j in 1:m) {
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  cen[j] = table(status)[1]/sum(table(status))
  fit_C5 = timeROC(T = Stime,delta = status,marker = M,cause = 1,
                   weighting = "cox",times = pt,ROC = F)
  #fit_C6 = timeROC(T = Stime,delta = status,marker = M,cause = 1,
  #                 weighting = "marginal",times = pt,ROC = F,iid = TRUE)
  CD5.20[j,] = fit_C5$AUC
  #CD6.20[j,] = fit_C6$AUC
}
mean(cen)
CD5.20
which(CD5.20<0.5)
which(is.na(CD5.20))
#CD6.20

##將AUC=NA的改成0.5
CD5.20[which(is.na(CD5.20))] = 0.5
##將AUC<0.5的改成0.5
CD5.20[which(CD5.20<0.5)] = 0.5

bn.CD5.20.nci.200 = data.frame(EST = round(apply(CD5.20,2,mean),3))
#bn.CD6.20.nci.200 = data.frame(EST = round(apply(C6,2,mean),3))

bn.CD5.SD.20.200 = sapply(1:9,function(i)sd(sapply(1:length(CD5.20[,1]),function(k) CD5.20[k,i])))
bn.CD5.SD.20.200 = data.frame(SD = round(bn.CD5.SD.20.200,3))
bn.CD5.SD.20.200

##40%
set.seed(2021)
n = 200;m = 1000
mean = c(0,0);sigma = matrix(c(1,-0.7,-0.7,1),byrow = T,nrow = 2)
meanlog = 0.36;sdlog = 1
pt = exp(seq(-2,2,0.5))
CD5.40 = matrix(rep(0,length(pt)*m),nrow = m,byrow = T)
##CD6.40 = matrix(rep(0,length(pt)*m),nrow = m,byrow = T)
cen = c()
for (j in 1:m) {
  s = rmvnorm(n,mean,sigma)
  T = exp(s[,2])
  M = s[,1]
  C = rlnorm(n, meanlog, sdlog)
  status = ifelse(T<C,1,0)
  Stime = ifelse(T<C,T,C)
  cen[j] = table(status)[1]/sum(table(status))
  fit_C5 = timeROC(T = Stime,delta = status,marker = M,cause = 1,
                   weighting = "cox",times = pt,ROC = F)
  ##fit_C6 = timeROC(T = Stime,delta = status,marker = M,cause = 1,
  ##                 weighting = "marginal",times = pt,ROC = F,iid = TRUE)
  CD5.40[j,] = fit_C5$AUC
  ##CD6.40[j,] = fit_C6$AUC
}
mean(cen)
CD5.40
which(CD5.40<0.5)
##CD6.40
timeROC
##將AUC=NA的改成0.5
CD5.40[which(is.na(CD5.40))] = 0.5
##將AUC<0.5的改成0.5
CD5.40[which(CD5.40<0.5)] = 0.5


bn.CD5.40.nci.200 = data.frame(EST = round(apply(CD5.40,2,mean),3))
#bn.CD6.40.nci.200 = data.frame(EST = round(apply(CD6.40,2,mean),3))

bn.CD5.SD.40.200 = sapply(1:9,function(i)sd(sapply(1:length(CD5.40[,1]),function(k) CD5.40[k,i])))
bn.CD5.SD.40.200 = data.frame(SD = round(bn.CD5.SD.40.200,3))
bn.CD5.SD.40.200
cbind(bn.CD5.40.nci.200[,1],bn.CD5.SD.40.200)

##CD6,with the Kaplan-Meier estimator for computing the weights.
ovarian_CD6<-timeROC(T = ovarian$futime,delta = ovarian$fustat,
                     marker = ovarian$age,cause = 1,
                     weighting = "marginal",
                     times = c(365,399,476,770,855),ROC = TRUE,iid = TRUE)
ovarian_CD6
plot(ovarian_CD6,time = 365)        

  
  
------------------------------------------------------------------------------------

a=risksetAUC(Stime = Stime,status = status,marker = M,predict.time = 2,
               method = "Cox",type = "b",tmax = 7.39)
CoxWeights(marker = M,Stime = Stime,status = status,predict.time = 2)

coxph(Surv(Stime,status)~M)
??coxphw
?coxph
?risksetROC
?risksetAUC
?SchoenSmooth
?CoxWeights
??survivalROC
llCoxReg(Stime = Stime,status = status,marker = M,span = 0.4)
summary(fit)
fit$AUC
exp(2)
train = cbind(Stime = Stime[1:250],status = status[1:250],M = M[1:250])
test = cbind(Stime = Stime[251:500],status = status[251:500],M = M[251:500])
test = data.frame(test)
train = data.frame(train)
train.fit <- coxph(Surv(Stime, status) ~ M,x=TRUE, y=TRUE, method="breslow",data = train)

lp <- predict(train.fit)
lpnew <- predict(train.fit, newdata = test)
Surv.rsp <- Surv(train$Stime, train$status)
Surv.rsp.new <- Surv(test$Stime, test$status)
times <- seq(0, 7.38, 0.1)                  
exp(2)
AUC_CD <- AUC.cd(Surv.rsp, Surv.rsp.new, lp, lpnew, times)
?AUC.cd
##iauc = 0.8728364 (The summary measure of AUC.)

plot(AUC_CD$times,AUC_CD$auc,type = "l",main = "CD4 under Cox model \n data =  ovarian",
     xlab = "times",ylab = "AUC",las = 1)






