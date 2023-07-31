library(parallel)
library(pbapply)
library(Metrics)
alpha = 1.5 ; lambda = 0.25 ; beta = 1

AUC.Weibull.PH.true = apply(cox.wei,2,mean) #真實值計算.R中的結果
Cindex.Weibull.PH.true = mean(PH.wei.C)

cpu.cores <- detectCores() #計算CPU個數
cl <- makeCluster(cpu.cores-1) #最後可選擇-1或-2，依照運算速度決定
clusterExport(cl, c("beta","alpha","lambda","MyrisksetAUC.wai.new","MyrisksetAUC.wai.old","EST.cindex.std",
                    "CoxWeights","AUC.Weibull.PH.true","Cindex.Weibull.PH.true","CoxWeights.new")) #需放入會用到的fun.跟參數
clusterEvalQ(cl, c(library(aftgee),library(survival),library(risksetROC),library(timereg))) #放入會用到的package

clusterSetRNGStream(cl, 01205) #seed    
ANS.PH.wei.20.200 = pblapply(X = lapply(1:500,function(x) c(200)),FUN = function(n){
  AUC.Weibull.PH.np.old = AUC.Weibull.PH.np.new = AUC.Weibull.PH.np.se = AUC.Weibull.PH.np.cp = c()
  z = rnorm(n)
  c = rexp(n,0.10) #0.10 <- 20% 0.27 <- 40%
  u = runif(n,0,1)
  t = (-log(u)/(lambda*exp(beta*z)))^(1/alpha)
  pt = exp(seq(-3,3,0.5))
  status = ifelse(t<c,1,0)
  Stime = ifelse(t<c,t,c)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  beta.ph <- coxph(Surv(Stime,status)~z)$coefficients
  fit.ph.old = MyrisksetAUC.wai.old( Stime , entry = NULL, status, beta.ph*z, "Cox"
                                     , tmax = max(utimes) , AFT.method.1 = "1",AFT.method.2 = "1"
                                     , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                     , weight = "rescale", plot=TRUE, type="l",xlab="Time", ylab="AUC")
  fit.ph.new = MyrisksetAUC.wai.new( Stime , entry = NULL, status, beta.ph*z, "Cox"
                                     , tmax = max(utimes) , AFT.method.1 = "1",AFT.method.2 = "1"
                                     , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                     , weight = "rescale", plot=TRUE, type="l",xlab="Time", ylab="AUC")
  for( r in 1:13){
    AUC.Weibull.PH.np.old[r] = fit.ph.old$AUC[which.min(abs(fit.ph.old$utimes-exp(-3.5+0.5*r)))]
    AUC.Weibull.PH.np.new[r] = fit.ph.new$AUC[which.min(abs(fit.ph.new$utimes-exp(-3.5+0.5*r)))]
  }
  Cindex.Weibull.PH.np.old = fit.ph.old$Cindex.old
  Cindex.Weibull.PH.np.new = fit.ph.new$Cindex.new
  A = data.frame(EST = c(AUC.Weibull.PH.np.old,Cindex.Weibull.PH.np.old),
                 Modify = c(AUC.Weibull.PH.np.new,Cindex.Weibull.PH.np.new))
  return(A)
},cl = cl)

PH.wei.20.200 = matrix(unlist(ANS.PH.wei.20.200,use.names = F),nrow = 500,byrow = T)
PH.wei.20.200.old = PH.wei.20.200[,1:14]
PH.wei.20.200.new = PH.wei.20.200[,15:28]

#計算MSE
mse.wei.PH.old.20 = c()
mse.wei.PH.new.20 = c()
for (i in 1:length(AUC.Weibull.PH.true)) {
  mse.wei.PH.old.20[i] = mse(AUC.Weibull.PH.true[i],PH.wei.20.200.old[,i])
  mse.wei.PH.new.20[i] = mse(AUC.Weibull.PH.true[i],PH.wei.20.200.new[,i])
}
plot(seq(-3,3,0.5), mse.wei.PH.old.20, type="l", xlab="Log time", ylab="MSE",lwd = 2,col = "red",main = "Weibull-Cox for 20% censoring")
lines(seq(-3,3,0.5), mse.wei.PH.new.20, type="l",lwd = 2,col = "blue")
legend("topleft", c("original","Modify"), cex=1, fill=c("red","blue"))

which(PH.wei.20.200.old < 0.5)
which(PH.wei.20.200.new < 0.5)
wei.PH.20.nci.200 = data.frame(Logtime = c(seq(-3,3,0.5),"C"),
                               True = c(round(AUC.Weibull.PH.true,3),round(Cindex.Weibull.PH.true,3)),
                               EST = round(apply(PH.wei.20.200.old,2,mean),3),
                               SD.old = round(apply(PH.wei.20.200.old,2,sd),3),
                               MSE.old = round(c(10*mse.wei.PH.old.20,10*mse(PH.wei.20.200.old[,14],Cindex.Weibull.PH.true)),3),
                               Modify = round(apply(PH.wei.20.200.new,2,mean),3),
                               SD.new = round(apply(PH.wei.20.200.new,2,sd),3),
                               MSE.new = round(c(10*mse.wei.PH.new.20,10*mse(PH.wei.20.200.new[,14],Cindex.Weibull.PH.true)),3))
wei.PH.20.nci.200
wei.PH.20.nci.200[c(2:12,14),]

#畫圖
plot(wei.PH.20.nci.200$Logtime[-14],wei.PH.20.nci.200$EST[-14],type = "l", xlab="Log time", ylab="AUC",lwd = 2,ylim=c(0.2,1.0),col = "blue",main = "Weibull-Cox for 20% censoring")
lines(wei.PH.20.nci.200$Logtime[-14], wei.PH.20.nci.200$Modify[-14], type="l",lwd = 2,col = "red")
lines(wei.PH.20.nci.200$Logtime[-14],wei.PH.20.nci.200$True[-14], type="l",lwd = 2,col = "black")
legend("topright", c("original","Modify","Target"), cex=1, fill=c("blue","red","black"))
abline(h=0.5)

#40%
clusterSetRNGStream(cl, 01205)     
ANS.PH.wei.40.200 = pblapply(X = lapply(1:500,function(x) c(200)),FUN = function(n){
  AUC.Weibull.PH.np.old = AUC.Weibull.PH.np.new = c()
  z = rnorm(n)
  c = rexp(n,0.27) #0.10 <- 20% 0.27 <- 40%
  u = runif(n,0,1)
  t = (-log(u)/(lambda*exp(beta*z)))^(1/alpha)
  pt = exp(seq(-3,3,0.5)) 
  status = ifelse(t<c,1,0)
  Stime = ifelse(t<c,t,c)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  beta.ph <- coxph(Surv(Stime,status)~z)$coefficients
  fit.ph.new = MyrisksetAUC.wai.new( Stime , entry = NULL, status, beta.ph*z, "Cox"
                                     , tmax = max(pt) , AFT.method.1 = "1",AFT.method.2 = "1"
                                     , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                     , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  fit.ph.old = MyrisksetAUC.wai.old( Stime , entry = NULL, status, beta.ph*z, "Cox"
                                     , tmax = max(pt) , AFT.method.1 = "1",AFT.method.2 = "1"
                                     , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                     , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  for( r in 1:13){
    AUC.Weibull.PH.np.old[r] = fit.ph.old$AUC[which.min(abs(fit.ph.old$utimes-exp(-3.5+0.5*r)))]
    AUC.Weibull.PH.np.new[r] = fit.ph.new$AUC[which.min(abs(fit.ph.new$utimes-exp(-3.5+0.5*r)))]
  }
  Cindex.Weibull.PH.np.old = fit.ph.old$Cindex.old
  Cindex.Weibull.PH.np.new = fit.ph.new$Cindex.new
  
  A = data.frame(EST = c(AUC.Weibull.PH.np.old,Cindex.Weibull.PH.np.old),
                 Modify = c(AUC.Weibull.PH.np.new,Cindex.Weibull.PH.np.new))
  return(A)
},cl = cl)

PH.wei.40.200 = matrix(unlist(ANS.PH.wei.40.200,use.names = F),nrow = 500,byrow = T)
PH.wei.40.200.old = PH.wei.40.200[,1:14]
PH.wei.40.200.new = PH.wei.40.200[,15:28]

mse.wei.PH.old.40 = c()
mse.wei.PH.new.40 = c()
for (i in 1:length(AUC.Weibull.PH.true)) {
  mse.wei.PH.old.40[i] = mse(AUC.Weibull.PH.true[i],PH.wei.40.200.old[,i])
  mse.wei.PH.new.40[i] = mse(AUC.Weibull.PH.true[i],PH.wei.40.200.new[,i])
}
plot(seq(-3,3,0.5), mse.wei.PH.old.40, type="l", xlab="Log time", ylab="MSE",lwd = 2,col = "red",main = "Weibull-Cox for 40% censoring")
lines(seq(-3,3,0.5), mse.wei.PH.new.40, type="l",lwd = 2,col = "blue")
legend("topleft", c("original","Modify"), cex=1, fill=c("red","blue"))

which(PH.wei.40.200.old < 0.5)
which(PH.wei.40.200.new < 0.5)
wei.PH.40.nci.200 = data.frame(Logtime = c(seq(-3,3,0.5),"C"),
                               True = c(round(AUC.Weibull.PH.true,3),round(Cindex.Weibull.PH.true,3)),
                               EST = round(apply(PH.wei.40.200.old,2,mean),3),
                               SD.old = round(apply(PH.wei.40.200.old,2,sd),3),
                               MSE.old = round(c(10*mse.wei.PH.old.40,10*mse(PH.wei.40.200.old[,14],Cindex.Weibull.PH.true)),3),
                               Modify = round(apply(PH.wei.40.200.new,2,mean),3),
                               SD.new = round(apply(PH.wei.40.200.new,2,sd),3),
                               MSE.new = round(c(10*mse.wei.PH.new.40,10*mse(PH.wei.40.200.new[,14],Cindex.Weibull.PH.true)),3))
wei.PH.40.nci.200
wei.PH.40.nci.200[c(2:12,14),]
plot(wei.PH.40.nci.200$Logtime[-14],wei.PH.40.nci.200$EST[-14],type = "l", xlab="Log time", ylab="AUC",lwd = 2,ylim=c(0.2,1.0),col = "blue",main = "Weibull-Cox for 40% censoring")
lines(wei.PH.40.nci.200$Logtime[-14], wei.PH.40.nci.200$Modify[-14], type="l",lwd = 2,col = "red")
lines(wei.PH.40.nci.200$Logtime[-14],wei.PH.40.nci.200$True[-14], type="l",lwd = 2,col = "black")
legend("topright", c("original","Modify","Target"), cex=1, fill=c("blue","red","black"))
abline(h=0.5)


stopCluster(cl)
