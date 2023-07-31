library(parallel)
library(pbapply)
alpha = 1.5 ; lambda = 0.05 ; beta = 1.0

AUC.Weibull.AFT.TRUE = apply(AFT.wei,2,mean)
Cindex.Weibull.AFT.TRUE = mean(AFT.wei.C)

cpu.cores.1 <- detectCores()
cl <- makeCluster(cpu.cores.1-2)
clusterExport(cl, c("MyrisksetAUC.wai.new","MyrisksetAUC.wai.old","AFTWeights.1","AFTWeights.2","EST.cindex.std",
                    "lambda","alpha","beta","AFTWeights.1.new","AFTWeights.2.new"))
clusterEvalQ(cl, c(library(muhaz),library(aftgee),library(survival),library(risksetROC)))
# 有信賴區間估計 ####
clusterSetRNGStream(cl, 0120)     
ANS.AFT.wei.20.200 = pblapply(X = lapply(1:500,function(x) c(200)),FUN = function(n){
  AUC.Weibull.AFT.np.old = AUC.Weibull.AFT.np.new = c()
  z = rnorm(n)
  u = runif(n,0,1)
  t = (-log(u)/(lambda*exp(alpha*beta*z)))^(1/alpha)
  c = rexp(n,0.02)  #設限率20% <- 0.02  40% <- 0.07
  pt = exp(seq(-5,5,0.25))
  Stime = ifelse(t<c,t,c)
  status = ifelse(t<c,1,0)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  beta.aft <- -aftsrr(Surv(Stime,status)~z)$beta
  fit.aft.np.old = MyrisksetAUC.wai.old(Stime , entry = NULL, status, beta.aft*z, "AFT"
                                        , max(pt) , AFT.method.1 = "1", AFT.method.2 = "1"
                                        , bw.1="1/2",bw.2 = "1/2", bw.mul.1=10, bw.mul.2=0.5, CI=F, confi=0.95
                                        , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  fit.aft.np.new = MyrisksetAUC.wai.new(Stime , entry = NULL, status, beta.aft*z, "AFT"
                                        , max(pt) , AFT.method.1 = "1", AFT.method.2 = "1"
                                        , bw.1="1/2",bw.2 = "1/2", bw.mul.1=10, bw.mul.2=0.5, CI=F, confi=0.95
                                        , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  for( r in 1:41){
    AUC.Weibull.AFT.np.old[r] = fit.aft.np.old$AUC[which.min(abs(fit.aft.np.old$utimes-exp(-5.25+0.25*r)))]
    AUC.Weibull.AFT.np.new[r] = fit.aft.np.new$AUC[which.min(abs(fit.aft.np.new$utimes-exp(-5.25+0.25*r)))]
  }
  Cindex.Weibull.AFT.np.old = fit.aft.np.old$Cindex.old
  Cindex.Weibull.AFT.np.new = fit.aft.np.new$Cindex.new
  
  A = data.frame(EST = c(AUC.Weibull.AFT.np.old,Cindex.Weibull.AFT.np.old),
                 Modify = c(AUC.Weibull.AFT.np.new,Cindex.Weibull.AFT.np.new))
  return(A)
},cl = cl)

AFT.wei.20.200 = matrix(unlist(ANS.AFT.wei.20.200,use.names = F),nrow = 500,byrow = T)
AFT.wei.20.200.old = AFT.wei.20.200[,1:42]
AFT.wei.20.200.new = AFT.wei.20.200[,43:84]

#計算MSE
mse.wei.AFT.old.20 = c()
mse.wei.AFT.new.20 = c()
for (i in 1:length(AUC.Weibull.AFT.TRUE)) {
  mse.wei.AFT.old.20[i] = mse(AUC.Weibull.AFT.TRUE[i],AFT.wei.20.200.old[,i])
  mse.wei.AFT.new.20[i] = mse(AUC.Weibull.AFT.TRUE[i],AFT.wei.20.200.new[,i])
}
plot(seq(-5,5,0.25), mse.wei.AFT.old.20, type="l", xlab="Log time", ylab="MSE",lwd = 2,col = "red",main = "Weibull-AFT for 20% censoring")
lines(seq(-5,5,0.25), mse.wei.AFT.new.20, type="l",lwd = 2,col = "blue")
legend("topleft", c("original","Modify"), cex=1, fill=c("red","blue"))

which(AFT.wei.20.200.old < 0.5)
which(AFT.wei.20.200.new < 0.5)

wei.AFT.20.nci.200 = data.frame(Logtime = c(seq(-5,5,0.25),"C"),
                                True = round(c(AUC.Weibull.AFT.TRUE,Cindex.Weibull.AFT.TRUE),3),
                                EST = round(apply(AFT.wei.20.200.old,2,mean),3),
                                SD.old = round(apply(AFT.wei.20.200.old,2,sd),3),
                                MSE.old = round(c(10*mse.wei.AFT.old.20,10*mse(AFT.wei.20.200.old[,42],Cindex.Weibull.AFT.TRUE)),3),
                                Modify = round(apply(AFT.wei.20.200.new,2,mean),3),
                                SD.new = round(apply(AFT.wei.20.200.new,2,sd),3),
                                MSE.new = round(c(10*mse.wei.AFT.new.20,10*mse(AFT.wei.20.200.new[,42],Cindex.Weibull.AFT.TRUE)),3))
wei.AFT.20.nci.200
wei.AFT.20.nci.200[c(17,21,25,29,33,35,37,38,39,42),]

#畫圖
plot(wei.AFT.20.nci.200$Logtime[-c(42)],wei.AFT.20.nci.200$EST[-c(42)],type = "l", xlab="Log time", ylab="AUC",lwd = 2,ylim=c(0.2,1.0),col = "blue",main = "Weibull-AFT for 20% censoring")
lines(wei.AFT.20.nci.200$Logtime[-c(42)],wei.AFT.20.nci.200$Modify[-c(42)], type="l",lwd = 2,col = "red")
lines(wei.AFT.20.nci.200$Logtime[-c(42)],wei.AFT.20.nci.200$True[-c(42)], type="l",lwd = 2,col = "black")
legend("topright", c("original","Modify","Target"), cex=1, fill=c("blue","red","black"))
abline(h=0.5)


#40%
clusterSetRNGStream(cl, 0120)     
ANS.AFT.wei.40.200 = pblapply(X = lapply(1:500,function(x) c(200)),FUN = function(n){
  AUC.Weibull.AFT.np.old = AUC.Weibull.AFT.np.new = c()
  z = rnorm(n)
  u = runif(n,0,1)
  t = (-log(u)/(lambda*exp(alpha*beta*z)))^(1/alpha)
  c = rexp(n,0.07)  #設限率20% <- 0.02  40% <- 0.07
  pt = exp(seq(-5,5,0.25))
  Stime = ifelse(t<c,t,c)
  status = ifelse(t<c,1,0)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  beta.aft <- -aftsrr(Surv(Stime,status)~z)$beta
  fit.aft.np.old = MyrisksetAUC.wai.old(Stime , entry = NULL, status, beta.aft*z, "AFT"
                                        , max(pt) , AFT.method.1 = "1", AFT.method.2 = "1"
                                        , bw.1="1/2",bw.2 = "1/2", bw.mul.1=10, bw.mul.2=0.5, CI=F, confi=0.95
                                        , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  fit.aft.np.new = MyrisksetAUC.wai.new(Stime , entry = NULL, status, beta.aft*z, "AFT"
                                        , max(pt) , AFT.method.1 = "1", AFT.method.2 = "1"
                                        , bw.1="1/2",bw.2 = "1/2", bw.mul.1=10, bw.mul.2=0.5, CI=F, confi=0.95
                                        , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  for( r in 1:41){
    AUC.Weibull.AFT.np.old[r] = fit.aft.np.old$AUC[which.min(abs(fit.aft.np.old$utimes-exp(-5.25+0.25*r)))]
    AUC.Weibull.AFT.np.new[r] = fit.aft.np.new$AUC[which.min(abs(fit.aft.np.new$utimes-exp(-5.25+0.25*r)))]
  }
  Cindex.Weibull.AFT.np.old = fit.aft.np.old$Cindex.old
  Cindex.Weibull.AFT.np.new = fit.aft.np.new$Cindex.new
  
  A = data.frame(EST = c(AUC.Weibull.AFT.np.old,Cindex.Weibull.AFT.np.old),
                 Modify = c(AUC.Weibull.AFT.np.new,Cindex.Weibull.AFT.np.new))
  return(A)
},cl = cl)

AFT.wei.40.200 = matrix(unlist(ANS.AFT.wei.40.200,use.names = F),nrow = 500,byrow = T)
AFT.wei.40.200.old = AFT.wei.40.200[,1:42]
AFT.wei.40.200.new = AFT.wei.40.200[,43:84]

mse.wei.AFT.old.40 = c()
mse.wei.AFT.new.40 = c()
for (i in 1:length(AUC.Weibull.AFT.TRUE)) {
  mse.wei.AFT.old.40[i] = mse(AUC.Weibull.AFT.TRUE[i],AFT.wei.40.200.old[,i])
  mse.wei.AFT.new.40[i] = mse(AUC.Weibull.AFT.TRUE[i],AFT.wei.40.200.new[,i])
}
plot(seq(-5,5,0.25), mse.wei.AFT.old.40, type="l", xlab="Log time", ylab="MSE",lwd = 2,col = "red",main = "Weibull-AFT for 40% censoring")
lines(seq(-5,5,0.25), mse.wei.AFT.new.40, type="l",lwd = 2,col = "blue")
legend("topleft", c("original","Modify"), cex=1, fill=c("red","blue"))

which(AFT.wei.40.200.old < 0.5)
which(AFT.wei.40.200.new < 0.5)
wei.AFT.40.nci.200 = data.frame(Logtime = c(seq(-5,5,0.25),"C"),
                                True = round(c(AUC.Weibull.AFT.TRUE,Cindex.Weibull.AFT.TRUE),3),
                                EST = round(apply(AFT.wei.40.200.old,2,mean),3),
                                SD.old = round(apply(AFT.wei.40.200.old,2,sd),3),
                                MSE.old = round(c(10*mse.wei.AFT.old.40,10*mse(AFT.wei.40.200.old[,42],Cindex.Weibull.AFT.TRUE)),3),
                                Modify = round(apply(AFT.wei.40.200.new,2,mean),3),
                                SD.new = round(apply(AFT.wei.40.200.new,2,sd),3),
                                MSE.new = round(c(10*mse.wei.AFT.new.40,10*mse(AFT.wei.40.200.new[,42],Cindex.Weibull.AFT.TRUE)),3))
wei.AFT.40.nci.200
wei.AFT.40.nci.200[c(17,21,25,29,33,35,37,38,39,42),]

plot(wei.AFT.40.nci.200$Logtime[-c(42)],wei.AFT.40.nci.200$EST[-c(42)],type = "l", xlab="Log time", ylab="AUC",lwd = 2,ylim=c(0.2,1.0),col = "blue",main = "Weibull-AFT for 40% censoring")
lines(wei.AFT.40.nci.200$Logtime[-c(42)],wei.AFT.40.nci.200$Modify[-c(42)], type="l",lwd = 2,col = "red")
lines(wei.AFT.40.nci.200$Logtime[-c(42)],wei.AFT.40.nci.200$True[-c(42)], type="l",lwd = 2,col = "black")
legend("topright", c("original","Modify","Target"), cex=1, fill=c("blue","red","black"))
abline(h=0.5)


stopCluster(cl)
