library(parallel)
library(pbapply)
mu = 1.0 ; sigma = 0.3 ; beta=1

AUC.lognormal.PO.TRUE = apply(PO.ln,2,mean)
Cindex.lognormal.PO.TRUE = round(mean(PO.ln.C),3)

cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores-2)
clusterExport(cl, c("beta","mu","sigma","MyrisksetAUC.wai.new","MyrisksetAUC.wai.old",
                    "EST.cindex.std","POWeights","POWeights.new"))
clusterEvalQ(cl, c(library(aftgee),library(survival),library(risksetROC),library(timereg)))
# 無信賴區間估計 ####
clusterSetRNGStream(cl, 0120)     
ANS.PO.ln.20.200 = pblapply(X = lapply(1:500,function(x) c(200)),FUN = function(n){
  AUC.lognormal.PO.np.old = AUC.lognormal.PO.np.new = c()
  u <- runif(n)
  z <- rnorm(n,0,1) # 不能改 配合參數模型真實TP
  t <- exp(mu+sigma*qnorm((1-(1-(1/u))*exp(beta*z))^(-1))) # lognormal
  c <- rexp(n,0.06)  #0.06 <- 20% 0.16 <- 40%
  pt = exp(seq(-2,2,0.25))
  status <- ifelse(t<c,1,0)
  Stime = ifelse(t<c,t,c)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  data = data.frame(time=rep(0,length(Stime)),time2=Stime,status = status,eta = z)
  fit = prop.odds(Event(time2,status)~eta,max.time=max(Stime),data = data,n.sim=100)
  gamma.t = as.numeric(fit$gamma)
  fit.po.old = MyrisksetAUC.wai.old( Stime , entry = NULL, status, gamma.t*z, "PO"
                                     , max(pt) , AFT.method.1 = "1",AFT.method.2 = "1"
                                     , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                     , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  fit.po.new = MyrisksetAUC.wai.new( Stime , entry = NULL, status, gamma.t*z, "PO"
                                     , max(pt) , AFT.method.1 = "1",AFT.method.2 = "1"
                                     , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                     , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  
  for( r in 1:17){
    AUC.lognormal.PO.np.old[r] = fit.po.old$AUC[which.min(abs(fit.po.old$utimes-exp(-2.25+0.25*r)))]
    AUC.lognormal.PO.np.new[r] = fit.po.new$AUC[which.min(abs(fit.po.new$utimes-exp(-2.25+0.25*r)))]
  }
  Cindex.lognormal.PO.np.old = fit.po.old$Cindex.old
  Cindex.lognormal.PO.np.new = fit.po.new$Cindex.new
  
  A = data.frame(EST = c(AUC.lognormal.PO.np.old,Cindex.lognormal.PO.np.old),
                 Modify = c(AUC.lognormal.PO.np.new,Cindex.lognormal.PO.np.new))
  return(A)
},cl = cl)

PO.ln.20.200 = matrix(unlist(ANS.PO.ln.20.200,use.names = F),nrow = 500,byrow = T)
PO.ln.20.200.old = PO.ln.20.200[,1:18]
PO.ln.20.200.new = PO.ln.20.200[,19:36]

#計算MSE
mse.ln.PO.old.20 = c()
mse.ln.PO.new.20 = c()
for (i in 1:length(AUC.lognormal.PO.TRUE)) {
  mse.ln.PO.old.20[i] = mse(AUC.lognormal.PO.TRUE[i],PO.ln.20.200.old[,i])
  mse.ln.PO.new.20[i] = mse(AUC.lognormal.PO.TRUE[i],PO.ln.20.200.new[,i])
}
plot(seq(-2,2,0.25), mse.ln.PO.old.20, type="l", xlab="Log time", ylab="MSE",lwd = 2,col = "red",main = "Lognormal-PO for 20% censoring")
lines(seq(-2,2,0.25), mse.ln.PO.new.20, type="l",lwd = 2,col = "blue")
legend("topleft", c("original","Modify"), cex=1, fill=c("red","blue"))

which(PO.ln.20.200.old < 0.5)
which(PO.ln.20.200.new < 0.5)
ln.PO.20.nci.200 = data.frame(Logtime = c(seq(-2,2,0.25),"C"),
                              True = c(round(AUC.lognormal.PO.TRUE,3),Cindex.lognormal.PO.TRUE),
                              EST = round(apply(PO.ln.20.200.old,2,mean),3),
                              SD.old = round(apply(PO.ln.20.200.old,2,sd),3),
                              MSE.old = round(c(10*mse.ln.PO.old.20,10*mse(PO.ln.20.200.old[,18],Cindex.lognormal.PO.TRUE)),3),
                              Modify = round(apply(PO.ln.20.200.new,2,mean),3),
                              SD.new = round(apply(PO.ln.20.200.new,2,sd),3),
                              MSE.new = round(c(10*mse.ln.PO.new.20,10*mse(PO.ln.20.200.new[,18],Cindex.lognormal.PO.TRUE)),3))
ln.PO.20.nci.200
ln.PO.20.nci.200[9:18,]

#畫圖
plot(ln.PO.20.nci.200$Logtime[-c(18)],ln.PO.20.nci.200$EST[-c(18)],type = "l", xlab="Log time", ylab="AUC",lwd = 2,ylim=c(0.0,1.0),col = "blue",main = "Lognormal-PO for 20% censoring")
lines(ln.PO.20.nci.200$Logtime[-c(18)],ln.PO.20.nci.200$Modify[-c(18)], type="l",lwd = 2,col = "red")
lines(ln.PO.20.nci.200$Logtime[-c(18)],ln.PO.20.nci.200$True[-c(18)], type="l",lwd = 2,col = "black")
legend("topright", c("original","Modify","Target"), cex=1, fill=c("blue","red","black"))
abline(h=0.5)

#40%
clusterSetRNGStream(cl, 0120)
ANS.PO.ln.40.200 = pblapply(X = lapply(1:500,function(x) c(200)),FUN = function(n){
  AUC.lognormal.PO.np.old = AUC.lognormal.PO.np.new = c()
  u <- runif(n)
  z <- rnorm(n,0,1) # 不能改 配合參數模型真實TP
  t <- exp(mu+sigma*qnorm((1-(1-(1/u))*exp(beta*z))^(-1))) # lognormal
  c <- rexp(n,0.16)  #0.06 <- 20% 0.16 <- 40%
  pt = exp(seq(-2,2,0.25))
  status <- ifelse(t<c,1,0)
  Stime = ifelse(t<c,t,c)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  data = data.frame(time=rep(0,length(Stime)),time2=Stime,status = status,eta = z)
  fit = prop.odds(Event(time2,status)~eta,max.time=max(Stime),data = data,n.sim=100)
  gamma.t = as.numeric(fit$gamma)
  fit.po.old = MyrisksetAUC.wai.old( Stime , entry = NULL, status, gamma.t*z, "PO"
                                     , max(pt) , AFT.method.1 = "1",AFT.method.2 = "1"
                                     , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                     , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  fit.po.new = MyrisksetAUC.wai.new( Stime , entry = NULL, status, gamma.t*z, "PO"
                                     , max(pt) , AFT.method.1 = "1",AFT.method.2 = "1"
                                     , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                     , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  
  for( r in 1:17){
    AUC.lognormal.PO.np.old[r] = fit.po.old$AUC[which.min(abs(fit.po.old$utimes-exp(-2.25+0.25*r)))]
    AUC.lognormal.PO.np.new[r] = fit.po.new$AUC[which.min(abs(fit.po.new$utimes-exp(-2.25+0.25*r)))]
  }
  Cindex.lognormal.PO.np.old = fit.po.old$Cindex.old
  Cindex.lognormal.PO.np.new = fit.po.new$Cindex.new
  
  A = data.frame(EST = c(AUC.lognormal.PO.np.old,Cindex.lognormal.PO.np.old),
                 Modify = c(AUC.lognormal.PO.np.new,Cindex.lognormal.PO.np.new))
  return(A)
},cl = cl)

PO.ln.40.200 = matrix(unlist(ANS.PO.ln.40.200,use.names = F),nrow = 500,byrow = T)
PO.ln.40.200.old = PO.ln.40.200[,1:18]
PO.ln.40.200.new = PO.ln.40.200[,19:36]

mse.ln.PO.old.40 = c()
mse.ln.PO.new.40 = c()
for (i in 1:length(AUC.lognormal.PO.TRUE)) {
  mse.ln.PO.old.40[i] = mse(AUC.lognormal.PO.TRUE[i],PO.ln.40.200.old[,i])
  mse.ln.PO.new.40[i] = mse(AUC.lognormal.PO.TRUE[i],PO.ln.40.200.new[,i])
}
plot(seq(-2,2,0.25), mse.ln.PO.old.40, type="l", xlab="Log time", ylab="MSE",lwd = 2,col = "red",main = "Lognormal-PO for 40% censoring")
lines(seq(-2,2,0.25), mse.ln.PO.new.40, type="l",lwd = 2,col = "blue")
legend("topleft", c("original","Modify"), cex=1, fill=c("red","blue"))

which(PO.ln.40.200.old < 0.5)
which(PO.ln.40.200.new < 0.5)
ln.PO.40.nci.200 = data.frame(Logtime = c(seq(-2,2,0.25),"C"),
                              True = c(round(AUC.lognormal.PO.TRUE,3),Cindex.lognormal.PO.TRUE),
                              EST = round(apply(PO.ln.40.200.old,2,mean),3),
                              SD.old = round(apply(PO.ln.40.200.old,2,sd),3),
                              MSE.old = round(c(10*mse.ln.PO.old.40,10*mse(PO.ln.40.200.old[,18],Cindex.lognormal.PO.TRUE)),3),
                              Modify = round(apply(PO.ln.40.200.new,2,mean),3),
                              SD.new = round(apply(PO.ln.40.200.new,2,sd),3),
                              MSE.new = round(c(10*mse.ln.PO.new.40,10*mse(PO.ln.40.200.new[,18],Cindex.lognormal.PO.TRUE)),3))
ln.PO.40.nci.200
ln.PO.40.nci.200[9:18,]

plot(ln.PO.40.nci.200$Logtime[-c(18)],ln.PO.40.nci.200$EST[-c(18)],type = "l", xlab="Log time", ylab="AUC",lwd = 2,ylim=c(0.0,1.0),col = "blue",main = "Lognormal-PO for 40% censoring")
lines(ln.PO.40.nci.200$Logtime[-c(18)],ln.PO.40.nci.200$Modify[-c(18)], type="l",lwd = 2,col = "red")
lines(ln.PO.40.nci.200$Logtime[-c(18)],ln.PO.40.nci.200$True[-c(18)], type="l",lwd = 2,col = "black")
legend("topright", c("original","Modify","Target"), cex=1, fill=c("blue","red","black"))
abline(h=0.5)


stopCluster(cl)
