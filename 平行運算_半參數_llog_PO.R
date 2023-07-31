library(parallel)
library(pbapply)
library(timereg)
beta = 1.0 ; mu = 0.50 ; sigma = 0.2 

AUC.Loglog.PO.TRUE = apply(PO.llog,2,mean)
Cindex.Loglog.PO.TRUE = mean(PO.llog.C)

cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores-2)
clusterExport(cl, c("beta","EST.cindex.std","mu","sigma","MyrisksetAUC.wai.old","MyrisksetAUC.wai.new",
                    "POWeights","POWeights.new"))
clusterEvalQ(cl, c(library(muhaz),library(aftgee),library(survival),library(risksetROC),
                   library(maxLik),library(timereg)))
clusterSetRNGStream(cl, 01205)
ANS.PO.llog.20.200 = pblapply(X = lapply(1:500,function(x) c(200)),FUN = function(n){
  AUC.Loglog.PO.old = AUC.Loglog.PO.new = c()
  z <- rnorm(n)
  u <- runif(n)
  t <- exp(sigma*(log((1/u)-1)-beta*z) + mu) # loglogistic
  c <- rexp(n,0.12)  #0.12 <- 20%; 0.3 <- 40%
  pt = exp(seq(-2,2,0.25))
  status = ifelse(t<c,1,0) ; Stime = ifelse(t<c,t,c)
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
    AUC.Loglog.PO.old[r] = fit.po.old$AUC[which.min(abs(fit.po.old$utimes-exp(-2.25+0.25*r)))]
    AUC.Loglog.PO.new[r] = fit.po.new$AUC[which.min(abs(fit.po.new$utimes-exp(-2.25+0.25*r)))]
  }
  Cindex.Loglog.PO.np.old = fit.po.old$Cindex.old
  Cindex.Loglog.PO.np.new = fit.po.new$Cindex.new

  A = data.frame(EST = c(AUC.Loglog.PO.old,Cindex.Loglog.PO.np.old),
                 Modify = c(AUC.Loglog.PO.new,Cindex.Loglog.PO.np.new))
  return(A)
},cl = cl)

PO.llog.20.200 = matrix(unlist(ANS.PO.llog.20.200,use.names = F),nrow = 500,byrow = T)
PO.llog.20.200.old = PO.llog.20.200[,1:18]
PO.llog.20.200.new = PO.llog.20.200[,19:36]

#計算MSE
mse.llog.PO.old.20 = c()
mse.llog.PO.new.20 = c()
for (i in 1:length(AUC.Loglog.PO.TRUE)) {
  mse.llog.PO.old.20[i] = mse(AUC.Loglog.PO.TRUE[i],PO.llog.20.200.old[,i])
  mse.llog.PO.new.20[i] = mse(AUC.Loglog.PO.TRUE[i],PO.llog.20.200.new[,i])
}
plot(seq(-2,2,0.25), mse.llog.PO.old.20, type="l", xlab="Log time", ylab="MSE",lwd = 2,col = "red",main = "Loglogistic-PO for 20% censoring")
lines(seq(-2,2,0.25), mse.llog.PO.new.20, type="l",lwd = 2,col = "blue")
legend("topleft", c("original","Modify"), cex=1, fill=c("red","blue"))

which(PO.llog.20.200.old < 0.5)
which(PO.llog.20.200.new < 0.5)
llog.PO.20.nci.200 = data.frame(Logtime = c(seq(-2,2,0.25),"C"),
                                True = round(c(AUC.Loglog.PO.TRUE,Cindex.Loglog.PO.TRUE),3),
                                EST = round(apply(PO.llog.20.200.old,2,mean),3),
                                SD.old = round(apply(PO.llog.20.200.old,2,sd),3),
                                MSE.old = round(c(10*mse.llog.PO.old.20,10*mse(PO.llog.20.200.old[,18],Cindex.Loglog.PO.TRUE)),3),
                                Modify = round(apply(PO.llog.20.200.new,2,mean),3),
                                SD.new = round(apply(PO.llog.20.200.new,2,sd),3),
                                MSE.new = round(c(10*mse.llog.PO.new.20,10*mse(PO.llog.20.200.new[,18],Cindex.Loglog.PO.TRUE)),3))
llog.PO.20.nci.200
llog.PO.20.nci.200[c(5,7,9:18),]

#畫圖
plot(llog.PO.20.nci.200$Logtime[-c(18)],llog.PO.20.nci.200$EST[-c(18)],type = "l", xlab="Log time", ylab="AUC",lwd = 2,ylim=c(0.0,1.0),col = "blue",main = "Loglogistic-PO for 20% censoring")
lines(llog.PO.20.nci.200$Logtime[-c(18)],llog.PO.20.nci.200$Modify[-c(18)], type="l",lwd = 2,col = "red")
lines(llog.PO.20.nci.200$Logtime[-c(18)],llog.PO.20.nci.200$True[-c(18)], type="l",lwd = 2,col = "black")
legend("topright", c("original","Modify","Target"), cex=1, fill=c("blue","red","black"))
abline(h=0.5)

#40%
clusterSetRNGStream(cl, 01205)
ANS.PO.llog.40.200 = pblapply(X = lapply(1:500,function(x) c(200)),FUN = function(n){
  AUC.Loglog.PO.old = AUC.Loglog.PO.new = c()
  z <- rnorm(n)
  u <- runif(n)
  t <- exp(sigma*(log((1/u)-1)-beta*z) + mu) # loglogistic
  c <- rexp(n,0.3)  #0.12 <- 20%; 0.3 <- 40%
  pt = exp(seq(-2,2,0.25))
  status = ifelse(t<c,1,0) ; Stime = ifelse(t<c,t,c)
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
    AUC.Loglog.PO.old[r] = fit.po.old$AUC[which.min(abs(fit.po.old$utimes-exp(-2.25+0.25*r)))]
    AUC.Loglog.PO.new[r] = fit.po.new$AUC[which.min(abs(fit.po.new$utimes-exp(-2.25+0.25*r)))]
  }
  Cindex.Loglog.PO.np.old = fit.po.old$Cindex.old
  Cindex.Loglog.PO.np.new = fit.po.new$Cindex.new
  
  A = data.frame(EST = c(AUC.Loglog.PO.old,Cindex.Loglog.PO.np.old),
                 Modify = c(AUC.Loglog.PO.new,Cindex.Loglog.PO.np.new))
  return(A)
},cl = cl)

PO.llog.40.200 = matrix(unlist(ANS.PO.llog.40.200,use.names = F),nrow = 500,byrow = T)
PO.llog.40.200.old = PO.llog.40.200[,1:18]
PO.llog.40.200.new = PO.llog.40.200[,19:36]

mse.llog.PO.old.40 = c()
mse.llog.PO.new.40 = c()
for (i in 1:length(AUC.Loglog.PO.TRUE)) {
  mse.llog.PO.old.40[i] = mse(AUC.Loglog.PO.TRUE[i],PO.llog.40.200.old[,i])
  mse.llog.PO.new.40[i] = mse(AUC.Loglog.PO.TRUE[i],PO.llog.40.200.new[,i])
}
plot(seq(-2,2,0.25), mse.llog.PO.old.40, type="l", xlab="Log time", ylab="MSE",lwd = 2,col = "red",main = "Loglogistic-PO for 40% censoring")
lines(seq(-2,2,0.25), mse.llog.PO.new.40, type="l",lwd = 2,col = "blue")
legend("topleft", c("original","Modify"), cex=1, fill=c("red","blue"))

which(PO.llog.40.200.old < 0.5)
which(PO.llog.40.200.new < 0.5)
llog.PO.40.nci.200 = data.frame(Logtime = c(seq(-2,2,0.25),"C"),
                                True = round(c(AUC.Loglog.PO.TRUE,Cindex.Loglog.PO.TRUE),3),
                                EST = round(apply(PO.llog.40.200.old,2,mean),3),
                                SD.old = round(apply(PO.llog.40.200.old,2,sd),3),
                                MSE.old = round(c(10*mse.llog.PO.old.40,10*mse(PO.llog.40.200.old[,18],Cindex.Loglog.PO.TRUE)),3),
                                Modify = round(apply(PO.llog.40.200.new,2,mean),3),
                                SD.new = round(apply(PO.llog.40.200.new,2,sd),3),
                                MSE.new = round(c(10*mse.llog.PO.new.40,10*mse(PO.llog.40.200.new[,18],Cindex.Loglog.PO.TRUE)),3))
llog.PO.40.nci.200
llog.PO.40.nci.200[c(5,7,9:18),]

plot(llog.PO.40.nci.200$Logtime[-c(18)],llog.PO.40.nci.200$EST[-c(18)],type = "l", xlab="Log time", ylab="AUC",lwd = 2,ylim=c(0.0,1.0),col = "blue",main = "Loglogistic-PO for 40% censoring")
lines(llog.PO.40.nci.200$Logtime[-c(18)],llog.PO.40.nci.200$Modify[-c(18)], type="l",lwd = 2,col = "red")
lines(llog.PO.40.nci.200$Logtime[-c(18)],llog.PO.40.nci.200$True[-c(18)], type="l",lwd = 2,col = "black")
legend("topright", c("original","Modify","Target"), cex=1, fill=c("blue","red","black"))
abline(h=0.5)


stopCluster(cl)
