library(parallel)
library(pbapply)
mu = 0.0 ; sigma = 1 ;beta <- 1

AUC.lognormal.PH.TRUE = apply(cox.ln,2,mean)
Cindex.lognormal.PH.TRUE = mean(PH.ln.C)

cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores-1)
clusterExport(cl, c("beta","mu","sigma","MyrisksetAUC.wai.new","MyrisksetAUC.wai.old","EST.cindex.std",
                    "CoxWeights","CoxWeights.new"))
clusterEvalQ(cl, c(library(aftgee),library(survival),library(risksetROC),library(timereg)))

# 有信賴區間估計 ####
clusterSetRNGStream(cl, 0120)     
ANS.PH.ln.20.200 = pblapply(X = lapply(1:500,function(x) c(200)),FUN = function(n){
  AUC.Lognormal.PH.np.new = AUC.Lognormal.PH.np.old = c()
  z = rnorm(n)
  u = runif(n)
  c = rexp(n,0.13) #0.13 <- 20% 0.3 <- 40%
  t = exp(sigma*qnorm(1-u^(exp(-beta*z)))+ mu )
  pt = exp(seq(-3,4,0.5))
  status = ifelse(t<c,1,0)
  Stime = ifelse(t<c,t,c)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  beta.ph <- coxph(Surv(Stime,status)~z)$coefficients
  fit.PH.new = MyrisksetAUC.wai.new(Stime , entry = NULL, status, beta.ph*z, "Cox"
                                    , tmax = max(pt) , AFT.method.1 = "1",AFT.method.2 = "1"
                                    , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                    , weight = "rescale", plot=TRUE, type="l",xlab="Time", ylab="AUC")
  
  fit.PH.old = MyrisksetAUC.wai.old(Stime , entry = NULL, status, beta.ph*z, "Cox"
                                    , tmax = max(pt) , AFT.method.1 = "1",AFT.method.2 = "1"
                                    , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                    , weight = "rescale", plot=TRUE, type="l",xlab="Time", ylab="AUC")
  for( r in 1:15){
    AUC.Lognormal.PH.np.old[r] = fit.PH.old$AUC[which.min(abs(fit.PH.old$utimes-exp(-3.5+0.5*r)))]
    AUC.Lognormal.PH.np.new[r] = fit.PH.new$AUC[which.min(abs(fit.PH.new$utimes-exp(-3.5+0.5*r)))]
  }
    Cindex.Lognormal.PH.np.old = fit.PH.old$Cindex.old
    Cindex.Lognormal.PH.np.new = fit.PH.new$Cindex.new
  A = data.frame(EST = c(AUC.Lognormal.PH.np.old,Cindex.Lognormal.PH.np.old),
                 Modify = c(AUC.Lognormal.PH.np.new,Cindex.Lognormal.PH.np.new))
  return(A)
},cl = cl)

PH.ln.20.200 = matrix(unlist(ANS.PH.ln.20.200,use.names = F),nrow = 500,byrow = T)
PH.ln.20.200.old = PH.ln.20.200[,1:16]
PH.ln.20.200.new = PH.ln.20.200[,17:32]
mse.ln.PH.old.20 = c()
mse.ln.PH.new.20 = c()
for (i in 1:length(AUC.lognormal.PH.TRUE)) {
  mse.ln.PH.old.20[i] = mse(AUC.lognormal.PH.TRUE[i],PH.ln.20.200.old[,i])
  mse.ln.PH.new.20[i] = mse(AUC.lognormal.PH.TRUE[i],PH.ln.20.200.new[,i])
}
plot(seq(-3,4,0.5), mse.ln.PH.old.20, type="l", xlab="Log time", ylab="MSE",lwd = 2,col = "red",main = "Lognormal-Cox for 20% censoring")
lines(seq(-3,4,0.5), mse.ln.PH.new.20, type="l",lwd = 2,col = "blue")
legend("topleft", c("original","Modify"), cex=1, fill=c("red","blue"))

which(PH.ln.20.200.old < 0.5)
which(PH.ln.20.200.new < 0.5)
ln.PH.20.nci.200 = data.frame(Logtime = c(seq(-3,4,0.5),"C"),
                              True = c(round(AUC.lognormal.PH.TRUE,3),round(Cindex.lognormal.PH.TRUE,3)),
                              EST = round(apply(PH.ln.20.200.old,2,mean),3),
                              SD.old = round(apply(PH.ln.20.200.old,2,sd),3),
                              MSE.old = round(c(10*mse.ln.PH.old.20,10*mse(PH.ln.20.200.old[,16],Cindex.lognormal.PH.TRUE)),3),
                              Modify = round(apply(PH.ln.20.200.new,2,mean),3),
                              SD.new = round(apply(PH.ln.20.200.new,2,sd),3),
                              MSE.new = round(c(10*mse.ln.PH.new.20,10*mse(PH.ln.20.200.new[,16],Cindex.lognormal.PH.TRUE)),3))
ln.PH.20.nci.200
ln.PH.20.nci.200[c(2:13,16),]
plot(ln.PH.20.nci.200$Logtime[-16],ln.PH.20.nci.200$EST[-16],type = "l", xlab="Log time", ylab="AUC",lwd = 2,ylim=c(0.4,1.0),col = "blue",main = "Lognormal-Cox for 20% censoring")
lines(ln.PH.20.nci.200$Logtime[-16], ln.PH.20.nci.200$Modify[-16], type="l",lwd = 2,col = "red")
lines(ln.PH.20.nci.200$Logtime[-16],ln.PH.20.nci.200$True[-16], type="l",lwd = 2,col = "black")
legend("topright", c("original","Modify","Target"), cex=1, fill=c("blue","red","black"))
abline(h=0.5)

#40%
clusterSetRNGStream(cl, 0120)     
ANS.PH.ln.40.200 = pblapply(X = lapply(1:500,function(x) c(200)),FUN = function(n){
  AUC.Lognormal.PH.np.new = AUC.Lognormal.PH.np.old = c()
  z = rnorm(n)
  u = runif(n)
  c = rexp(n,0.3) #0.13 <- 20% 0.3 <- 40%
  t = exp(sigma*qnorm(1-u^(exp(-beta*z)))+ mu )
  pt = exp(seq(-4,3,0.5))
  status = ifelse(t<c,1,0)
  Stime = ifelse(t<c,t,c)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  beta.ph <- coxph(Surv(Stime,status)~z)$coefficients
  fit.PH.old = MyrisksetAUC.wai.old( Stime , entry = NULL, status, beta.ph*z, "Cox"
                                     , tmax = max(pt) , AFT.method.1 = "1",AFT.method.2 = "1"
                                     , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                     , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  
  fit.PH.new = MyrisksetAUC.wai.new( Stime , entry = NULL, status, beta.ph*z, "Cox"
                                     , tmax = max(pt) , AFT.method.1 = "1",AFT.method.2 = "1"
                                     , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                     , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  for( r in 1:15){
    AUC.Lognormal.PH.np.old[r] = fit.PH.old$AUC[which.min(abs(fit.PH.old$utimes-exp(-3.5+0.5*r)))]
    AUC.Lognormal.PH.np.new[r] = fit.PH.new$AUC[which.min(abs(fit.PH.new$utimes-exp(-3.5+0.5*r)))]
  }
  Cindex.Lognormal.PH.np.old = fit.PH.old$Cindex.old
  Cindex.Lognormal.PH.np.new = fit.PH.new$Cindex.new
  A = data.frame(EST = c(AUC.Lognormal.PH.np.old,Cindex.Lognormal.PH.np.old),
                 Modify = c(AUC.Lognormal.PH.np.new,Cindex.Lognormal.PH.np.new))
  return(A)
},cl = cl)

PH.ln.40.200 = matrix(unlist(ANS.PH.ln.40.200,use.names = F),nrow = 500,byrow = T)
PH.ln.40.200.old = PH.ln.40.200[,1:16]
PH.ln.40.200.new = PH.ln.40.200[,17:32]
mse.ln.PH.old.40 = c()
mse.ln.PH.new.40 = c()
for (i in 1:length(AUC.lognormal.PH.TRUE)) {
  mse.ln.PH.old.40[i] = mse(AUC.lognormal.PH.TRUE[i],PH.ln.40.200.old[,i])
  mse.ln.PH.new.40[i] = mse(AUC.lognormal.PH.TRUE[i],PH.ln.40.200.new[,i])
}
plot(seq(-3,4,0.5), mse.ln.PH.old.40, type="l", xlab="Log time", ylab="MSE",lwd = 2,col = "red",main = "Lognormal-Cox for 40% censoring")
lines(seq(-3,4,0.5), mse.ln.PH.new.40, type="l",lwd = 2,col = "blue")
legend("topleft", c("original","Modify"), cex=1, fill=c("red","blue"))
which(PH.ln.40.200.old < 0.5)
which(PH.ln.40.200.new < 0.5)

ln.PH.40.nci.200 = data.frame(Logtime = c(seq(-3,4,0.5),"C"),
                              True = c(round(AUC.lognormal.PH.TRUE,3),round(Cindex.lognormal.PH.TRUE,3)),
                              EST = round(apply(PH.ln.40.200.old,2,mean),3),
                              SD.old = round(apply(PH.ln.40.200.old,2,sd),3),
                              MSE.old = round(c(10*mse.ln.PH.old.40,10*mse(PH.ln.40.200.old[,16],Cindex.lognormal.PH.TRUE)),3),
                              Modify = round(apply(PH.ln.40.200.new,2,mean),3),
                              SD.new = round(apply(PH.ln.40.200.new,2,sd),3),
                              MSE.new = round(c(10*mse.ln.PH.new.40,10*mse(PH.ln.40.200.new[,16],Cindex.lognormal.PH.TRUE)),3))
ln.PH.40.nci.200
ln.PH.40.nci.200[c(2:13,16),]
plot(ln.PH.40.nci.200$Logtime[-16],ln.PH.40.nci.200$EST[-16],type = "l", xlab="Log time", ylab="AUC",lwd = 2,ylim=c(0.4,1.0),col = "blue",main = "Lognormal-Cox for 40% censoring")
lines(ln.PH.40.nci.200$Logtime[-16], ln.PH.40.nci.200$Modify[-16], type="l",lwd = 2,col = "red")
lines(ln.PH.40.nci.200$Logtime[-16],ln.PH.40.nci.200$True[-16], type="l",lwd = 2,col = "black")
legend("topright", c("original","Modify","Target"), cex=1, fill=c("blue","red","black"))
abline(h=0.5)

stopCluster(cl)
