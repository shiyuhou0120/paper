library(parallel)
library(pbapply)
mu = 1.0 ; sigma = 2.5 ;beta = 1

AUC.Lognormal.AFT.TRUE.value = apply(AFT.ln,2,mean)
Cindex.Lognormal.AFT.TRUE.value = mean(AFT.ln.C)

cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores-2)
clusterExport(cl, c("beta","mu","sigma","MyrisksetAUC.wai.new","MyrisksetAUC.wai.old","EST.cindex.std",
                    "AFTWeights.1","AFTWeights.2","AFTWeights.1.new","AFTWeights.2.new"))
clusterEvalQ(cl, c(library(aftgee),library(survival),library(risksetROC),library(timereg)))

# 無信賴區間估計 ####
clusterSetRNGStream(cl, 0120)     
ANS.AFT.ln.20.200 = pblapply(X = lapply(1:500,function(x) c(200)),FUN = function(n){
  AUC.Lognormal.AFT.np.old = AUC.Lognormal.AFT.np.new = c()
  z = rnorm(n)  #固定
  c = rexp(n,0.02)  # 20% <- 0.02  ,  40% <- 0.13
  u = runif(n,0,1)
  t = exp(sigma*qnorm(1-u)-beta*z+mu)
  pt = exp(seq(-6,6,0.5))
  status = ifelse(t<c,1,0)
  Stime = ifelse(t<c,t,c)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  beta.aft <- -aftsrr(Surv(Stime,status)~z)$beta
  fit.AFT.np.old = MyrisksetAUC.wai.old(Stime , entry = NULL, status, beta.aft*z, "AFT" 
                                        , max(pt) , AFT.method.1 = "1", AFT.method.2 = "1"
                                        , bw.1="1/2",bw.2 = "1/2", bw.mul.1=10, bw.mul.2=0.5, CI=F, confi=0.95
                                        , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  fit.AFT.np.new = MyrisksetAUC.wai.new(Stime , entry = NULL, status, beta.aft*z, "AFT" 
                                        , max(pt) , AFT.method.1 = "1", AFT.method.2 = "1"
                                        , bw.1="1/2",bw.2 = "1/2", bw.mul.1=10, bw.mul.2=0.5, CI=F, confi=0.95
                                        , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  
  for( r in 1:25){
    AUC.Lognormal.AFT.np.old[r] = fit.AFT.np.old$AUC[which.min(abs(fit.AFT.np.old$utimes-exp(-6.5+0.5*r)))]
    AUC.Lognormal.AFT.np.new[r] = fit.AFT.np.new$AUC[which.min(abs(fit.AFT.np.new$utimes-exp(-6.5+0.5*r)))]
  }
  Cindex.Lognormal.AFT.np.old = fit.AFT.np.old$Cindex.old
  Cindex.Lognormal.AFT.np.new = fit.AFT.np.new$Cindex.new

  A = data.frame(EST = c(AUC.Lognormal.AFT.np.old,Cindex.Lognormal.AFT.np.old),
                 Modify = c(AUC.Lognormal.AFT.np.new,Cindex.Lognormal.AFT.np.new))
  return(A)
},cl = cl)

AFT.ln.20.200 = matrix(unlist(ANS.AFT.ln.20.200,use.names = F),nrow = 500,byrow = T)
AFT.ln.20.200.old = AFT.ln.20.200[,1:26]
AFT.ln.20.200.new = AFT.ln.20.200[,27:52]

#計算MSE
mse.ln.AFT.old.20 = c()
mse.ln.AFT.new.20 = c()
for (i in 1:length(AUC.Lognormal.AFT.TRUE.value)) {
  mse.ln.AFT.old.20[i] = mse(AUC.Lognormal.AFT.TRUE.value[i],AFT.ln.20.200.old[,i])
  mse.ln.AFT.new.20[i] = mse(AUC.Lognormal.AFT.TRUE.value[i],AFT.ln.20.200.new[,i])
}
plot(seq(-6,6,0.5), mse.ln.AFT.old.20, type="l", xlab="Log time", ylab="MSE",lwd = 2,col = "red",main = "Lognormal-AFT for 20% censoring")
lines(seq(-6,6,0.5), mse.ln.AFT.new.20, type="l",lwd = 2,col = "blue")
legend("topleft", c("original","Modify"), cex=1, fill=c("red","blue"))

which(AFT.ln.20.200.old < 0.5)
which(AFT.ln.20.200.new < 0.5)
ln.AFT.20.nci.200 = data.frame(Logtime = c(seq(-6,6,0.5),"C"),
                               True = round(c(apply(AFT.ln,2,mean),mean(AFT.ln.C)),3),
                               EST = round(apply(AFT.ln.20.200.old,2,mean),3),
                               SD.old = round(apply(AFT.ln.20.200.old,2,sd),3),
                               MSE.old = round(c(10*mse.ln.AFT.old.20,10*mse(AFT.ln.20.200.old[,26],Cindex.Lognormal.AFT.TRUE.value)),3),
                               Modify = round(apply(AFT.ln.20.200.new,2,mean),3),
                               SD.new = round(apply(AFT.ln.20.200.new,2,sd),3),
                               MSE.new = round(c(10*mse.ln.AFT.new.20,10*mse(AFT.ln.20.200.new[,26],Cindex.Lognormal.AFT.TRUE.value)),3))
ln.AFT.20.nci.200
ln.AFT.20.nci.200[c(seq(3,23,2),26),]

#畫圖
plot(ln.AFT.20.nci.200$Logtime[-c(24:26)],ln.AFT.20.nci.200$EST[-c(24:26)],type = "l", xlab="Log time", ylab="AUC",lwd = 2,ylim=c(0.2,1.0),col = "blue",main = "Lognormal-AFT for 20% censoring")
lines(ln.AFT.20.nci.200$Logtime[-c(24:26)],ln.AFT.20.nci.200$Modify[-c(24:26)], type="l",lwd = 2,col = "red")
lines(ln.AFT.20.nci.200$Logtime[-c(24:26)],ln.AFT.20.nci.200$True[-c(24:26)], type="l",lwd = 2,col = "black")
legend("topright", c("original","Modify","Target"), cex=1, fill=c("blue","red","black"))
abline(h=0.5)


#40%
clusterSetRNGStream(cl, 0120) 
ANS.AFT.ln.40.200 = pblapply(X = lapply(1:500,function(x) c(200)),FUN = function(n){
  AUC.Lognormal.AFT.np.old = AUC.Lognormal.AFT.np.new = c()
  z = rnorm(n)  #固定
  c = rexp(n,0.13)  # 20% <- 0.02  ,  40% <- 0.13
  u = runif(n,0,1)
  t = exp(sigma*qnorm(1-u)-beta*z+mu)
  pt = exp(seq(-6,6,0.5))
  status = ifelse(t<c,1,0)
  Stime = ifelse(t<c,t,c)
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  beta.aft <- -aftsrr(Surv(Stime,status)~z)$beta
  fit.AFT.np.old = MyrisksetAUC.wai.old(Stime , entry = NULL, status, beta.aft*z, "AFT" 
                                        , max(pt) , AFT.method.1 = "1", AFT.method.2 = "1"
                                        , bw.1="1/2",bw.2 = "1/2", bw.mul.1=10, bw.mul.2=0.5, CI=F, confi=0.95
                                        , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  fit.AFT.np.new = MyrisksetAUC.wai.new(Stime , entry = NULL, status, beta.aft*z, "AFT" 
                                        , max(pt) , AFT.method.1 = "1", AFT.method.2 = "1"
                                        , bw.1="1/2",bw.2 = "1/2", bw.mul.1=10, bw.mul.2=0.5, CI=F, confi=0.95
                                        , weight = "rescale", plot=F, type="l",xlab="Time", ylab="AUC")
  for( r in 1:25){
    AUC.Lognormal.AFT.np.old[r] = fit.AFT.np.old$AUC[which.min(abs(fit.AFT.np.old$utimes-exp(-6.5+0.5*r)))]
    AUC.Lognormal.AFT.np.new[r] = fit.AFT.np.new$AUC[which.min(abs(fit.AFT.np.new$utimes-exp(-6.5+0.5*r)))]
  }
  Cindex.Lognormal.AFT.np.old = fit.AFT.np.old$Cindex.old
  Cindex.Lognormal.AFT.np.new = fit.AFT.np.new$Cindex.new
  
  A = data.frame(EST = c(AUC.Lognormal.AFT.np.old,Cindex.Lognormal.AFT.np.old),
                 Modify = c(AUC.Lognormal.AFT.np.new,Cindex.Lognormal.AFT.np.new))
  return(A)
},cl = cl)

AFT.ln.40.200 = matrix(unlist(ANS.AFT.ln.40.200,use.names = F),nrow = 500,byrow = T)
AFT.ln.40.200.old = AFT.ln.40.200[,1:26]
AFT.ln.40.200.new = AFT.ln.40.200[,27:52]

mse.ln.AFT.old.40 = c()
mse.ln.AFT.new.40 = c()
for (i in 1:length(AUC.Lognormal.AFT.TRUE.value)) {
  mse.ln.AFT.old.40[i] = mse(AUC.Lognormal.AFT.TRUE.value[i],AFT.ln.40.200.old[,i])
  mse.ln.AFT.new.40[i] = mse(AUC.Lognormal.AFT.TRUE.value[i],AFT.ln.40.200.new[,i])
}
plot(seq(-6,6,0.5), mse.ln.AFT.old.40, type="l", xlab="Log time", ylab="MSE",lwd = 2,col = "red",main = "Lognormal-AFT for 40% censoring")
lines(seq(-6,6,0.5), mse.ln.AFT.new.40, type="l",lwd = 2,col = "blue")
legend("topleft", c("original","Modify"), cex=1, fill=c("red","blue"))

which(AFT.ln.40.200.old < 0.5)
which(AFT.ln.40.200.new < 0.5)

ln.AFT.40.nci.200 = data.frame(Logtime = c(seq(-6,6,0.5),"C"),
                               True = round(c(apply(AFT.ln,2,mean),mean(AFT.ln.C)),3),
                               EST = round(apply(AFT.ln.40.200.old,2,mean),3),
                               SD.old = round(apply(AFT.ln.40.200.old,2,sd),3),
                               MSE.old = round(c(10*mse.ln.AFT.old.40,10*mse(AFT.ln.40.200.old[,26],Cindex.Lognormal.AFT.TRUE.value)),3),
                               Modify = round(apply(AFT.ln.40.200.new,2,mean),3),
                               SD.new = round(apply(AFT.ln.40.200.new,2,sd),3),
                               MSE.new = round(c(10*mse.ln.AFT.new.40,10*mse(AFT.ln.40.200.new[,26],Cindex.Lognormal.AFT.TRUE.value)),3))
ln.AFT.40.nci.200
ln.AFT.40.nci.200[c(seq(3,23,2),26),]

plot(ln.AFT.40.nci.200$Logtime[-c(24:26)],ln.AFT.40.nci.200$EST[-c(24:26)],type = "l", xlab="Log time", ylab="AUC",lwd = 2,ylim=c(0.2,1.0),col = "blue",main = "Lognormal-AFT for 40% censoring")
lines(ln.AFT.40.nci.200$Logtime[-c(24:26)],ln.AFT.40.nci.200$Modify[-c(24:26)], type="l",lwd = 2,col = "red")
lines(ln.AFT.40.nci.200$Logtime[-c(24:26)],ln.AFT.40.nci.200$True[-c(24:26)], type="l",lwd = 2,col = "black")
legend("topright", c("original","Modify","Target"), cex=1, fill=c("blue","red","black"))
abline(h=0.5)


stopCluster(cl)
