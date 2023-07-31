library(parallel)
library(pbapply)

#計算AUC真實值
#AUC.true wei_PH
alpha = 1.5 ; lambda = 0.25 ; beta = 1
cox.wei = matrix(0,nrow = 1000,ncol = 13)
for (i in 1:1000) {
  cox.wei[i,] = sapply(1:13,function(r) Cox.dist.Weibull.true(exp(-3.5+0.5*r),alpha,lambda,beta)$AUC)
}
round(apply(cox.wei,2,mean),3)
apply(cox.wei,2,sd)

#AUC.true ln_PH
mu = 0.0 ; sigma = 1 ;beta <- 1
cox.ln = matrix(0,nrow = 1000,ncol = 15)
for (i in 1:1000) {
  cox.ln[i,] = sapply(1:15,function(r) Cox.dist.lognormal.true(exp(-3.5+0.5*r),mu,sigma,beta)$AUC)
}
round(apply(cox.ln,2,mean),3)
apply(cox.ln,2,sd)

#AUC.true ln_AFT
mu = 1.0 ; sigma = 2.5 ;beta = 1
AFT.ln = matrix(0,nrow = 100,ncol = 25)
for (i in 1:100) {
  AFT.ln[i,] = sapply(1:25,function(r) AFT.dist.lognormal.true(exp(-6.5+0.5*r) , mu ,sigma, beta)$AUC)
}
round(apply(AFT.ln,2,mean),3)
apply(AFT.ln,2,sd)

#AUC.true wei_AFT
alpha = 1.5 ; lambda = 0.05 ; beta = 1.0
AFT.wei = matrix(0,nrow = 100,ncol = 41)
for (i in 1:100) {
  AFT.wei[i,] = sapply(1:41,function(r) AFT.dist.Weibull.true(exp(-5.25+0.25*r) , alpha , lambda, beta)$AUC)
}
round(apply(AFT.wei,2,mean),3)
apply(AFT.wei,2,sd)

#AUC.true ln_PO
mu = 1.0 ; sigma = 0.3 ; beta=1
PO.ln = matrix(0,nrow = 100,ncol = 17)
for (i in 1:100) {
  PO.ln[i,] = sapply(1:17,function(r) PO.dist.lognormal.true(exp(-2.25+0.25*r),mu,sigma,beta)$AUC )
}
round(apply(PO.ln,2,mean),3)
apply(PO.ln,2,sd)

#AUC.true llog_PO
beta = 1.0 ; mu = 0.50 ; sigma = 0.2 
PO.llog = matrix(0,nrow = 100,ncol = 17)
for (i in 1:100) {
  PO.llog[i,] = sapply(1:17,function(r) PO.dist.loglogistic.true(exp(-2.25+0.25*r) , mu ,sigma, beta)$AUC)
}
round(apply(PO.llog,2,mean),3)
apply(PO.llog,2,sd)


#計算C真實值
#C.true wei_PH
alpha = 1.5 ; lambda = 0.25 ; beta = 1
cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores-1)
clusterExport(cl, c("beta","alpha","lambda","MyrisksetAUC.wai.new","MyrisksetAUC.wai.old","EST.cindex.std",
                    "CoxWeights","AUC.Weibull.PH.true","Cindex.Weibull.PH.true","MyrisksetAUC.Weibull.True"))
clusterEvalQ(cl, c(library(aftgee),library(survival),library(risksetROC),library(timereg)))
clusterSetRNGStream(cl, 0120)
C.PH.wei.true = pblapply(X = lapply(1:100,function(x) c(200)),FUN = function(n){
  C.true = c()
  C.true = MyrisksetAUC.Weibull.True(alpha,lambda,beta,0,1,method = "Cox")$Cindex
  return(C.true)
},cl = cl)

PH.wei.C = matrix(unlist(C.PH.wei.true,use.names = F),nrow = 100,byrow = T)
mean(PH.wei.C)
sd(PH.wei.C)
stopCluster(cl)

#C.true ln_PH
mu = 0.0 ; sigma = 1 ;beta <- 1
cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores-1)
clusterExport(cl, c("beta","mu","sigma","MyrisksetAUC.wai.new","MyrisksetAUC.wai.old","EST.cindex.std",
                    "CoxWeights","AUC.lognormal.PH.TRUE","Cindex.lognormal.PH.TRUE","MyrisksetAUC.lognormal.True"))
clusterEvalQ(cl, c(library(aftgee),library(survival),library(risksetROC),library(timereg)))
clusterSetRNGStream(cl, 0120)
C.ln.wei.true = pblapply(X = lapply(1:100,function(x) c(200)),FUN = function(n){
  C.true = c()
  C.true = MyrisksetAUC.lognormal.True(mu,sigma,beta,mu.z=0,tau=1,method = "Cox")$Cindex
  return(C.true)
},cl = cl)

PH.ln.C = matrix(unlist(C.ln.wei.true,use.names = F),nrow = 100,byrow = T)
round(mean(PH.ln.C),3)
sd(PH.ln.C)
stopCluster(cl)

#C.true ln_AFT
mu = 1.0 ; sigma = 2.5 ;beta = 1
cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores-1)
clusterExport(cl, c("beta","mu","sigma","MyrisksetAUC.wai.new","MyrisksetAUC.wai.old","EST.cindex.std","MyrisksetAUC.lognormal.True",
                    "AFTWeights.1","AFTWeights.2","AUC.Lognormal.AFT.TRUE.value",'Cindex.Lognormal.AFT.TRUE.value'))
clusterEvalQ(cl, c(library(aftgee),library(survival),library(risksetROC),library(timereg)))
clusterSetRNGStream(cl, 0122)
C.ln.AFT.true = pblapply(X = lapply(1:100,function(x) c(200)),FUN = function(n){
  C.true = c()
  C.true = MyrisksetAUC.lognormal.True(mu,sigma,beta,mu.z=0,1,method = "AFT")$Cindex
  return(C.true)
},cl = cl)

AFT.ln.C = matrix(unlist(C.ln.AFT.true,use.names = F),nrow = 100,byrow = T)
mean(AFT.ln.C)
sd(AFT.ln.C)
stopCluster(cl)

#C.true ln_PO
mu = 1.0 ; sigma = 0.3 ; beta=1
cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores-1)
clusterExport(cl, c("beta","mu","sigma","MyrisksetAUC.wai.new","MyrisksetAUC.wai.old","MyrisksetAUC.lognormal.True",
                    "EST.cindex.std","POWeights","AUC.lognormal.PO.TRUE","Cindex.lognormal.PO.TRUE"))
clusterEvalQ(cl, c(library(aftgee),library(survival),library(risksetROC),library(timereg)))
clusterSetRNGStream(cl, 0120)
C.ln.PO.true = pblapply(X = lapply(1:100,function(x) c(200)),FUN = function(n){
  C.true = c()
  C.true = MyrisksetAUC.lognormal.True(mu,sigma,beta,mu.z=0,1,method = "PO")$Cindex
  return(C.true)
},cl = cl)

PO.ln.C = matrix(unlist(C.ln.PO.true,use.names = F),nrow = 100,byrow = T)
mean(PO.ln.C)
sd(PO.ln.C)
stopCluster(cl)

#C.ture wei_AFT
alpha = 1.5 ; lambda = 0.05 ; beta = 1.0
cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores-1)
clusterExport(cl, c("MyrisksetAUC.wai.new","MyrisksetAUC.wai.old","AFTWeights.1","AFTWeights.2","EST.cindex.std",
                    "lambda","alpha","beta","AUC.Weibull.AFT.TRUE.value","MyrisksetAUC.Weibull.True"))#"Cindex.Weibull.AFT.value"))
clusterEvalQ(cl, c(library(muhaz),library(aftgee),library(survival),library(risksetROC)))
clusterSetRNGStream(cl, 0120)
C.wei.AFT.true = pblapply(X = lapply(1:100,function(x) c(200)),FUN = function(n){
  C.true = c()
  C.true = MyrisksetAUC.Weibull.True(alpha,lambda,beta,0,1,method = "AFT")$Cindex
  return(C.true)
},cl = cl)

AFT.wei.C = matrix(unlist(C.wei.AFT.true,use.names = F),nrow = 100,byrow = T)
round(mean(AFT.wei.C),3)
sd(AFT.wei.C)
stopCluster(cl)

#C.true llog_PO
beta = 1.0 ; mu = 0.50 ; sigma = 0.2 
cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores-1)
clusterExport(cl, c("beta","EST.cindex.std","mu","sigma","MyrisksetAUC.wai.old","MyrisksetAUC.wai.new",
                    "AUC.Loglog.PO.TRUE","Cindex.Loglog.PO.TRUE","POWeights","AFTWeights.1","MyrisksetAUC.loglogistic.True"))
clusterEvalQ(cl, c(library(muhaz),library(aftgee),library(survival),library(risksetROC),
                   library(maxLik),library(timereg)))
clusterSetRNGStream(cl, 0120)
C.llog.PO.true = pblapply(X = lapply(1:100,function(x) c(200)),FUN = function(n){
  C.true = c()
  C.true = MyrisksetAUC.loglogistic.True(mu,sigma,beta,mu.z=0,1,method = "PO")$Cindex
  return(C.true)
},cl = cl)

PO.llog.C = matrix(unlist(C.llog.PO.true,use.names = F),nrow = 100,byrow = T)
mean(PO.llog.C)
sd(PO.llog.C)
stopCluster(cl)











