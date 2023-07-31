library(survival)
library(risksetROC)
library(muhaz)
library(survPresmooth)
library(bshazard)
library(MASS)
library(mvtnorm)
library(aftgee)
library(maxLik)
library(timereg)

risksetROC.new = function (Stime, entry = NULL, status, marker, predict.time, 
          method = "Cox", span = NULL, order = 1, window = "asymmetric", 
          prop = 0.5, plot = TRUE, type = "l", xlab = "FP", ylab = "TP", ...) 
{
  mChoice <- match(method, c("Cox", "LocalCox", 
                             "Schoenfeld"))
  if (is.na(mChoice)) {
    cat("error in method choice")
    stop(0)
  }
  if (is.null(span) & ((mChoice == 2) || (mChoice == 3))) {
    cat("Need span for methods = \"LocalCox\" or \"Schoenfeld\" \n")
    stop(0)
  }
  p = 1
  eta = marker
  if (length(entry) == 0) {
    entry = rep(0, NROW(Stime))
  }
  time = entry
  time2 = Stime
  if (mChoice == 1) {
    fit = coxph(Surv(time, time2, status) ~ eta)
    new.eta = eta * fit$coefficients
  }
  if (mChoice == 2) {
    nt = length(unique(time2[status == 1]))
    grid.t = time2
    nd = prop
    delta = abs(grid.t - predict.time)/(max(time2[status == 1]) - min(time2[status == 1]))
    keep = delta <= nd
    bfnx.ll = llCoxReg(entry = time[keep == 1], Stime = time2[keep == 1], status = status[keep == 1], 
                       marker = eta[keep == 1], span = span, p = p, window = window)
    gamma.t = bfnx.ll$beta[NROW(bfnx.ll$time[bfnx.ll$time <= predict.time]), ]
    new.eta = eta * gamma.t[1]
  }
  if (mChoice == 3) {
    fit = coxph(Surv(time, time2, status) ~ eta)
    bfnx.SS = SchoenSmooth(fit = fit, Stime = time2, status = status, span = span, order = order)
    gamma.t = bfnx.SS$beta[NROW(bfnx.SS$time[bfnx.SS$time <= predict.time])]
    new.eta = eta * gamma.t
  }
  out = CoxWeights.new(marker = new.eta, entry = time, status = status, 
                       predict.time = predict.time, Stime = time2)
  if (plot == TRUE) {
    plot(out$FP, out$TP, type = type, xlab = xlab, ylab = ylab, 
         ...)
    abline(c(0, 0), c(1, 1))
  }
  return(out)
}

#### PH 模型估計 AUC (舊版)
CoxWeights <- function (marker, Stime, status, predict.time,entry = NULL) {
  eta <- marker
  Target <- predict.time  # 欲估計的時間點
  if (length(entry) == 0) {
    entry = rep(0, NROW(Stime))  
    # 若沒給定起始時間， 則起始設為0
  }
  at.risk <- ((Stime >= Target) & (entry <= Target))  # 風險指標(仍承受風險)
  the.eta <- eta[at.risk] #仍承受風險之觀測者marker
  n <- length(the.eta)
  the.dead <- (Stime == Target) & (status == 1)
  the.dead <- the.dead[at.risk]
  n.dead <- sum(the.dead)
  p0 <- rep(1/(n - n.dead), n)
  if (n.dead > 0) 
    p0[the.dead] <- 0
  p1 <- exp(the.eta)/sum(exp(the.eta))
  ooo <- order(the.eta)
  eta.out <- the.eta[ooo]
  TP <- c(1, 1 - cumsum(p1[ooo]), 0)  # True Positive
  FP <- c(1, 1 - cumsum(p0[ooo]), 0)  # False Positive
  dFP <- abs(FP[-1] - FP[-length(FP)])
  aTP <- 0.5 * (TP[-1] + TP[-length(TP)])
  area <- sum( dFP * aTP )  # 使用梯形法求得 AUC
  out <- list( marker = eta.out , TP = TP, FP = FP, AUC = area )
  return(out)
}
#### PH 模型估計 AUC (修正後)
CoxWeights.new <- function (marker, Stime, status, predict.time,entry = NULL) {
  eta <- marker
  Target <- predict.time  # 欲估計的時間點
  if (length(entry) == 0) {
    entry = rep(0, NROW(Stime))  
    # 若沒給定起始時間， 則起始設為0
  }
  at.risk <- ((Stime >= Target) & (entry <= Target))  # 風險指標(仍承受風險)
  the.eta <- eta[at.risk] #仍承受風險之觀測者marker
  n <- length(the.eta)
  the.dead <- (Stime == Target) & (status == 1)
  the.dead <- the.dead[at.risk]
  n.dead <- sum(the.dead)
  p0 <- rep(1/(n - n.dead), n)
  if (n.dead > 0) 
    p0[the.dead] <- 0
  p1 <- exp(the.eta)/sum(exp(the.eta))
  ooo <- order(the.eta)
  eta.out <- the.eta[ooo]
  TP <- c(1, 1 - cumsum(p1[ooo]), 0)  # True Positive
  FP <- c(1, 1 - cumsum(p0[ooo]), 0)  # False Positive
  a = which(FP[-1] == FP[-length(FP)])
  a = a + 1    #將重複點拿掉
  TP = TP[-a]
  FP = FP[-a]
  #plot(FP,TP)
  dFP <- abs(FP[-1] - FP[-length(FP)])
  aTP <- 0.5 * (TP[-1] + TP[-length(TP)])
  area <- sum( dFP * aTP )  # 使用梯形法求得 AUC
  out <- list( marker = eta.out , TP = TP, FP = FP, AUC = area )
  return(out)
}
#### AFT 模型估計 AUC
##  修正 Kernel smooth 方法估計 AUC (舊版)
AFTWeights.1 <- function(marker , Stime , status , predict.time , entry=NULL,bw = "opt"){
  Target <- predict.time
  eta <- marker
  if(length(entry) == 0){
    entry = rep(0,NROW(Stime))
  }
  u <- exp( eta)*Stime
  find <- exp( eta)*predict.time
  # 修正後的風險函數估計
  haz.fun <- function(u){
    r = log(Stime)+eta ; n <- length(eta)
    if(bw == "opt"){
      bw1 = ((8*sqrt(2)/3)^(0.2)) * sd((log(Stime)+eta)[which(status==1)]) * (n^(-0.2))
      bw2 = (4^(1/3)) * sd(log(Stime)+eta) * (n^(-1/3))
    }
    else if(bw == "1/2"){
      bw1 = sd((log(Stime)+eta)[which(status==1)]) * (n^(-1/2))
      bw2 = sd(log(Stime)+eta) * (n^(-1/2))
    }
    else if(bw == "1/5"){
      bw1 = sd((log(Stime)+eta)[which(status==1)]) * (n^(-1/5))
      bw2 = sd(log(Stime)+eta) * (n^(-1/5))
    }
    else if(bw == "1/7"){
      bw1 = sd((log(Stime)+eta)[which(status==1)]) * (n^(-1/7))
      bw2 = sd(log(Stime)+eta) * (n^(-1/7))
    }
    else if((bw == "1/9")){
      bw1 = sd((log(Stime)+eta)[which(status==1)]) * (n^(-1/9))
      bw2 = sd(log(Stime)+eta) * (n^(-1/9))
    }
    temp1 = temp2 = 0
    for(ii in 1:length(Stime)){
      temp1 = temp1 + (status[ii]*dnorm((r[ii]-log(u))/bw1)/bw1)
      temp2 = temp2 + pnorm((r[ii]-log(u))/bw2)*u
    }
    return(as.numeric(temp1/temp2))
  }
  hazard_u <- sapply(1:length(eta),function(k){
    haz.fun(find[k])
  })
  at.risk <- ((Stime >=Target)&(entry <=Target))
  the.eta <- eta[ at.risk ]   #風 險 集 合 中 的 共 變 數 值
  the.haz <- hazard_u[ at.risk ] #風 險 集 合 中 的 基 準 風 險 函 數
  n <- length(the.eta)
  the.dead <- (Stime==Target)&(status==1)
  the.dead <- the.dead[ at.risk ]
  n.dead <- sum(the.dead)
  p0 <- rep( 1/(n-n.dead), n )
  if( n.dead > 0 ) {
    p0[ the.dead ] <- 0
  }
  p1 <- exp(the.eta)*the.haz/sum(exp(the.eta)*the.haz)
  ooo <- order(the.eta)
  eta.out <- the.eta[ooo]
  TP <- c(1, 1 - cumsum(p1[ooo]), 0)  # True Positive
  FP <- c(1, 1 - cumsum(p0[ooo]), 0)  # False Positive
  dFP <- abs(FP[-1] - FP[-length(FP)])
  aTP <- 0.5 * (TP[-1] + TP[-length(TP)])
  area <- sum( dFP * aTP )
  out <- list( marker = eta.out , TP = TP, FP = FP, AUC = area )
  return(out)
}

##  修正 Kernel smooth 方法估計 AUC (修正後)
AFTWeights.1.new <- function(marker , Stime , status , predict.time , entry=NULL,bw = "opt"){
  Target <- predict.time
  eta <- marker
  if(length(entry) == 0){
    entry = rep(0,NROW(Stime))
  }
  u <- exp( eta)*Stime
  find <- exp( eta)*predict.time
  # 修正後的風險函數估計
  haz.fun <- function(u){
    r = log(Stime)+eta ; n <- length(eta)
    if(bw == "opt"){
      bw1 = ((8*sqrt(2)/3)^(0.2)) * sd((log(Stime)+eta)[which(status==1)]) * (n^(-0.2))
      bw2 = (4^(1/3)) * sd(log(Stime)+eta) * (n^(-1/3))
    }
    else if(bw == "1/2"){
      bw1 = sd((log(Stime)+eta)[which(status==1)]) * (n^(-1/2))
      bw2 = sd(log(Stime)+eta) * (n^(-1/2))
    }
    else if(bw == "1/5"){
      bw1 = sd((log(Stime)+eta)[which(status==1)]) * (n^(-1/5))
      bw2 = sd(log(Stime)+eta) * (n^(-1/5))
    }
    else if(bw == "1/7"){
      bw1 = sd((log(Stime)+eta)[which(status==1)]) * (n^(-1/7))
      bw2 = sd(log(Stime)+eta) * (n^(-1/7))
    }
    else if((bw == "1/9")){
      bw1 = sd((log(Stime)+eta)[which(status==1)]) * (n^(-1/9))
      bw2 = sd(log(Stime)+eta) * (n^(-1/9))
    }
    temp1 = temp2 = 0
    for(ii in 1:length(Stime)){
      temp1 = temp1 + (status[ii]*dnorm((r[ii]-log(u))/bw1)/bw1)
      temp2 = temp2 + pnorm((r[ii]-log(u))/bw2)*u
    }
    return(as.numeric(temp1/temp2))
  }
  hazard_u <- sapply(1:length(eta),function(k){
    haz.fun(find[k])
  })
  at.risk <- ((Stime >=Target)&(entry <=Target))
  the.eta <- eta[ at.risk ]   #風 險 集 合 中 的 共 變 數 值
  the.haz <- hazard_u[ at.risk ] #風 險 集 合 中 的 基 準 風 險 函 數
  n <- length(the.eta)
  the.dead <- (Stime==Target)&(status==1)
  the.dead <- the.dead[ at.risk ]
  n.dead <- sum(the.dead)
  p0 <- rep( 1/(n-n.dead), n )
  if( n.dead > 0 ) {
    p0[ the.dead ] <- 0
  }
  p1 <- exp(the.eta)*the.haz/sum(exp(the.eta)*the.haz)
  ooo <- order(the.eta)
  eta.out <- the.eta[ooo]
  TP <- c(1, 1 - cumsum(p1[ooo]), 0)  # True Positive
  FP <- c(1, 1 - cumsum(p0[ooo]), 0)  # False Positive
  a = which(FP[-1] == FP[-length(FP)])
  a = a + 1 
  TP = TP[-a]
  FP = FP[-a]
  #plot(FP,TP)
  dFP <- abs(FP[-1] - FP[-length(FP)])
  aTP <- 0.5 * (TP[-1] + TP[-length(TP)])
  area <- sum( dFP * aTP )
  out <- list( marker = eta.out , TP = TP, FP = FP, AUC = area )
  return(out)
}

##  原始 Kernel smooth 方法估計 AUC (舊版)
AFTWeights.2 <- function(marker , Stime , status , predict.time , entry=NULL,bw.mul=5){
  Target <- predict.time
  eta <- marker
  if(length(entry) == 0){
    entry = rep(0,NROW(Stime))
  }
  find <- exp( eta)*predict.time
  u <- exp( eta)*Stime
  fit.mu.1 <- muhaz(u,status ,n.est.grid=10000, bw.method="g")
  bw <- fit.mu.1$pin$bw.smooth
  newBW.1 <- bw*bw.mul
  fit.mu.2 <- muhaz(u,status ,n.est.grid=10000, bw.smooth = newBW.1)
  hazard_u <- sapply(1:length(Stime),function(i) fit.mu.2$haz.est[which.min(abs(fit.mu.2$est.grid - find[i]))] )
  at.risk <- ((Stime >=Target)&(entry <=Target))
  the.eta <- eta[ at.risk ]   #風 險 集 合 中 的 共 變 數 值
  the.haz <- hazard_u[ at.risk ] #風 險 集 合 中 的 基 準 風 險 函 數
  n <- length(the.eta)
  the.dead <- (Stime==Target)&(status==1)
  the.dead <- the.dead[ at.risk ]
  n.dead <- sum(the.dead)
  p0 <- rep( 1/(n-n.dead), n )
  if( n.dead > 0 ) {
    p0[ the.dead ] <- 0
  }
  p1 <- exp(the.eta)*the.haz/sum(exp(the.eta)*the.haz)
  ooo <- order(the.eta)
  eta.out <- the.eta[ooo]
  TP <- c(1, 1 - cumsum(p1[ooo]), 0)  # True Positive
  FP <- c(1, 1 - cumsum(p0[ooo]), 0)  # False Positive
  dFP <- abs(FP[-1] - FP[-length(FP)])
  aTP <- 0.5 * (TP[-1] + TP[-length(TP)])
  area <- sum( dFP * aTP )
  out <- list( marker = eta.out , TP = TP, FP = FP, AUC = area )
  return(out)
}

##  原始 Kernel smooth 方法估計 AUC (修正後)
AFTWeights.2.new <- function(marker , Stime , status , predict.time , entry=NULL,bw.mul=5){
  Target <- predict.time
  eta <- marker
  if(length(entry) == 0){
    entry = rep(0,NROW(Stime))
  }
  find <- exp( eta)*predict.time
  u <- exp( eta)*Stime
  fit.mu.1 <- muhaz(u,status ,n.est.grid=10000, bw.method="g")
  bw <- fit.mu.1$pin$bw.smooth
  newBW.1 <- bw*bw.mul
  fit.mu.2 <- muhaz(u,status ,n.est.grid=10000, bw.smooth = newBW.1)
  hazard_u <- sapply(1:length(Stime),function(i) fit.mu.2$haz.est[which.min(abs(fit.mu.2$est.grid - find[i]))] )
  at.risk <- ((Stime >=Target)&(entry <=Target))
  the.eta <- eta[ at.risk ]   #風 險 集 合 中 的 共 變 數 值
  the.haz <- hazard_u[ at.risk ] #風 險 集 合 中 的 基 準 風 險 函 數
  n <- length(the.eta)
  the.dead <- (Stime==Target)&(status==1)
  the.dead <- the.dead[ at.risk ]
  n.dead <- sum(the.dead)
  p0 <- rep( 1/(n-n.dead), n )
  if( n.dead > 0 ) {
    p0[ the.dead ] <- 0
  }
  p1 <- exp(the.eta)*the.haz/sum(exp(the.eta)*the.haz)
  ooo <- order(the.eta)
  eta.out <- the.eta[ooo]
  TP <- c(1, 1 - cumsum(p1[ooo]), 0)  # True Positive
  FP <- c(1, 1 - cumsum(p0[ooo]), 0)  # False Positive
  a = which(FP[-1] == FP[-length(FP)])
  a = a + 1 
  TP = TP[-a]
  FP = FP[-a]
  #plot(FP,TP)
  dFP <- abs(FP[-1] - FP[-length(FP)])
  aTP <- 0.5 * (TP[-1] + TP[-length(TP)])
  area <- sum( dFP * aTP )
  out <- list( marker = eta.out , TP = TP, FP = FP, AUC = area )
  return(out)
}

#### PO 模型估計 AUC (舊版)
POWeights <- function(marker , Stime , status , predict.time , entry=NULL){
  if(length(entry) == 0){
    entry = rep(0,NROW(Stime))
  }
  Target = predict.time
  eta = marker
  data = data.frame(time=entry,time2=Stime,status = status,eta = eta)
  fit.PO = prop.odds(Event(time2,status)~eta,max.time=max(Stime),data = data,n.sim=100)
  uncen_time = fit.PO$cum[,1]
  cum_odds_0_uncen = fit.PO$cum[,2]  
  Gt = cum_odds_0_uncen[which.min(abs(Target-uncen_time))]  # 基線勝算估計
  at.risk <- ((Stime >= Target)&(entry <= Target))
  the.eta <- marker[ at.risk ]   # 風險集合中的共變數值
  n <- length(the.eta)
  the.dead <- (Stime==Target)&(status==1)
  # 在預測時間上發生事件的個數
  the.dead <- the.dead[ at.risk ]
  n.dead <- sum(the.dead)
  p0 <- rep( 1/(n-n.dead), n )
  if( n.dead > 0 ) {
    p0[ the.dead ] <- 0
  }
  p1 <- (1/(exp(-the.eta)+Gt))/sum(1/(exp(-the.eta)+Gt))
  ooo <- order(the.eta)
  eta.out <- the.eta[ooo]
  TP <- c(1, 1 - cumsum(p1[ooo]), 0)  # True Positive
  FP <- c(1, 1 - cumsum(p0[ooo]), 0)  # False Positive
  dFP <- abs(FP[-1] - FP[-length(FP)])
  aTP <- 0.5 * (TP[-1] + TP[-length(TP)])
  area <- sum( dFP * aTP )
  out <- list( marker = eta.out , TP = TP, FP = FP, AUC = area)
  return(out)
}

#### PO 模型估計 AUC (修正後)
POWeights.new <- function(marker , Stime , status , predict.time , entry=NULL){
  if(length(entry) == 0){
    entry = rep(0,NROW(Stime))
  }
  Target = predict.time
  eta = marker
  data = data.frame(time=entry,time2=Stime,status = status,eta = eta)
  fit.PO = prop.odds(Event(time2,status)~eta,max.time=max(Stime),data = data,n.sim=100)
  uncen_time = fit.PO$cum[,1]
  cum_odds_0_uncen = fit.PO$cum[,2]  
  Gt = cum_odds_0_uncen[which.min(abs(Target-uncen_time))]  # 基線勝算估計
  at.risk <- ((Stime >= Target)&(entry <= Target))
  the.eta <- marker[ at.risk ]   # 風險集合中的共變數值
  n <- length(the.eta)
  the.dead <- (Stime==Target)&(status==1)
  # 在預測時間上發生事件的個數
  the.dead <- the.dead[ at.risk ]
  n.dead <- sum(the.dead)
  p0 <- rep( 1/(n-n.dead), n )
  if( n.dead > 0 ) {
    p0[ the.dead ] <- 0
  }
  p1 <- (1/(exp(-the.eta)+Gt))/sum(1/(exp(-the.eta)+Gt))
  ooo <- order(the.eta)
  eta.out <- the.eta[ooo]
  TP <- c(1, 1 - cumsum(p1[ooo]), 0)  # True Positive
  FP <- c(1, 1 - cumsum(p0[ooo]), 0)  # False Positive
  a = which(FP[-1] == FP[-length(FP)])
  a = a + 1 
  TP = TP[-a]
  FP = FP[-a]
  #plot(FP,TP)
  dFP <- abs(FP[-1] - FP[-length(FP)])
  aTP <- 0.5 * (TP[-1] + TP[-length(TP)])
  area <- sum( dFP * aTP )
  out <- list( marker = eta.out , TP = TP, FP = FP, AUC = area)
  return(out)
}

#### 計算各時間下 AUC 與 Concordance (舊版)
MyrisksetAUC.wai.old <- function( Stime , entry = NULL, status, marker, method 
                              , tmax , AFT.method.1 = "1",AFT.method.2 = "1"
                              , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                              , weight = "rescale", plot=TRUE, type="l",xlab="Time", ylab="AUC"){
  # Stime:受試者觀測時間
  # entry:受試者進入試驗時間
  # status:受試者設限指標
  # marker:共變數以迴歸係數加權後總和
  # method:選擇模型 (Cox、AFT、PO)
  # tmax: 估計 Concordance 的最大時間範圍
  # 若使用AFT模型分析，需選擇是否使用新的 kernel smooth 方法，"1":修正後，"2":2019方法
  # AFT.method.1: 點估計方法選擇 ，AFT.method.2: 信賴區間估計方法選擇 
  # 新的 kernel smooth 方法需要給定帶寬 bw.1、bw.2 預設為"opt"，可設定"opt"、"1/2"、"1/5"、"1/7"、"1/9"
  # 原始方法需要給定帶寬 bw.mul.1、bw.mul.2，預設為 5
  # CI:是否需要區間估計，預設為 F ， confi: 信賴區間信心水準
  # weight: 一致性指標加權方式，可選 "rescale"、"conditional" 
  # plot: 是否需要劃出 AUC 變化趨勢
  mChoice <- match(method , c("Cox", "AFT","PO")) # 選取使用的模型
  if (is.na(mChoice)) {
    cat("error_in_method_choice")
    stop(0) 
  }
  if (length(entry) == 0) {
    time = rep(0, length(Stime))
  }
  else {
    time = entry
  }
  eta <- marker
  time2 = Stime
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  km.out = weightedKM(Stime = Stime , status = status , entry = time)#計算 survival function
  AUC <- rep( NA, length(utimes) )
  # 計算出事件發生時間下的 AUC
  for (i in 1:length(utimes)) {
    new.eta <- marker
    if (mChoice == 1 ){ #選用Cox模型
      out <- CoxWeights(marker=new.eta, Stime=Stime, status=status,
                        predict.time=utimes[i], entry=entry)
    }
    if (mChoice == 2 ){ # 選用AFT模型  
      if(AFT.method.1=="1"){ # 使用修正後的 kernel smooth 方法
        out <- AFTWeights.1(new.eta , Stime , status , utimes[i],bw = bw.1, entry=entry)
      }
      else if(AFT.method.1=="2"){   # 使用原始的 kernel smooth 方法 (林冠廷,2019)
        out <- AFTWeights.2(new.eta , Stime , status , utimes[i],bw.mul = bw.mul.1, entry=entry)
      }
    }
    if (mChoice == 3 ){ #選用PO模型
      out <- POWeights(new.eta , Stime , status , utimes[i], entry=entry)
    }
    AUC[i] <- out$AUC
  }
  Cindex = IntegrateAUC(AUC , utimes , km.out$survival , tmax, weight = weight) # concordance 估計
  confi_value = abs(qnorm((1-confi)/2))
  if(CI == T){ # AUC 與 Concordance 標準差估計
    Std_all = EST.cindex.std(marker=eta, Stime=Stime, status=status, utimes=utimes,
                             St = km.out$survival, method = method,AFT.method = AFT.method.2, tmax, 
                             weight = weight, entry=entry, bw=bw.2, bw.mul=bw.mul.2)
    AUC_std <- Std_all$AUC_std ; Cindex_std <- Std_all$Cindex_std
  }else{
    AUC_std <- 0 ; Cindex_std <- 0
  }
  AUC_low <- AUC-confi_value*AUC_std
  AUC_upp <- AUC+confi_value*AUC_std
  AUC_confi = rbind(AUC_low,AUC_upp)
  Cindex_confi = c(Cindex-confi_value*Cindex_std,Cindex+confi_value*Cindex_std)
  if(plot==TRUE)
  {
    # plot(utimes, AUC, type=type, xlim=c(min(utimes), tmax+1),
    #      ylim=c(0.4,1.0), xlab=xlab, ylab=ylab,...)
    plot(utimes, AUC, type=type, xlim=c(0, 10),
         ylim=c(0.0,1.0), xlab=xlab, ylab=ylab,lwd = 2)
    abline(h=0.5)
  }
  #計 算 出 在 tmax 時 間 點 下 的 concordance
  return(out = list(utimes = utimes , St = km.out$survival ,
                    AUC = AUC, AUC_std = AUC_std, AUC_confi = AUC_confi ,
                    Cindex.old = Cindex, Cindex_std = Cindex_std, Cindex_confi = Cindex_confi
  ))
}

#### 計算各時間下 AUC 與 Concordance (修正後)
MyrisksetAUC.wai.new <- function( Stime , entry = NULL, status, marker, method 
                                  , tmax , AFT.method.1 = "1",AFT.method.2 = "1"
                                  , bw.1="opt",bw.2="opt", bw.mul.1=5,bw.mul.2=5, CI=F, confi=0.95
                                  , weight = "rescale", plot=TRUE, type="l",xlab="Time", ylab="AUC"){
  # Stime:受試者觀測時間
  # entry:受試者進入試驗時間
  # status:受試者設限指標
  # marker:共變數以迴歸係數加權後總和
  # method:選擇模型 (Cox、AFT、PO)
  # tmax: 估計 Concordance 的最大時間範圍
  # 若使用AFT模型分析，需選擇是否使用新的 kernel smooth 方法，"1":修正後，"2":2019方法
  # AFT.method.1: 點估計方法選擇 ，AFT.method.2: 信賴區間估計方法選擇 
  # 新的 kernel smooth 方法需要給定帶寬 bw.1、bw.2 預設為"opt"，可設定"opt"、"1/2"、"1/5"、"1/7"、"1/9"
  # 原始方法需要給定帶寬 bw.mul.1、bw.mul.2，預設為 5
  # CI:是否需要區間估計，預設為 F ， confi: 信賴區間信心水準
  # weight: 一致性指標加權方式，可選 "rescale"、"conditional" 
  # plot: 是否需要劃出 AUC 變化趨勢
  mChoice <- match(method , c("Cox", "AFT","PO")) # 選取使用的模型
  if (is.na(mChoice)) {
    cat("error_in_method_choice")
    stop(0) 
  }
  if (length(entry) == 0) {
    time = rep(0, length(Stime))
  }
  else {
    time = entry
  }
  eta <- marker
  time2 = Stime
  utimes = unique(Stime[status == 1]) # 事件發生時間
  utimes = utimes[order(utimes)]
  km.out = weightedKM(Stime = Stime , status = status , entry = time)#計算 survival function
  AUC <- rep( NA, length(utimes) )
  # 計算出事件發生時間下的 AUC
  for (i in 1:length(utimes)) {
    new.eta <- marker
    if (mChoice == 1 ){ #選用Cox模 型
      out <- CoxWeights.new(marker=new.eta, Stime=Stime, status=status,
                        predict.time=utimes[i], entry=entry)
    }
    if (mChoice == 2 ){ # 選用AFT模 型  
      if(AFT.method.1=="1"){ # 使用修正後的 kernel smooth 方法
        out <- AFTWeights.1.new(new.eta , Stime , status , utimes[i],bw = bw.1, entry=entry)
      }
      else if(AFT.method.1=="2"){   # 使用原始的 kernel smooth 方法 (林冠廷,2019)
        out <- AFTWeights.2.new(new.eta , Stime , status , utimes[i],bw.mul = bw.mul.1, entry=entry)
      }
    }
    if (mChoice == 3 ){ #選用PO模 型
      out <- POWeights.new(new.eta , Stime , status , utimes[i], entry=entry)
    }
    AUC[i] <- out$AUC
  }
  Cindex = IntegrateAUC(AUC , utimes , km.out$survival , tmax, weight = weight) # concordance 估計
  confi_value = abs(qnorm((1-confi)/2))
  if(CI == T){ # AUC 與 Concordance 標準差估計
    Std_all = EST.cindex.std(marker=eta, Stime=Stime, status=status, utimes=utimes,
                             St = km.out$survival, method = method,AFT.method = AFT.method.2, tmax, 
                             weight = weight, entry=entry, bw=bw.2, bw.mul=bw.mul.2)
    AUC_std <- Std_all$AUC_std ; Cindex_std <- Std_all$Cindex_std
  }else{
    AUC_std <- 0 ; Cindex_std <- 0
  }
  AUC_low <- AUC-confi_value*AUC_std
  AUC_upp <- AUC+confi_value*AUC_std
  AUC_confi = rbind(AUC_low,AUC_upp)
  Cindex_confi = c(Cindex-confi_value*Cindex_std,Cindex+confi_value*Cindex_std)
  if(plot==TRUE)
  {
    # plot(utimes, AUC, type=type, xlim=c(min(utimes), tmax+1),
    #      ylim=c(0.4,1.0), xlab=xlab, ylab=ylab,...)
    plot(utimes, AUC, type=type, xlim=c(0, 10),
         ylim=c(0.0,1.0), xlab=xlab, ylab=ylab,lwd = 2)
    abline(h=0.5)
  }
  #計 算 出 在 tmax 時 間 點 下 的 concordance
  return(out = list(utimes = utimes , St = km.out$survival ,
                    AUC = AUC, AUC_std = AUC_std, AUC_confi = AUC_confi ,
                    Cindex.new = Cindex, Cindex_std = Cindex_std, Cindex_confi = Cindex_confi
  ))
}



#### 計算各時間下 AUC 與 Concordance 的標準差
EST.cindex.std <- function(marker, Stime, status, utimes, St, 
                           method = "Cox",AFT.method = "1", tmax, 
                           weight = "rescale",entry = NULL, bw="opt", bw.mul=5){
  mChoice <- match(method , c("Cox", "AFT","PO"))
  wChoice <- match(weight, c("rescale", "conditional"))
  eta = marker
  ft <- rep(NA, length(St))
  ft[1] <- 1 - St[1]  
  for(j in 2:length(St)) ft[j] <- St[j - 1] - St[j]
  mIndex <- length(utimes[utimes <= tmax])
  www <- 2 * ft * St
  wTotal <- sum(www[1:mIndex])
  Smax <- St[min(mIndex + 1, length(St))]
  if (wChoice == 1) 
       w = 2 * ft * St/wTotal
  else w = 2 * ft * (St - Smax)/((1 - Smax)^2)
  if (length(entry) == 0) {
    entry = rep(0, NROW(Stime))
  }
  ## PO 模型需要基線勝算估計
  if(mChoice == 3){
    eta = marker
    data = data.frame(time=entry,time2=Stime,status = status,eta = eta)
    fit.PO = prop.odds(Event(time2,status)~eta,max.time=max(Stime),data = data,n.sim=100)
  }
  # 根據 2.3 節，將所有可能的 Y_{p_k}^{t_i} 組合列出
  Y_t_p  = lapply(1:length(utimes),function(i){
    Risk <- ifelse(((Stime >= utimes[i]) & (entry <= utimes[i])),1,0)
    at.risk <- ((Stime >= utimes[i]) & (entry <= utimes[i]))
    the.eta <- eta[at.risk]
    eta.out <- sort(the.eta)
    if(mChoice==1){
      Y_p.1 <- sapply(1:(length(eta.out)),function(j) {
        index = ifelse(eta >= eta.out[j],1,0)
        Y_star = (exp(eta))*Risk*index
        return(Y_star)
      })
    }
    else if(mChoice==2){
      find <- exp( eta)*utimes[i]
      if(AFT.method=="1"){
        haz.fun <- function(u){
          r = log(Stime)+eta ; n <- length(eta)
          if(bw == "opt"){
            bw1 = ((8*sqrt(2)/3)^(0.2)) * sd((log(Stime)+eta)[which(status==1)]) * (n^(-0.2))
            bw2 = (4^(1/3)) * sd(log(Stime)+eta) * (n^(-1/3))
          }
          else if(bw == "1/2"){
            bw1 = sd((log(Stime)+eta)[which(status==1)]) * (n^(-1/2))
            bw2 = sd(log(Stime)+eta) * (n^(-1/2))
          }
          else if(bw == "1/5"){
            bw1 = sd((log(Stime)+eta)[which(status==1)]) * (n^(-1/5))
            bw2 = sd(log(Stime)+eta) * (n^(-1/5))
          }
          else if(bw == "1/7"){
            bw1 = sd((log(Stime)+eta)[which(status==1)]) * (n^(-1/7))
            bw2 = sd(log(Stime)+eta) * (n^(-1/7))
          }
          else if((bw == "1/9")){
            bw1 = sd((log(Stime)+eta)[which(status==1)]) * (n^(-1/9))
            bw2 = sd(log(Stime)+eta) * (n^(-1/9))
          }
          temp1 = temp2 = 0
          for(ii in 1:length(Stime)){
            temp1 = temp1 + (status[ii]*dnorm((r[ii]-log(u))/bw1)/bw1)
            temp2 = temp2 + pnorm((r[ii]-log(u))/bw2)*u
          }
          return(as.numeric(temp1/temp2))
        }
        hazard_u <- sapply(1:length(eta),function(k)haz.fun(find[k]))
      }
      else if(AFT.method=="2"){
        u <- exp( eta)*Stime
        fit.mu.1 <- muhaz(u,status ,n.est.grid=10000, bw.method="g")
        bw <- fit.mu.1$pin$bw.smooth
        newBW.1 <- bw*bw.mul
        fit.mu.2 <- muhaz(u,status ,n.est.grid=10000, bw.smooth = newBW.1)
        hazard_u <- sapply(1:length(eta),function(k) fit.mu.2$haz.est[which.min(abs(fit.mu.2$est.grid - find[k]))] )
      }
      Y_p.1 <- sapply(1:(length(eta.out)),function(j) {
        index = ifelse(eta>=eta.out[j],1,0)
        Y_star = (exp(eta)*hazard_u)*Risk*index
        return(Y_star)
      })
    }
    else if(mChoice==3){
      uncen_time = fit.PO$cum[,1]
      cum_odds_0_uncen = fit.PO$cum[,2]
      Gt = cum_odds_0_uncen[which.min(abs(utimes[i]-uncen_time))]
      Y_p.1 <- sapply(1:(length(eta.out)),function(j) {
        index = ifelse(eta>=eta.out[j],1,0)
        Y_star = (1/(exp(-eta)+Gt))*Risk*index
        return(Y_star)
      })
    }
  })
  # 各時間點當下仍存活的個數
  L_in_time_i = sapply(1:length(utimes),function(i){dim(Y_t_p[[i]])[2]})
  # 為了取不同時間點與切點下 Y_{p_k}^{t_i} 的共變異數矩陣資訊，需要確認在共變異數矩陣中的位置
  start       = c(1,(cumsum(L_in_time_i)+1)[-length(L_in_time_i)])
  end         = cumsum(L_in_time_i)
  # 將 Y_t_p 轉為向量表示
  Y_t_p = matrix(as.numeric(unlist(Y_t_p)),length(eta))
  # 使用樣本平均數與樣本共變異數矩陣替代未知的母體平均數與母體共變異數矩陣
  mu_all = colMeans(Y_t_p) ; Sigma_all = cov(Y_t_p)  
  ## 生成使用 Delta method 所需的 g'(u)
  gu     = matrix(0,length(mu_all),length(mu_all))
  for(i in 1:length(utimes)){
    gu_sub = matrix(0,L_in_time_i[i],L_in_time_i[i])
    diag(gu_sub) = 1/mu_all[start[i]]
    if( i == length(utimes)){
      gu_sub[,1] = -mu_all[start[i]]/(mu_all[start[i]]^2)
    }else{
      gu_sub[,1] = -mu_all[start[i]:end[i]]/(mu_all[start[i]]^2)
    }
    gu_sub[1,1] = 0 ##
    gu[start[i]:end[i],start[i]:end[i]] = gu_sub
  }
  d          = sapply(1:length(utimes),function(i){
    at.risk <- ((Stime >= utimes[i]) & (entry <= utimes[i]))
    the.eta <- eta[at.risk] 
    n <- length(the.eta) 
    the.dead <- (Stime == utimes[i]) & (status == 1)
    the.dead <- the.dead[at.risk]
    n.dead <- sum(the.dead)
    d.i <- 1/(n - n.dead)
    min(1,d.i)
  })
  Sigma_time = gu %*% Sigma_all %*% t(gu)
  ## 估計各時間點的AUC標準差
  AUC_std    = sapply(1:length(utimes),function(i){
    sqrt(sum(d[i]*d[i]*Sigma_time[start[i]:end[i],start[i]:end[i]])/length(eta))
  })
  ## 估計 Concordance 標準差
  ccc = 0
  for(i in 1:length(utimes[utimes <= tmax])){
  for(j in 1:length(utimes[utimes <= tmax])){
    ccc = ccc + sum(w[i]*w[j]*d[i]*d[j]*Sigma_time[start[i]:end[i],start[j]:end[j]])/length(eta)
  }
  }
  Cindex_std = sqrt(ccc)
  return(list(Cindex_std = Cindex_std ,AUC_std = AUC_std))
}
# 利用長期追蹤模型建構完整共變量歷史
findding_eta <- function(predict.time, Stime, bi, id){
  Mt = sapply(1:length(unique(id)),function(i){
    time = min(predict.time,Stime[i])
    vec.time = sapply(1:dim(bi)[2],function(j){
      time^(j-1)
    })
    bi[i,]%*%t(t(vec.time))
  })
  return(as.numeric(Mt))
}

###### True value 
#### AUC
## Weibull
Cox.dist.Weibull.true <- function (predict.time , alpha , lambda, beta, mu.z=0){
  Target = predict.time
  f1 = function(x){exp(beta*x-lambda*(Target^alpha)*exp(beta*x)-0.5*(x-mu.z)^2)}
  f2 = function(x){exp(-0.5*(x-mu.z)^2-lambda*(Target^alpha)*exp(beta*x))}
  c = rnorm(1000)
  FPR = rep(0,1000);TPR = rep(0,1000)
  marker = sort(c)
  for (i in 1:1000) {
    if(beta>0){
      tp.top = integrate(f1,marker[i]/beta,Inf)$value
      fp.top = integrate(f2,marker[i]/beta,Inf)$value
    }else{
      tp.top = integrate(f1,-Inf,marker[i]/beta)$value
      fp.top = integrate(f2,-Inf,marker[i]/beta)$value
    }
    tp.low = integrate(f1,-Inf,Inf)$value
    fp.low = integrate(f2,-Inf,Inf)$value
    TPR[i] = tp.top/tp.low
    FPR[i] = fp.top/fp.low
  }
  TPR = c(1,TPR,0) ;FPR = c(1,FPR,0) 
  f=approxfun(FPR,TPR)
  area <- integrate(f,0,1)$value
  plot(FPR,TPR)
  out = list(AUC = area,
             TPR = TPR,
             FPR = FPR)
  return(out)
}

AFT.dist.Weibull.true <- function (predict.time , alpha , lambda, beta, mu.z=0){
  Target = predict.time
  f1 = function(x){exp(alpha*beta*x-(lambda*Target^alpha)*exp(alpha*beta*x)-0.5*(x-mu.z)^2)}
  f2 = function(x){exp(-0.5*(x-mu.z)^2-lambda*(Target^alpha)*exp(alpha*beta*x))}
  c = rnorm(1000)
  FPR = rep(0,1000);TPR = rep(0,1000)
  marker = sort(c)
  for (i in 1:1000) {
    if(beta>0){
      tp.top = integrate(f1,marker[i]/beta,Inf)$value
      fp.top = integrate(f2,marker[i]/beta,Inf)$value
    }else{
      tp.top = integrate(f1,-Inf,marker[i]/beta)$value
      fp.top = integrate(f2,-Inf,marker[i]/beta)$value
    }
    tp.low = integrate(f1,-Inf,Inf)$value
    fp.low = integrate(f2,-Inf,Inf)$value
    TPR[i] = tp.top/tp.low
    FPR[i] = fp.top/fp.low
  }
  TPR = c(1,TPR,0) ;FPR = c(1,FPR,0) 
  f=approxfun(FPR,TPR)
  area <- integrate(f,0,1)$value
  plot(FPR,TPR)
  out = list(AUC = area)
  return(out)
}

PO.dist.Weibull.true <- function (predict.time , alpha , lambda, beta, mu.z=0){
  Target = predict.time
  f1 = function(x){exp(beta*x-0.5*(x-mu.z)^2)/((1+(exp(lambda*Target^alpha)-1)*exp(beta*x))^2)}
  f2 = function(x){exp(-0.5*(x-mu.z)^2)/(1+(exp(lambda*Target^alpha)-1)*exp(beta*x))}
  c = rnorm(1000)
  FPR = rep(0,1000);TPR = rep(0,1000)
  marker = sort(c)
  for (i in 1:1000) {
    if(beta>0){
      tp.top = integrate(f1,marker[i]/beta,Inf)$value
      fp.top = integrate(f2,marker[i]/beta,Inf)$value
    }else{
      tp.top = integrate(f1,-Inf,marker[i]/beta)$value
      fp.top = integrate(f2,-Inf,marker[i]/beta)$value
    }
    tp.low = integrate(f1,-Inf,Inf)$value
    fp.low = integrate(f2,-Inf,Inf)$value
    TPR[i] = tp.top/tp.low
    FPR[i] = fp.top/fp.low
  }
  TPR = c(1,TPR,0) ;FPR = c(1,FPR,0) 
  f=approxfun(FPR,TPR)
  area <- integrate(f,0,1)$value
  plot(FPR,TPR)
  out = list(AUC = area)
  return(out)
}
## Lognormal
Cox.dist.lognormal.true <- function(predict.time , mu ,sigma, beta, mu.z=0){
  Target = predict.time
  f1 = function(x){exp(beta*x-0.5*(x-mu.z)^2)*(1-pnorm((log(Target)-mu)/sigma))^(exp(beta*x))}
  f2 = function(x){exp(-0.5*(x-mu.z)^2)*(1-pnorm((log(Target)-mu)/sigma))^(exp(beta*x))}
  c = rnorm(1000,0,1)
  FPR = rep(0,1000);TPR = rep(0,1000)
  marker = sort(c)
  for (i in 1:1000) {
    if(beta>0){
      tp.top = integrate(f1,marker[i]/beta,Inf)$value
      fp.top = integrate(f2,marker[i]/beta,Inf)$value
    }else{
      tp.top = integrate(f1,-Inf,marker[i]/beta)$value
      fp.top = integrate(f2,-Inf,marker[i]/beta)$value
    }
    tp.low = integrate(f1,-Inf,Inf)$value
    fp.low = integrate(f2,-Inf,Inf)$value
    TPR[i] = tp.top/tp.low
    FPR[i] = fp.top/fp.low
  }
  TPR = c(1,TPR,0) ;FPR = c(1,FPR,0) 
  f=approxfun(FPR,TPR)
  area <- integrate(f,0,1)$value
  plot(FPR,TPR)
  out = list(AUC = area)
  return(out)
}

AFT.dist.lognormal.true <- function(predict.time , mu ,sigma, beta, mu.z=0){
  Target = predict.time
  f1 = function(x){exp(-0.5*(x-mu.z)^2)*dnorm((log(Target)+beta*x-mu)/sigma)}
  f2 = function(x){exp(-0.5*(x-mu.z)^2)*(1-pnorm((log(Target)+beta*x-mu)/sigma))}
  c = rnorm(1000)
  FPR = rep(0,1000);TPR = rep(0,1000)
  marker = sort(c)
  for (i in 1:1000) {
    if(beta>0){
      tp.top = integrate(f1,marker[i]/beta,Inf)$value
      fp.top = integrate(f2,marker[i]/beta,Inf)$value
    }else{
      tp.top = integrate(f1,-Inf,marker[i]/beta)$value
      fp.top = integrate(f2,-Inf,marker[i]/beta)$value
    }
    tp.low = integrate(f1,-Inf,Inf)$value
    fp.low = integrate(f2,-Inf,Inf)$value
    TPR[i] = tp.top/tp.low
    FPR[i] = fp.top/fp.low
  }
  TPR = c(1,TPR,0) ;FPR = c(1,FPR,0) 
  f=approxfun(FPR,TPR)
  area <- integrate(f,0,1)$value
  plot(FPR,TPR)
  out = list(AUC = area)
  return(out)
}

PO.dist.lognormal.true <- function(predict.time , mu ,sigma, beta, mu.z=0){
  Target = predict.time
  f1 = function(x){exp(beta*x-0.5*(x-mu.z)^2)/((1+(pnorm((log(Target)-mu)/sigma)/(1-pnorm((log(Target)-mu)/sigma)))*exp(beta*x))^2)}
  f2 = function(x){exp(-0.5*(x-mu.z)^2)/(1+(pnorm((log(Target)-mu)/sigma)/(1-pnorm((log(Target)-mu)/sigma)))*exp(beta*x))}
  c = rnorm(1000)
  FPR = rep(0,1000);TPR = rep(0,1000)
  marker = sort(c)
  for (i in 1:1000) {
    if(beta>0){
      tp.top = integrate(f1,marker[i]/beta,Inf)$value
      fp.top = integrate(f2,marker[i]/beta,Inf)$value
    }else{
      tp.top = integrate(f1,-Inf,marker[i]/beta)$value
      fp.top = integrate(f2,-Inf,marker[i]/beta)$value
    }
    tp.low = integrate(f1,-Inf,Inf)$value
    fp.low = integrate(f2,-Inf,Inf)$value
    TPR[i] = tp.top/tp.low
    FPR[i] = fp.top/fp.low
  }
  TPR = c(1,TPR,0) ;FPR = c(1,FPR,0) 
  f=approxfun(FPR,TPR)
  area <- integrate(f,0,1)$value
  plot(FPR,TPR)
  out = list(AUC = area)
  return(out)
}
## loglogistic
Cox.dist.loglogistic.true <- function(predict.time , mu ,sigma, beta, mu.z=0){
  Target = predict.time
  f1 = function(x){exp(-0.5*(x-mu.z)^2+beta*x)*(1+exp((log(Target)-mu)/sigma))^(-exp(beta*x))}
  f2 = function(x){exp(-0.5*(x-mu.z)^2)*(1+exp((log(Target)-mu)/sigma))^(-exp(beta*x))}
  c = rnorm(1000,0,1)
  FPR = rep(0,1000);TPR = rep(0,1000)
  marker = sort(c)
  for (i in 1:1000) {
    if(beta>0){
      tp.top = integrate(f1,marker[i]/beta,Inf)$value
      fp.top = integrate(f2,marker[i]/beta,Inf)$value
    }else{
      tp.top = integrate(f1,-Inf,marker[i]/beta)$value
      fp.top = integrate(f2,-Inf,marker[i]/beta)$value
    }
    tp.low = integrate(f1,-Inf,Inf)$value
    fp.low = integrate(f2,-Inf,Inf)$value
    TPR[i] = tp.top/tp.low
    FPR[i] = fp.top/fp.low
  }
  TPR = c(1,TPR,0) ;FPR = c(1,FPR,0) 
  f=approxfun(FPR,TPR)
  area <- integrate(f,0,1)$value
  plot(FPR,TPR)
  out = list(AUC = area)
  return(out)
}

AFT.dist.loglogistic.true <- function(predict.time , mu ,sigma, beta, mu.z=0){
  Target = predict.time
  f1 = function(x){exp(-0.5*(x-mu.z)^2)*exp((log(Target)+beta*x-mu)/sigma)/((1+exp((log(Target)+beta*x-mu)/sigma))^2)}
  f2 = function(x){exp(-0.5*(x-mu.z)^2)/(1+exp((log(Target)+beta*x-mu)/sigma))}
  c = runif(1000,-20,20)
  FPR = rep(0,1000);TPR = rep(0,1000)
  marker = sort(c)
  for (i in 1:1000) {
    if(beta>0){
      tp.top = integrate(f1,marker[i]/beta,Inf)$value
      fp.top = integrate(f2,marker[i]/beta,Inf)$value
    }else{
      tp.top = integrate(f1,-Inf,marker[i]/beta)$value
      fp.top = integrate(f2,-Inf,marker[i]/beta)$value
    }
    tp.low = integrate(f1,-Inf,Inf)$value
    fp.low = integrate(f2,-Inf,Inf)$value
    TPR[i] = tp.top/tp.low
    FPR[i] = fp.top/fp.low
  }
  TPR = c(1,TPR,0) ;FPR = c(1,FPR,0) 
  f=approxfun(FPR,TPR)
  area <- integrate(f,0,1)$value
  plot(FPR,TPR)
  out = list(AUC = area)
  return(out)
}

PO.dist.loglogistic.true  <- function(predict.time , mu ,sigma, beta, mu.z=0){
  Target = predict.time
  f1 = function(x){exp(beta*x-0.5*(x-mu.z)^2)/((1+(exp((log(Target)-mu)/sigma))*exp(beta*x))^2)}
  f2 = function(x){exp(-0.5*(x-mu.z)^2)/(1+(exp((log(Target)-mu)/sigma))*exp(beta*x))}
  c = rnorm(1000)
  FPR = rep(0,1000);TPR = rep(0,1000)
  marker = sort(c)
  for (i in 1:1000) {
    if(beta>0){
      tp.top = integrate(f1,marker[i]/beta,Inf)$value
      fp.top = integrate(f2,marker[i]/beta,Inf)$value
    }else{
      tp.top = integrate(f1,-Inf,marker[i]/beta)$value
      fp.top = integrate(f2,-Inf,marker[i]/beta)$value
    }
    tp.low = integrate(f1,-Inf,Inf)$value
    fp.low = integrate(f2,-Inf,Inf)$value
    TPR[i] = tp.top/tp.low
    FPR[i] = fp.top/fp.low
  }
  TPR = c(1,TPR,0) ;FPR = c(1,FPR,0) 
  f=approxfun(FPR,TPR)
  area <- integrate(f,0,1)$value
  plot(FPR,TPR)
  out = list(AUC = area)
  return(out)
}
#### Cindex
## Weibull
MyrisksetAUC.Weibull.True <- function(alpha,lambda,beta, mu.z=0,tau,method = "Cox"){
  u = runif(1000)
  z <- rnorm(1000,mu.z,1)
  if(method=="Cox"){
    t0 = (-log(u)/(lambda*exp(beta*z)))^(1/alpha)
  }
  else if(method == "AFT"){
    t0 = (-log(u)/(lambda*exp(alpha*beta*z)))^(1/alpha)
  }
  else if(method == "PO"){
    t0 = ((log((1/u)-1+exp(beta*z))-beta*z)/lambda)^(1/alpha)
  }
  tmax = quantile(t0,tau)
  t0 = t0[which(t0 <= tmax)]
  AUC = rep(0,length(t0))
  for(k in 1:length(t0)){
    if(method=="Cox"){
      f1 = function(x){exp(beta*x-lambda*(t0[k]^alpha)*exp(beta*x)-0.5*(x-mu.z)^2)}
      f2 = function(x){exp(-0.5*(x-mu.z)^2-lambda*(t0[k]^alpha)*exp(beta*x))}
    }
    else if(method == "AFT"){
      f1 = function(x){exp(alpha*beta*x-(lambda*t0[k]^alpha)*exp(alpha*beta*x)-0.5*(x-mu.z)^2)}
      f2 = function(x){exp(-0.5*(x-mu.z)^2-lambda*(t0[k]^alpha)*exp(alpha*beta*x))}
    }
    else if(method == "PO"){
      f1 = function(x){exp(beta*x-0.5*x^2)/((1+(exp(lambda*t0[k]^alpha)-1)*exp(beta*x))^2)}
      f2 = function(x){exp(-0.5*x^2)/(1+(exp(lambda*t0[k]^alpha)-1)*exp(beta*x))}
    }
    c = rnorm(1000)
    FPR = rep(0,1000);TPR = rep(0,1000)
    marker = sort(c)
    for(i in 1:1000){
      if(beta>0){
        tp.top = integrate(f1,marker[i]/beta,Inf)$value
        fp.top = integrate(f2,marker[i]/beta,Inf)$value
      }else{
        tp.top = integrate(f1,-Inf,marker[i]/beta)$value
        fp.top = integrate(f2,-Inf,marker[i]/beta)$value
      }
      tp.low = integrate(f1,-Inf,Inf)$value
      fp.low = integrate(f2,-Inf,Inf)$value
      TPR[i] = tp.top/tp.low
      FPR[i] = fp.top/fp.low
    }
    TPR = c(1,TPR,0) ;FPR = c(1,FPR,0) 
    f = approxfun(FPR,TPR)
    AUC[k] = integrate(f,0,1)$value
  }
  AUC.weight = AUC*2*exp(-lambda*t0^alpha)/(1-exp(-2*lambda*tmax^alpha))
  Cindex = mean(AUC.weight)
  return(out = list(Cindex = Cindex))
}
## Lognormal
MyrisksetAUC.lognormal.True <- function(mu,sigma,beta, mu.z=0,tau,method = "Cox"){
  u = runif(1000)
  z <- rnorm(1000,mu.z,1)
  if(method=="Cox"){
    t0 = exp(sigma*qnorm(1-u^(exp(-beta*z)))+ mu)
  }
  else if(method == "AFT"){
    t0 = exp(sigma*qnorm(1-u)-beta*z+mu)
  }
  else if(method == "PO"){
    t0 = exp(mu+sigma*qnorm((1-(1-(1/u))*exp(beta*z))^(-1)))
  }
  tmax = quantile(t0,tau)
  t0 = t0[which(t0 <= tmax)]
  AUC = rep(0,length(t0))
  for(k in 1:length(t0)){
    if(method=="Cox"){
      f1 = function(x){exp(beta*x-0.5*(x-mu.z)^2)*(1-pnorm((log(t0[k])-mu)/sigma))^(exp(beta*x))}
      f2 = function(x){exp(-0.5*(x-mu.z)^2)*(1-pnorm((log(t0[k])-mu)/sigma))^(exp(beta*x))}
    }
    else if(method == "AFT"){
      f1 = function(x){exp(-0.5*(x-mu.z)^2)*dnorm((log(t0[k])+beta*x-mu)/sigma)}
      f2 = function(x){exp(-0.5*(x-mu.z)^2)*(1-pnorm((log(t0[k])+beta*x-mu)/sigma))}
    }
    else if(method == "PO"){
      f1 = function(x){exp(beta*x-0.5*(x-mu.z)^2)/((1+(pnorm((log(t0[k])-mu)/sigma)/(1-pnorm((log(t0[k])-mu)/sigma)))*exp(beta*x))^2)}
      f2 = function(x){exp(-0.5*(x-mu.z)^2)/(1+(pnorm((log(t0[k])-mu)/sigma)/(1-pnorm((log(t0[k])-mu)/sigma)))*exp(beta*x))}
    }
    c = rnorm(1000,0,1)
    FPR = rep(0,1000);TPR = rep(0,1000)
    marker = sort(c)
    for(i in 1:length(marker)){
      if(beta>0){
        tp.top = integrate(f1,marker[i]/beta,Inf)$value
        fp.top = integrate(f2,marker[i]/beta,Inf)$value
      }else{
        tp.top = integrate(f1,-Inf,marker[i]/beta)$value
        fp.top = integrate(f2,-Inf,marker[i]/beta)$value
      }
      tp.low = integrate(f1,-Inf,Inf)$value
      fp.low = integrate(f2,-Inf,Inf)$value
      TPR[i] = tp.top/tp.low
      FPR[i] = fp.top/fp.low
    }
    TPR = c(1,TPR,0) ;FPR = c(1,FPR,0) 
    f = approxfun(FPR,TPR)
    AUC[k] = integrate(f,0,1)$value
  }
  AUC.weight = AUC*2*(1-pnorm((log(t0)-mu)/sigma))/(1-(1-pnorm((log(tmax)-mu)/sigma))^2)
  Cindex = mean(AUC.weight)
  return(out = list(Cindex = Cindex))
}
## loglogistic
MyrisksetAUC.loglogistic.True <- function(mu,sigma,beta, mu.z=0,tau,method = "Cox"){
  u = runif(1000)
  z <- rnorm(1000,mu.z,1)
  if(method=="Cox"){
    t0 = exp(sigma*log(exp(-log(u)*exp(-beta*z))-1)+mu)
  }
  else if(method == "AFT"){
    t0 = exp(sigma*log((1/u)-1) - beta*z  + mu)
  }
  else if(method == "PO"){
    t0 = exp(sigma*(log((1/u)-1)-beta*z) + mu) 
  }
  tmax = quantile(t0,tau)
  t0 = t0[which(t0 <= tmax)]
  t0=sort(t0)
  AUC = rep(0,length(t0))
  for(k in 1:length(t0)){
    if(method=="Cox"){
      f1 = function(x){exp(-0.5*(x-mu.z)^2+beta*x)*(1+exp((log(t0[k])-mu)/sigma))^(-exp(beta*x))}
      f2 = function(x){exp(-0.5*(x-mu.z)^2)*(1+exp((log(t0[k])-mu)/sigma))^(-exp(beta*x))}
    }
    else if(method == "AFT"){
      f1 = function(x){((1+exp((log(t0[k])-mu)/sigma))^(-exp(beta*x)))*exp(beta*x-0.5*(x-mu.z)^2)}
      f2 = function(x){exp(-0.5*(x-mu.z)^2)*(exp(-exp(beta*x)*(log(1+exp((log(t0[k])-mu)/sigma)))))}
    }
    else if(method == "PO"){
      f1 = function(x){exp(-0.5*(x-mu.z)^2 + beta*x)/((1+(exp((log(t0[k])-mu)/sigma))*exp(beta*x))^2)}
      f2 = function(x){exp(-0.5*(x-mu.z)^2)/(1+(exp((log(t0[k])-mu)/sigma))*exp(beta*x))}
    }
    c = rnorm(1000,0,1)
    FPR = rep(0,1000);TPR = rep(0,1000)
    marker = sort(c)
    for(i in 1:length(c)){
      if(beta>0){
        tp.top = integrate(f1,marker[i]/beta,Inf)$value
        fp.top = integrate(f2,marker[i]/beta,Inf)$value
      }else{
        tp.top = integrate(f1,-Inf,marker[i]/beta)$value
        fp.top = integrate(f2,-Inf,marker[i]/beta)$value
      }
      tp.low = integrate(f1,-Inf,Inf)$value
      fp.low = integrate(f2,-Inf,Inf)$value
      TPR[i] = tp.top/tp.low
      FPR[i] = fp.top/fp.low
    }
    TPR = c(1,TPR,0) ;FPR = c(1,FPR,0) 
    f = approxfun(FPR,TPR)
    AUC[k] = integrate(f,0,1)$value
  }
  AUC.weight = AUC*2*((1+exp((log(t0)-mu)/sigma))^-1)/(1-(1/(1+exp((log(tmax)-mu)/sigma)))^2)
  Cindex = mean(AUC.weight)
  return(out = list(Cindex = Cindex))
}
