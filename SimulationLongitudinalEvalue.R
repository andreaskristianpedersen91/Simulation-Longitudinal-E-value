library(forecast)
library(dplyr)
##Normal
set.seed(04022025)
Simulation<-function(N,rep){
V<-replicate(rep,expr = {
  time<-as.integer(runif(1000000, min=1, max=N))
  Sim<-as.data.frame(time)
  Sim$year<-floor(time/120)
  Sim$month<-Sim$time-120*Sim$year
  Sim$intervention<-0
  Sim$intervention[time>N/2]<-1
  Sim$U<-rnorm(1000000, 1+2*Sim$intervention, 1)
  Sim$outcome<-rnorm(1000000, 5+5*Sim$U, 10)
  Sim$random<-0
  for (j in 0:120){
    Sim$random[Sim$month==j]<-rnorm(1,0,2)
  }
  Sim$outcome<-Sim$outcome+Sim$random
  TimeseriesSim<-aggregate(Sim, by=list(time), FUN="mean")
  XY<-lm(outcome~intervention,data=TimeseriesSim)
  UY<-lm(outcome~U, data=TimeseriesSim)
  UX<-glm(U~intervention, data = TimeseriesSim)
  mean<-c(XY$coefficients[2],UY$coefficients[2],(UX$coefficients[2]+UX$coefficients[1])/(UX$coefficients[1]))
  Prediction<-as.vector(TimeseriesSim$outcome[TimeseriesSim$time<31])
  EFmin<-min(abs(as.vector(TimeseriesSim$outcome[TimeseriesSim$time>29])-Prediction))
  c(XY$coefficients[2],UY$coefficients[2],(UX$coefficients[2]+UX$coefficients[1])/(UX$coefficients[1]), EFmin)
})
Vnew<-as.data.frame(t(V))
names(Vnew)<-c("mean EF_{XY}^{obs}","EF_{UY}", "EF_{UX}", "min EF_{XY}^{obs}")
Vnew$kontrolmean<-Vnew$`mean EF_{XY}^{obs}`/(Vnew$`EF_{UX}`-1)
Vnew$kontrolmin<-Vnew$`min EF_{XY}^{obs}`/(Vnew$`EF_{UX}`-1)
Vnew$biasmean<-Vnew$kontrolmean-Vnew$`EF_{UY}`
Vnew$biasmin<-Vnew$kontrolmin-Vnew$`EF_{UY}`
Aggregated<-c(mean(Vnew$`mean EF_{XY}^{obs}`), mean(Vnew$`min EF_{XY}^{obs}`), mean(Vnew$`EF_{UY}`),mean(Vnew$`EF_{UX}`),mean(Vnew$biasmean), sd(Vnew$biasmean)/sqrt(length(Vnew$biasmean)),mean(Vnew$biasmin),sd(Vnew$biasmin)/sqrt(length(Vnew$biasmin)))
V<-replicate(rep,expr = {
  outcome<-as.vector(arima.sim(model=list(order=c(1,1,1), ar=0.7, ma=3), n=N-1, sd=1))
  time<-seq(1,N)
  Sim<-as.data.frame(cbind(time,outcome))
  Sim$intervention<-0
  Sim$intervention[time>N/2]<-1
  Sim$U<-rnorm(N, 1+2*Sim$intervention, 1)
  Sim$outcome<-Sim$outcome+Sim$U
  XY<-auto.arima(as.vector(Sim$outcome),d=1,D=0, max.p = 1,max.q = 1, seasonal = FALSE,xreg = as.vector(Sim$intervention))
  XYmin<-auto.arima(as.vector(Sim[Sim$time<N/2,]$outcome),d=1,D=0, max.p = 1,max.q = 1, seasonal = FALSE)
  UY<-auto.arima(as.vector(Sim$outcome),d=1,D=0, max.p = 1,max.q = 1, seasonal = FALSE,xreg = as.vector(Sim$U))
  V<-as.data.frame(forecast(XYmin, h=N/2))
  V$time<-seq(N/2+1,N)
  Simfinal<-merge(V,Sim, by="time",all=TRUE)
  Simfinal$EF<-abs(Simfinal$outcome-Simfinal$`Point Forecast`)
  UX<-glm(U~intervention, data = Simfinal)
  c(XY[["coef"]][["xreg"]], UY[["coef"]][["xreg"]],(UX$coefficients[2]+UX$coefficients[1])/(UX$coefficients[1]),min(Simfinal$EF, na.rm = TRUE))
})
Vnew<-as.data.frame(t(V))
names(Vnew)<-c("mean EF_{XY}^{obs}","EF_{UY}", "EF_{UX}", "min EF_{XY}^{obs}")
 
Vnew$kontrolmean<-Vnew$`mean EF_{XY}^{obs}`/(Vnew$`EF_{UX}`-1)
Vnew$kontrolmin<-Vnew$`min EF_{XY}^{obs}`/(Vnew$`EF_{UX}`-1)
Vnew$biasmean<-Vnew$kontrolmean-Vnew$`EF_{UY}`
Vnew$biasmin<-Vnew$kontrolmin-Vnew$`EF_{UY}`
ARIMA<-c(mean(Vnew$`mean EF_{XY}^{obs}`), mean(Vnew$`min EF_{XY}^{obs}`), mean(Vnew$`EF_{UY}`),mean(Vnew$`EF_{UX}`),mean(Vnew$biasmean), sd(Vnew$biasmean)/sqrt(length(Vnew$biasmean)),mean(Vnew$biasmin),sd(Vnew$biasmin)/sqrt(length(Vnew$biasmin)))
 
Simulationtablenormal<-as.data.frame(rbind(Aggregated,ARIMA))
names(Simulationtablenormal)<-c("mean EF_{XY}^{obs}","min EF_{XY}^{obs}","EF_{UY}", "EF_{UX}","Bias mean approach","SE Bias mean approach","Bias min approach", "SE Bias min approach")
##Poisson

V<-replicate(rep,expr = {
  time<-as.integer(runif(1000000, min=1, max=N))
  Sim<-as.data.frame(time)
  Sim$year<-floor(time/12)
  Sim$month<-Sim$time-12*Sim$year
  Sim$intervention<-0
  Sim$intervention[time>30]<-1
  Sim$U<-rpois(1000000, 1+2*Sim$intervention)
  Sim$outcome<-rnorm(1000000, 5+5*Sim$U, 10)
  Sim$random<-0
  for (j in 0:12){
    Sim$random[Sim$month==j]<-rnorm(1,0,2)
  }
  Sim$outcome<-Sim$outcome+Sim$random
  TimeseriesSim<-aggregate(Sim, by=list(time), FUN="mean")
  XY<-lm(outcome~intervention,data=TimeseriesSim)
  UY<-lm(outcome~U, data=TimeseriesSim)
  UX<-glm(U~intervention, data = TimeseriesSim, family = poisson)
  Prediction<-as.vector(TimeseriesSim$outcome[TimeseriesSim$time<N/2+1])
  EFmin<-min(abs(as.vector(TimeseriesSim$outcome[TimeseriesSim$time>N/2-1])-Prediction))
  c(XY$coefficients[2],UY$coefficients[2],exp(UX$coefficients[2]), EFmin)
})
Vnew<-as.data.frame(t(V))
names(Vnew)<-c("mean EF_{XY}^{obs}","EF_{UY}", "EF_{UX}", "min EF_{XY}^{obs}")
Vnew$kontrolmean<-Vnew$`mean EF_{XY}^{obs}`/(Vnew$`EF_{UX}`-1)
Vnew$kontrolmin<-Vnew$`min EF_{XY}^{obs}`/(Vnew$`EF_{UX}`-1)
Vnew$biasmean<-Vnew$kontrolmean-Vnew$`EF_{UY}`
Vnew$biasmin<-Vnew$kontrolmin-Vnew$`EF_{UY}`
Aggregated<-c(mean(Vnew$`mean EF_{XY}^{obs}`), mean(Vnew$`min EF_{XY}^{obs}`), mean(Vnew$`EF_{UY}`),mean(Vnew$`EF_{UX}`),mean(Vnew$biasmean), sd(Vnew$biasmean)/sqrt(length(Vnew$biasmean)),mean(Vnew$biasmin),sd(Vnew$biasmin)/sqrt(length(Vnew$biasmin)))

V<-replicate(rep,expr = {
  outcome<-as.vector(arima.sim(model=list(order=c(1,1,1), ar=0.7, ma=3), n=N-1, sd=1))
  time<-seq(1,N)
  Sim<-as.data.frame(cbind(time,outcome))
  Sim$intervention<-0
  Sim$intervention[time>N/2]<-1
  Sim$U<-rpois(N, 1+2*Sim$intervention)
  Sim$outcome<-Sim$outcome+Sim$U
  XY<-auto.arima(as.vector(Sim$outcome),d=1,D=0, max.p = 1,max.q = 1, seasonal = FALSE,xreg = as.vector(Sim$intervention))
  XYmin<-auto.arima(as.vector(Sim[Sim$time<N/2,]$outcome),d=1,D=0, max.p = 1,max.q = 1, seasonal = FALSE)
  UY<-auto.arima(as.vector(Sim$outcome),d=1,D=0, max.p = 1,max.q = 1, seasonal = FALSE,xreg = as.vector(Sim$U))
  V<-as.data.frame(forecast(XYmin, h=N/2))
  V$time<-seq(N/2+1,N)
  Simfinal<-merge(V,Sim, by="time",all=TRUE)
  Simfinal$EF<-abs(Simfinal$outcome-Simfinal$`Point Forecast`)
  UX<-glm(U~intervention, data = Simfinal, family = poisson)
  c(XY[["coef"]][["xreg"]], UY[["coef"]][["xreg"]],exp(UX$coefficients[2]),min(Simfinal$EF, na.rm = TRUE))
})
Vnew<-as.data.frame(t(V))
names(Vnew)<-c("mean EF_{XY}^{obs}","EF_{UY}", "EF_{UX}", "min EF_{XY}^{obs}")

Vnew$kontrolmean<-Vnew$`mean EF_{XY}^{obs}`/(Vnew$`EF_{UX}`-1)
Vnew$kontrolmin<-Vnew$`min EF_{XY}^{obs}`/(Vnew$`EF_{UX}`-1)
Vnew$biasmean<-Vnew$kontrolmean-Vnew$`EF_{UY}`
Vnew$biasmin<-Vnew$kontrolmin-Vnew$`EF_{UY}`
ARIMA<-c(mean(Vnew$`mean EF_{XY}^{obs}`), mean(Vnew$`min EF_{XY}^{obs}`), mean(Vnew$`EF_{UY}`),mean(Vnew$`EF_{UX}`),mean(Vnew$biasmean), sd(Vnew$biasmean)/sqrt(length(Vnew$biasmean)),mean(Vnew$biasmin),sd(Vnew$biasmin)/sqrt(length(Vnew$biasmin)))

Simulationtablepoisson<-as.data.frame(rbind(Aggregated,ARIMA))
names(Simulationtablepoisson)<-c("mean EF_{XY}^{obs}","min EF_{XY}^{obs}","EF_{UY}", "EF_{UX}","Bias mean approach","SE Bias mean approach","Bias min approach", "SE Bias min approach")

##Binomial

V<-replicate(rep,expr = {
  U<-rbinom(1000000,size=1, prob=0.50)
  outcome<-rnorm(1000000, 0, 10)
  outcome<-outcome+5*U
  intervention<-rbinom(n=1000000, size=1, prob = 0.25+0.5*U)
  Sim<-as.data.frame(cbind(outcome,U, intervention))
  Sim$time<-0
  Sim[Sim$intervention==0,]$time<-as.integer(runif(nrow(Sim[Sim$intervention==0,]), min=1, max=N/2))
  Sim[Sim$intervention==1,]$time<-as.integer(runif(nrow(Sim[Sim$intervention==1,]), min=N/2+1, max=N))
  Sim$year<-floor(Sim$time/12)
  Sim$month<-Sim$time-12*Sim$year
  lm(outcome~intervention,data=Sim)
  TimeseriesSim<-Sim %>% group_by(time) %>% summarise(outcome=mean(outcome), intervention=mean(intervention), U=mean(U)) %>% as.data.frame()
  XY<-lm(outcome~intervention,data=TimeseriesSim)
  UY<-lm(outcome~U, data=TimeseriesSim)
  UX<-glm(U~intervention, data = TimeseriesSim, family = binomial(link = "logit"))
  Prediction<-as.vector(TimeseriesSim$outcome[TimeseriesSim$time<N/2+1])
  EFmin<-min(abs(as.vector(TimeseriesSim$outcome[TimeseriesSim$time>N/2-1])-Prediction))
  p1<-exp(UX$coefficients[2]+UX$coefficients[1])/(exp(UX$coefficients[2]+UX$coefficients[1])+1)
  p0<-exp(UX$coefficients[1])/(exp(UX$coefficients[1])+1)
  c(XY$coefficients[2],UY$coefficients[2],p1/p0, EFmin)
})
Vnew<-as.data.frame(t(V))
names(Vnew)<-c("mean EF_{XY}^{obs}","EF_{UY}", "EF_{UX}", "min EF_{XY}^{obs}")
Vnew$kontrolmean<-Vnew$`mean EF_{XY}^{obs}`/(Vnew$`EF_{UX}`-1)
Vnew$kontrolmin<-Vnew$`min EF_{XY}^{obs}`/(Vnew$`EF_{UX}`-1)
Vnew$biasmean<-Vnew$kontrolmean-Vnew$`EF_{UY}`
Vnew$biasmin<-Vnew$kontrolmin-Vnew$`EF_{UY}`
Aggregated<-c(mean(Vnew$`mean EF_{XY}^{obs}`), mean(Vnew$`min EF_{XY}^{obs}`), mean(Vnew$`EF_{UY}`),mean(Vnew$`EF_{UX}`),mean(Vnew$biasmean), sd(Vnew$biasmean)/sqrt(length(Vnew$biasmean)),mean(Vnew$biasmin),sd(Vnew$biasmin)/sqrt(length(Vnew$biasmin)))

V<-replicate(rep,expr = {
  outcome<-as.vector(arima.sim(model=list(order=c(1,1,1), ar=0.7, ma=3), n=N-1, sd=1))
  time<-seq(1,N)
  Sim<-as.data.frame(cbind(time,outcome))
  Sim$intervention<-0
  Sim$intervention[time>N/2]<-1
  Sim$U<-rbinom(n=N, size=1,prob=0.25+0.50*Sim$intervention)
  Sim$outcome<-Sim$outcome+5*Sim$U
  XY<-auto.arima(as.vector(Sim$outcome),d=1,D=0, max.p = 1,max.q = 1, seasonal = FALSE,xreg = as.vector(Sim$intervention))
  XYmin<-auto.arima(as.vector(Sim[Sim$time<N/2,]$outcome),d=1,D=0, max.p = 1,max.q = 1, seasonal = FALSE)
  UY<-auto.arima(as.vector(Sim$outcome),d=1,D=0, max.p = 1,max.q = 1, seasonal = FALSE,xreg = as.vector(Sim$U))
  V<-as.data.frame(forecast(XYmin, h=N/2))
  V$time<-seq(N/2+1,N)
  Simfinal<-merge(V,Sim, by="time",all=TRUE)
  Simfinal$EF<-abs(Simfinal$outcome-Simfinal$`Point Forecast`)
  UX<-glm(U~intervention, data = Simfinal, family = binomial(link = "logit"))
  p1<-exp(UX$coefficients[2]+UX$coefficients[1])/(exp(UX$coefficients[2]+UX$coefficients[1])+1)
  p0<-exp(UX$coefficients[1])/(exp(UX$coefficients[1])+1)
  c(XY[["coef"]][["xreg"]], UY[["coef"]][["xreg"]],p1/p0,min(Simfinal$EF, na.rm = TRUE))
})

Vnew<-as.data.frame(t(V))
names(Vnew)<-c("mean EF_{XY}^{obs}","EF_{UY}", "EF_{UX}", "min EF_{XY}^{obs}")

Vnew$kontrolmean<-Vnew$`mean EF_{XY}^{obs}`/(Vnew$`EF_{UX}`-1)
Vnew$kontrolmin<-Vnew$`min EF_{XY}^{obs}`/(Vnew$`EF_{UX}`-1)
Vnew$biasmean<-Vnew$kontrolmean-Vnew$`EF_{UY}`
Vnew$biasmin<-Vnew$kontrolmin-Vnew$`EF_{UY}`
ARIMA<-c(mean(Vnew$`mean EF_{XY}^{obs}`), mean(Vnew$`min EF_{XY}^{obs}`), mean(Vnew$`EF_{UY}`),mean(Vnew$`EF_{UX}`),mean(Vnew$biasmean), sd(Vnew$biasmean)/sqrt(length(Vnew$biasmean)),mean(Vnew$biasmin),sd(Vnew$biasmin)/sqrt(length(Vnew$biasmin)))

Simulationtablebinomial<-as.data.frame(rbind(Aggregated,ARIMA))
names(Simulationtablebinomial)<-c("mean EF_{XY}^{obs}","min EF_{XY}^{obs}","EF_{UY}", "EF_{UX}","Bias mean approach","SE Bias mean approach","Bias min approach", "SE Bias min approach")
Simulationtable$N<-N
names(Simulationtablebinomial)<-c("mean EF_{XY}^{obs}","min EF_{XY}^{obs}","EF_{UY}", "EF_{UX}","Bias mean approach","SE Bias mean approach","Bias min approach", "SE Bias min approach")

Simulationtable<-rbind(Simulationtablenormal,Simulationtablepoisson,Simulationtablebinomial)
Simulationtable$N<-N
Simulationtable$CIllmeanbias<-Simulationtable$`Bias mean approach`-Simulationtable$`SE Bias mean approach`*1.96
Simulationtable$CIulmeanbias<-Simulationtable$`Bias mean approach`+Simulationtable$`SE Bias mean approach`*1.96
Simulationtable$`Bias mean approach 95CI`<-paste(round(Simulationtable$`Bias mean approach`, digits = 3),"(", round(Simulationtable$CIllmeanbias, digits = 3), ";", round(Simulationtable$CIulmeanbias, digits = 3), ")")

Simulationtable$CIllminbias<-Simulationtable$`Bias min approach`-Simulationtable$`SE Bias min approach`*1.96
Simulationtable$CIulminbias<-Simulationtable$`Bias min approach`+Simulationtable$`SE Bias min approach`*1.96
Simulationtable$`Bias min approach 95CI`<-paste(round(Simulationtable$`Bias min approach`, digits = 3),"(", round(Simulationtable$CIllminbias, digits = 3), ";", round(Simulationtable$CIulminbias, digits = 3), ")")

Simulationtable<-as.data.frame(subset(Simulationtable,select=c("N","mean EF_{XY}^{obs}","min EF_{XY}^{obs}","EF_{UY}","EF_{UX}", "Bias mean approach 95CI", "Bias min approach 95CI")))
Simulationtable$`mean EF_{XY}^{obs}`<-round(Simulationtable$`mean EF_{XY}^{obs}`, digits=3)
Simulationtable$`min EF_{XY}^{obs}`<-round(Simulationtable$`min EF_{XY}^{obs}`, digits=3)
Simulationtable$`EF_{UY}`<-round(Simulationtable$`EF_{UY}`, digits=3)
Simulationtable$`EF_{UX}`<-round(Simulationtable$`EF_{UX}`, digits=3)
Simulationtable
}


Simulation100<-Simulation(100,1000)
Simulation500<-Simulation(500,1000)
Simulation1000<-Simulation(1000,1000)
kableExtra::kable(rbind(Simulation100,Simulation500,Simulation1000), format = "latex")