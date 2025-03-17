library(forecast)
library(dplyr)
##Normal
set.seed(04022025)
V<-replicate(1000,expr = {
  time<-as.integer(runif(1000000, min=1, max=600))
  Sim<-as.data.frame(time)
  Sim$year<-floor(time/120)
  Sim$month<-Sim$time-120*Sim$year
  Sim$intervention<-0
  Sim$intervention[time>300]<-1
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
Aggregated<-c(mean(Vnew$`mean EF_{XY}^{obs}`), mean(Vnew$`min EF_{XY}^{obs}`), mean(Vnew$`EF_{UY}`),mean(Vnew$`EF_{UX}`),mean(Vnew$biasmean),mean(Vnew$biasmin))

V<-replicate(1000,expr = {
  outcome<-as.vector(arima.sim(model=list(order=c(1,1,1), ar=0.7, ma=3), n=599, sd=1))
  time<-seq(1,600)
  Sim<-as.data.frame(cbind(time,outcome))
  Sim$intervention<-0
  Sim$intervention[time>300]<-1
  Sim$U<-rnorm(600, 1+2*Sim$intervention, 1)
  Sim$outcome<-Sim$outcome+Sim$U
  XY<-auto.arima(as.vector(Sim$outcome),d=1,D=0, max.p = 1,max.q = 1, seasonal = FALSE,xreg = as.vector(Sim$intervention))
  XYmin<-auto.arima(as.vector(Sim[Sim$time<300,]$outcome),d=1,D=0, max.p = 1,max.q = 1, seasonal = FALSE)
  UY<-auto.arima(as.vector(Sim$outcome),d=1,D=0, max.p = 1,max.q = 1, seasonal = FALSE,xreg = as.vector(Sim$U))
  V<-as.data.frame(forecast(XYmin, h=300))
  V$time<-seq(301,600)
  Simfinal<-merge(V,Sim, by="time",all=TRUE)
  Simfinal$EF<-abs(Simfinal$outcome-Simfinal$`Point Forecast`)
  UX<-glm(U~intervention, data = Simfinal)
  c(XY[["coef"]][["xreg"]], UY[["coef"]][["xreg"]],(UX$coefficients[2]+UX$coefficients[1])/(UX$coefficients[1]),min(Simfinal$EF, na.rm = TRUE))
})
Vnew<-as.data.frame(t(V))
View(Vnew)
names(Vnew)<-c("mean EF_{XY}^{obs}","EF_{UY}", "EF_{UX}", "min EF_{XY}^{obs}")

Vnew$kontrolmean<-Vnew$`mean EF_{XY}^{obs}`/(Vnew$`EF_{UX}`-1)
Vnew$kontrolmin<-Vnew$`min EF_{XY}^{obs}`/(Vnew$`EF_{UX}`-1)
Vnew$biasmean<-Vnew$kontrolmean-Vnew$`EF_{UY}`
Vnew$biasmin<-Vnew$kontrolmin-Vnew$`EF_{UY}`
ARIMA<-c(mean(Vnew$`mean EF_{XY}^{obs}`), mean(Vnew$`min EF_{XY}^{obs}`), mean(Vnew$`EF_{UY}`),mean(Vnew$`EF_{UX}`),mean(Vnew$biasmean),mean(Vnew$biasmin))

Simulationtablenormal<-as.data.frame(rbind(Aggregated,ARIMA))
names(Simulationtablenormal)<-c("mean EF_{XY}^{obs}","min EF_{XY}^{obs}","EF_{UY}", "EF_{UX}","Bias mean approach", "Bias min approach")

##Poisson

V<-replicate(1000,expr = {
  time<-as.integer(runif(1000000, min=1, max=600))
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
  Prediction<-as.vector(TimeseriesSim$outcome[TimeseriesSim$time<301])
  EFmin<-min(abs(as.vector(TimeseriesSim$outcome[TimeseriesSim$time>299])-Prediction))
  c(XY$coefficients[2],UY$coefficients[2],exp(UX$coefficients[2]), EFmin)
})
Vnew<-as.data.frame(t(V))
names(Vnew)<-c("mean EF_{XY}^{obs}","EF_{UY}", "EF_{UX}", "min EF_{XY}^{obs}")

Vnew$kontrolmean<-Vnew$`mean EF_{XY}^{obs}`/(Vnew$`EF_{UX}`-1)
Vnew$kontrolmin<-Vnew$`min EF_{XY}^{obs}`/(Vnew$`EF_{UX}`-1)
Vnew$biasmean<-Vnew$kontrolmean-Vnew$`EF_{UY}`
Vnew$biasmin<-Vnew$kontrolmin-Vnew$`EF_{UY}`
Aggregated<-c(mean(Vnew$`mean EF_{XY}^{obs}`), mean(Vnew$`min EF_{XY}^{obs}`), mean(Vnew$`EF_{UY}`),mean(Vnew$`EF_{UX}`),mean(Vnew$biasmean),mean(Vnew$biasmin))

V<-replicate(1000,expr = {
  outcome<-as.vector(arima.sim(model=list(order=c(1,1,1), ar=0.7, ma=3), n=599, sd=1))
  time<-seq(1,600)
  Sim<-as.data.frame(cbind(time,outcome))
  Sim$intervention<-0
  Sim$intervention[time>300]<-1
  Sim$U<-rpois(600, 1+2*Sim$intervention)
  Sim$outcome<-Sim$outcome+Sim$U
  XY<-auto.arima(as.vector(Sim$outcome),d=1,D=0, max.p = 1,max.q = 1, seasonal = FALSE,xreg = as.vector(Sim$intervention))
  XYmin<-auto.arima(as.vector(Sim[Sim$time<300,]$outcome),d=1,D=0, max.p = 1,max.q = 1, seasonal = FALSE)
  UY<-auto.arima(as.vector(Sim$outcome),d=1,D=0, max.p = 1,max.q = 1, seasonal = FALSE,xreg = as.vector(Sim$U))
  V<-as.data.frame(forecast(XYmin, h=300))
  V$time<-seq(301,600)
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
ARIMA<-c(mean(Vnew$`mean EF_{XY}^{obs}`), mean(Vnew$`min EF_{XY}^{obs}`), mean(Vnew$`EF_{UY}`),mean(Vnew$`EF_{UX}`),mean(Vnew$biasmean),mean(Vnew$biasmin))

Simulationtablepoisson<-as.data.frame(rbind(Aggregated,ARIMA))
names(Simulationtablepoisson)<-c("mean EF_{XY}^{obs}","min EF_{XY}^{obs}","EF_{UY}", "EF_{UX}","Bias mean approach", "Bias min approach")
View(Simulationtablepoisson)

##Binomial

V<-replicate(1000,expr = {
  U<-rbinom(1000000,size=1, prob=0.50)
  outcome<-rnorm(1000000, 0, 10)
  outcome<-outcome+5*U
  intervention<-rbinom(n=1000000, size=1, prob = 0.25+0.5*U)
  Sim<-as.data.frame(cbind(outcome,U, intervention))
  Sim$time<-0
  Sim[Sim$intervention==0,]$time<-as.integer(runif(nrow(Sim[Sim$intervention==0,]), min=1, max=300))
  Sim[Sim$intervention==1,]$time<-as.integer(runif(nrow(Sim[Sim$intervention==1,]), min=301, max=600))
  Sim$year<-floor(Sim$time/12)
  Sim$month<-Sim$time-12*Sim$year
  lm(outcome~intervention,data=Sim)
  TimeseriesSim<-Sim %>% group_by(time) %>% summarise(outcome=mean(outcome), intervention=mean(intervention), U=mean(U)) %>% as.data.frame()
  XY<-lm(outcome~intervention,data=TimeseriesSim)
  UY<-lm(outcome~U, data=TimeseriesSim)
  UX<-glm(U~intervention, data = TimeseriesSim, family = binomial(link = "logit"))
  Prediction<-as.vector(TimeseriesSim$outcome[TimeseriesSim$time<301])
  EFmin<-min(abs(as.vector(TimeseriesSim$outcome[TimeseriesSim$time>299])-Prediction))
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
Aggregated<-c(mean(Vnew$`mean EF_{XY}^{obs}`), mean(Vnew$`min EF_{XY}^{obs}`), mean(Vnew$`EF_{UY}`),mean(Vnew$`EF_{UX}`),mean(Vnew$biasmean),mean(Vnew$biasmin))

V<-replicate(1000,expr = {
  outcome<-as.vector(arima.sim(model=list(order=c(1,1,1), ar=0.7, ma=3), n=599, sd=1))
  time<-seq(1,600)
  Sim<-as.data.frame(cbind(time,outcome))
  Sim$intervention<-0
  Sim$intervention[time>300]<-1
  Sim$U<-rbinom(n=600, size=1,prob=0.25+0.50*Sim$intervention)
  Sim$outcome<-Sim$outcome+5*Sim$U
  XY<-auto.arima(as.vector(Sim$outcome),d=1,D=0, max.p = 1,max.q = 1, seasonal = FALSE,xreg = as.vector(Sim$intervention))
  XYmin<-auto.arima(as.vector(Sim[Sim$time<300,]$outcome),d=1,D=0, max.p = 1,max.q = 1, seasonal = FALSE)
  UY<-auto.arima(as.vector(Sim$outcome),d=1,D=0, max.p = 1,max.q = 1, seasonal = FALSE,xreg = as.vector(Sim$U))
  V<-as.data.frame(forecast(XYmin, h=300))
  V$time<-seq(301,600)
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
ARIMA<-c(mean(Vnew$`mean EF_{XY}^{obs}`), mean(Vnew$`min EF_{XY}^{obs}`), mean(Vnew$`EF_{UY}`),mean(Vnew$`EF_{UX}`),mean(Vnew$biasmean),mean(Vnew$biasmin))

Simulationtablebinomial<-as.data.frame(rbind(Aggregated,ARIMA))
names(Simulationtablebinomial)<-c("mean EF_{XY}^{obs}","min EF_{XY}^{obs}","EF_{UY}", "EF_{UX}","Bias mean approach", "Bias min approach")
View(Simulationtablebinomial)
Simulationtable<-rbind(Simulationtablenormal,Simulationtablepoisson,Simulationtablebinomial)

kableExtra::kable(Simulationtable,format = "latex")