library(LongitudinalEvalue)
library(yuima)
library(dplyr)
library(ggplot2)
library(kableExtra)
library(parallel)
library(future.apply)
plan(multisession, workers = 10)




Simulation<-function(EFXUsim, EFUYsim, volatilitysim, rateEFXUsim=1,mu1,sigma1,mu2=0,sigma2=0, mu3=0,sigma3=0,mu4=0,sigma4=0,Ito){
SimulatingconfounderIto1<-function(mu,sigma){
mod <- setModel(drift = mu, diffusion = sigma, state.variable = "Y", solve.variable = "Y", time.variable = "t", xinit = 1)
samp <- setSampling(Terminal = 100, n = 100)
Simulation <- simulate(mod, sampling = samp, true.parameter = c(mu,sigma))
Simulation<-as.data.frame(cbind(Simulation@data@original.data,seq(1,101)))
names(Simulation)<-c("Confounder", "time")
Simulation
}

SimulatingconfounderIto2<-function(mu1,sigma1,mu2,sigma2){
  mod <- setModel(drift = mu1, diffusion = sigma1, state.variable = "Y", solve.variable = "Y", time.variable = "t", xinit = 1)
  samp <- setSampling(Terminal = 50, n = 50)
  Simulation <- simulate(mod, sampling = samp, true.parameter = c(mu1,sigma1))
  Simulation1<-as.data.frame(cbind(Simulation@data@original.data,seq(0,50)))
  names(Simulation1)<-c("Confounder", "time")

  mod <- setModel(drift = mu2, diffusion = sigma2, state.variable = "Y", solve.variable = "Y", time.variable = "t", xinit = Simulation1[Simulation1$time==50,]$Confounder)
  samp <- setSampling(Terminal = 100, n = 50)
  Simulation <- simulate(mod, sampling = samp, true.parameter = c(mu2,sigma2))
  Simulation2<-as.data.frame(cbind(Simulation@data@original.data,seq(51,101)))
  names(Simulation2)<-c("Confounder", "time")
  rbind(Simulation1,Simulation2)
}



SimulatingconfounderIto3<-function(mu1,sigma1,mu2,sigma2, mu3,sigma3){
 mod <- setModel(drift = mu1, diffusion = sigma1, state.variable = "Y", solve.variable = "Y", time.variable = "t", xinit = 1)
  samp <- setSampling(Terminal = 33, n = 33)
  Simulation <- simulate(mod, sampling = samp, true.parameter = c(mu1,sigma1))
  Simulation1<-as.data.frame(cbind(Simulation@data@original.data,seq(0,33)))
  names(Simulation1)<-c("Confounder", "time")

  mod <- setModel(drift = mu2, diffusion = sigma2, state.variable = "Y", solve.variable = "Y", time.variable = "t", xinit = Simulation1[Simulation1$time==33,]$Confounder)
  samp <- setSampling(Terminal = 33, n = 33)
  Simulation <- simulate(mod, sampling = samp, true.parameter = c(mu2,sigma2))
  Simulation2<-as.data.frame(cbind(Simulation@data@original.data,seq(34,67)))
  names(Simulation2)<-c("Confounder", "time")

  mod <- setModel(drift = mu3, diffusion = sigma3, state.variable = "Y", solve.variable = "Y", time.variable = "t", xinit = Simulation2[Simulation2$time==67,]$Confounder)
  samp <- setSampling(Terminal = 32, n = 32)
  Simulation <- simulate(mod, sampling = samp, true.parameter = c(mu3,sigma3))
  Simulation3<-as.data.frame(cbind(Simulation@data@original.data,seq(68,100)))
  names(Simulation3)<-c("Confounder", "time")

  rbind(Simulation1,Simulation2,Simulation3)
}



SimulatingconfounderIto4<-function(mu1,sigma1,mu2,sigma2, mu3,sigma3,mu4,sigma4){
  mod <- setModel(drift = mu1, diffusion = sigma1, state.variable = "Y", solve.variable = "Y", time.variable = "t", xinit = 1)
  samp <- setSampling(Terminal = 25, n = 25)
  Simulation <- simulate(mod, sampling = samp, true.parameter = c(mu1,sigma1))
  Simulation1<-as.data.frame(cbind(Simulation@data@original.data,seq(0,25)))
  names(Simulation1)<-c("Confounder", "time")

  mod <- setModel(drift = mu2, diffusion = sigma2, state.variable = "Y", solve.variable = "Y", time.variable = "t", xinit = Simulation1[Simulation1$time==25,]$Confounder)
  samp <- setSampling(Terminal = 25, n = 25)
  Simulation <- simulate(mod, sampling = samp, true.parameter = c(mu2,sigma2))
  Simulation2<-as.data.frame(cbind(Simulation@data@original.data,seq(26,51)))
  names(Simulation2)<-c("Confounder", "time")

  mod <- setModel(drift = mu3, diffusion = sigma3, state.variable = "Y", solve.variable = "Y", time.variable = "t", xinit = Simulation2[Simulation2$time==51,]$Confounder)
  samp <- setSampling(Terminal = 25, n = 25)
  Simulation <- simulate(mod, sampling = samp, true.parameter = c(mu3,sigma3))
  Simulation3<-as.data.frame(cbind(Simulation@data@original.data,seq(52,77)))
  names(Simulation3)<-c("Confounder", "time")

  mod <- setModel(drift = mu4, diffusion = sigma4, state.variable = "Y", solve.variable = "Y", time.variable = "t", xinit = Simulation3[Simulation3$time==77,]$Confounder)
  samp <- setSampling(Terminal = 22, n = 22)
  Simulation <- simulate(mod, sampling = samp, true.parameter = c(mu4,sigma4))
  Simulation4<-as.data.frame(cbind(Simulation@data@original.data,seq(78,100)))
  names(Simulation4)<-c("Confounder", "time")

  rbind(Simulation1,Simulation2,Simulation3,Simulation4)
}



if (Ito==1) {
  Confounder<-SimulatingconfounderIto1(mu1,sigma1)
  mu_t<-mu1
  sigma_t<-sigma1
  minmu_t<-mu1-0.1
  minsigma_t<-sigma1-0.1
  maxmu_t<-mu1+0.1
  maxsigma_t<-sigma1+0.1
} else if (Ito==2) {
  Confounder<-SimulatingconfounderIto2(mu1,sigma1,mu2,sigma2)
  mu_t<-mean(c(mu1,mu2))
  sigma_t<-mean(c(sigma1,sigma2))
  minmu_t<-min(c(mu1,mu2))
  minsigma_t<-min(c(sigma1,sigma2))
  maxmu_t<-max(c(mu1,mu2))
  maxsigma_t<-max(c(sigma1,sigma2))
} else if (Ito==3) {
  Confounder<-SimulatingconfounderIto3(mu1,sigma1,mu2,sigma2,mu3,sigma3)
  mu_t<-mean(c(mu1,mu2,mu3))
  sigma_t<-mean(c(sigma1,sigma2,sigma3))
  minmu_t<-min(c(mu1,mu2,mu3))
  minsigma_t<-min(c(sigma1,sigma2,sigma3))
  maxmu_t<-max(c(mu1,mu2,mu3))
  maxsigma_t<-max(c(sigma1,sigma2,sigma3))
} else {
  Confounder<-SimulatingconfounderIto4(mu1,sigma1,mu2,sigma2,mu3,sigma3,mu4,sigma4)
  mu_t<-mean(c(mu1,mu2,mu3,mu4))
  sigma_t<-mean(c(sigma1,sigma2,sigma3,sigma4))
  minmu_t<-min(c(mu1,mu2,mu3,mu4))
  minsigma_t<-min(c(sigma1,sigma2,sigma3,sigma4))
  maxmu_t<-max(c(mu1,mu2,mu3,mu4))
  maxsigma_t<-max(c(sigma1,sigma2,sigma3,sigma4))

}
Simulation<-Confounder
Simulation$U_t0<-Simulation$Confounder
Simulation$U_t1<-EFXUsim*Simulation$Confounder

Simulation$Y0<-Simulation$U_t0+rnorm(nrow(Simulation),0,1)
Simulation$Y1<-EFUYsim*Simulation$U_t1+rnorm(nrow(Simulation),0,1)
Simulation$ATEobserved<-Simulation$Y1-Simulation$Y0
Simulation$time<-1:nrow(Simulation)

drift <-"adiff"
diffusion <-" sigmadiff"
Delta<-1/228
mod <-setModel(drift=drift , diffusion =diffusion , state.var="ATEobserved",time.var = "time", solve.var = "ATEobserved",xinit = 0.05)
model <-setYuima(model=mod , data= setData(zoo(Simulation$ATEobserved, order.by
                                               = Simulation$time),delta=Delta))
fit <-qmle(model , start=list(adiff =mu_t,sigmadiff =sigma_t) ,lower=list(adiff =minmu_t,sigmadiff =minsigma_t) ,upper=list(adiff =maxmu_t,sigmadiff =maxsigma_t))

U_tpara<-SDEconfounderYuima1d(fit , ATE=0, EFXU =EFXUsim , EFUY =EFXUsim , volatility =0,
                     rateEFXU =0,initial.value = 1)
Paradata<-as.data.frame(U_tpara@data)
names(Paradata)<-c("time","meanConfounderpara", "CIllpara","CIulpara")

names(Results_final)<-c("time", "ATE", "Confounder_true", "Y0", "Y1", "Confounder_para","CIll_para", "CIul_para")

Results_final$RMSEAparaperc<-((Results_final$Confounder_para-Results_final$Confounder_true)/Results_final$Confounder_true)^2
Results_final$RMSEApara<-((Results_final$Confounder_para-Results_final$Confounder_true))^2
Results_final$Abovepara<-ifelse(Results_final$Confounder_true>Results_final$CIll_para,1,0)
Results_final$Belowpara<-ifelse(Results_final$Confounder_true<Results_final$CIul_para,1,0)

Results_final$withinpara<-min(Results_final$Abovepara,Results_final$Belowpara)

Simulation_results<-Results_final %>% subset(select=c("RMSEApara","RMSEAparaperc","withinpara")) %>% summarise(RMSEApara=sqrt(mean(RMSEApara,na.rm=TRUE)), RMSEAparaperc=sqrt(mean(RMSEAparaperc,na.rm=TRUE)), Coveragepara=mean(withinpara,na.rm=TRUE)) %>% as.vector()
c(Ito, EFXUsim, EFUYsim, volatilitysim,Simulation_results)
}



# Total number of replications
n_reps <- 1000
# Progress checkpoints

library(future.apply)

library(progressr)

V<-future_replicate(n_reps,expr = {

  
  # Define parameter grid manually
  param_grid <- expand.grid(
    Ito=1:4,
    EFUYsim = 1:3,
    EFXUsim = 2:4,
    volatilitysim = c(0.1,1,2)
  )
  
  # Run simulations over the grid
  results <- unlist(
    apply(param_grid, 1, function(row) {
      Simulation(
        EFXUsim = as.numeric(row["EFXUsim"]),
        EFUYsim = as.numeric(row["EFUYsim"]),
        volatilitysim = as.numeric(row["volatilitysim"]),
        mu1 = 0.1, sigma1 = 0.1,
        mu2 = 0.2, sigma2 = 0.2,
        mu3 = 0.3, sigma3 = 0.3,
        mu4 = 0.4, sigma4 = 0.4,
        Ito = as.integer(row["Ito"])
      )
      
      
    })
  )
  
  results
  
  
})


Vnew<-as.data.frame(t(V))
View(Vnew)

names(Vnew) <- make.unique(names(Vnew))


TableIto <- Vnew %>%
  summarise(across(everything(), mean, na.rm = TRUE))
View(TableIto)
TableIto <- as.data.frame(matrix(TableIto, ncol = 10, byrow = TRUE))
names(TableIto)<-c("Ito", "EFXU", "EFUY", "volatility", "RMSEA parametric", "RMSEA non parametric","RMSEA parametric (%)", "RMSEA non parametric (%)", "Coverage parametric", "Coverage non parametric")

TableIto1<-TableIto %>% filter(Ito==1)
TableIto2<-TableIto %>% filter(Ito==2)
TableIto3<-TableIto %>% filter(Ito==3)
TableIto4<-TableIto %>% filter(Ito==4)




# Example data frames (2 rows each, 9 columns)

# Add an ID column for grouping
TableIto1$group <- "Ito 1"
TableIto2$group <- "Ito 2"
TableIto3$group <- "Ito 3"
TableIto4$group <- "Ito 4"

# Combine vertically
combined_TableIto <- bind_rows(TableIto1, TableIto2, TableIto3, TableIto4)

# Remove the group column from display, but keep for grouping
combined_TableIto_display <- combined_TableIto %>% select(-group)

# Now make the table and add group headers by row index
n_rows <- nrow(TableIto1)  # r
