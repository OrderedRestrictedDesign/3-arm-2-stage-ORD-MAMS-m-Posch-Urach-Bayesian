rm(list=ls())

####### Call for functions ####

source("functionsORDMAMSUPBAYESIAN_3arms2stages.R")

##### PARAMETERS for the functions #########################

theta2power <- 120 #theta of interest to power the design

interim <- c(0.5) #at which proportion of the first sample size we do the interim analysis 

ratio <- c(1,2) #allocation ratio sample size at different stages: n1=r[1], n2=r[2]*n1,...

rho <- c(1,1,1) #sample size in each treatment divided by sample size in the control: rho[1]=n11, rho[2]=n21, rho[3]=n01

name <- "triangular" #shape of the critical bound

stage <- 2 #number of stages

arms <- 3 #number of arms in the trial - including control arm

sigma <- 340 #known standard deviation

alpha <- 0.025 #alpha-level of the test

beta <- 0.2 #type II error

prec <- 4 #number of digits for the critical bounds

power <- "reject all" # type of power for the design: "reject all" or "reject at least one"


# Parameters and simulation scenario

seed <- 64736 #set seed to simulate the data

nsim <- 5000#10^4#10^6 #number of simulations

# Scenarios

trt1 <- c(0,0.5,0.5,0.5,0.5,0.5,0.5)

trt2 <- c(0,0,0.1,0.2,0.3,0.4,0.5)

trt1 <- round(trt1,2)

trt2 <- round(trt2,2)

# simulation scenarios for Section 6.3

scen <- cbind(1,1+trt1,1+trt2)

scenario <- matrix(data = scen,
                   nrow = nrow(scen),
                   ncol = ncol(scen))

scenario <- cbind(scenario, scenario[,2]-scenario[,1], scenario[,3]-scenario[,1])

colnames(scenario) <- c("mu0", "mu1", "mu2", "theta1", "theta2")

scenario <- as.data.frame(scenario)

# Simulation scenarios for Section 5 - case study

diff <- 120 # clinically relevant difference

FEV <- 2747 # common baseline mean FEV1: (2821*134 +2680*125 +2736*138)/(134+125+138)

plac <- rep(FEV, 7)

trt2 <- c(FEV-100,FEV-300,FEV-100,FEV,FEV,FEV+diff,FEV+diff)

trt1 <-  c(FEV-50,FEV,FEV,FEV,FEV+diff,FEV+diff, FEV)

scen <- cbind(plac,trt1,trt2)

scenario <- matrix(data = scen,
                   nrow = nrow(scen),
                   ncol = ncol(scen))

scenario <- cbind(scenario, scenario[,2]-scenario[,1], scenario[,3]-scenario[,1])

colnames(scenario) <- c("placebo",
                        "highdose",
                        "lowdose",
                        "thetaH",
                        "thetaL")

scenario <- as.data.frame(scenario)

# Bayesian design scenarios

grid <- c(0,20,40,60,80,120)

grid2 <- expand.grid(grid,grid,grid)+489

scenario <- unique(data.frame(grid2))

colnames(scenario) <- c("mu0s", "mu1s", "mu2s")

library(dplyr)

scenario <- scenario %>% dplyr::filter(
  mu0s ==489 & (mu1s==609 | mu2s == 609 | (mu1s==489 & mu2s == 489))
)

scenario <- cbind(scenario, scenario[,2]-scenario[,1], scenario[,3]-scenario[,1])

colnames(scenario) <- c("placebo",
                        "highdose",
                        "lowdose",
                        "thetaH",
                        "thetaL")

scenario <- as.data.frame(scenario)

##### 3 arms - 2 stage MAMS(m) design #################

func <- boundaries_3armJstagem(theta = rep(theta2power, arms-1), 
                               sigma = sigma,
                               interimN = stage,
                               alpha = alpha,
                               beta = beta,
                               r = ratio,
                               rho = rho,
                               prec = prec,
                               arms = arms,
                               ushape = name,
                               lshape = name)


popall <- func$`sample size per arm per stage`

# CALL the function to simulate the data

summary_AllProm <- simulAll(scenario, nsim, stage, 
                            popall, ratio, rho, sigma, seed= seed, func=func)

summary_AllProm

##### 3 arms - 2 stages Ordered Restricted Design ORD ############

ushape1 = ushape2 = name

fun <- boundariesSample_3armsORD_differentbounds(theta = rep(theta2power, arms-1), 
                                                 sigma = sigma,
                                                 interimN = stage,
                                                 alpha = alpha,
                                                 beta = beta,
                                                 r = ratio,
                                                 rho=rho,
                                                 prec = prec,
                                                 arms = arms,
                                                 ushape1 = ushape1,
                                                 lshape1 = ushape1,
                                                 ushape2 = ushape2,
                                                 lshape2 = ushape2,
                                                 power = power)

fun


popord <- fun$`sample size per arm per stage`

# CALL the function to simulate the data

data_ORD <- simul(scenario, nsim, stage, popord, ratio, rho, sigma, seed= seed, fun = fun)

data_ORD

##### Hierarchical 3-arm 2-stage design (S. Urach and M. Posch (2016)) ############

results <- boundpowerUrach3arm2stage(theta = c(theta2power,theta2power),
                                     sigma = sigma,
                                     interimN = stage,
                                     alpha = alpha,
                                     beta = beta,
                                     r = ratio,
                                     prec = prec,
                                     ushape = name,
                                     lshape = name)
results

u1 <- results$u1
u2 <- results$u2
l1 <- results$l1

v1 <- results$v1
v2 <- results$v2

popord <- results$`sample size per arm per stage`

##### Simulations #####

# Save rejections of H01 and H02

rej_h01andh02 <- list()

rej_h01andh02[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                             ncol=nsim)

rej_h01andh02[[2]] <- matrix(data=FALSE, nrow=nrow(scenario),
                             ncol=nsim)

# Save rejections of H01 

rej_h01 <- list()

rej_h01[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                       ncol=nsim)

rej_h01[[2]] <- matrix(data=FALSE, nrow=nrow(scenario),
                       ncol=nsim)

# Save rejections of H02

rej_h02 <- list()

rej_h02[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                       ncol=nsim)

rej_h02[[2]] <- matrix(data=FALSE, nrow=nrow(scenario),
                       ncol=nsim)


bothineff_all  <- matrix(data=FALSE, nrow=nrow(scenario),
                         ncol=nsim)

prop_rej <- list()

prop_rej[[1]] <- matrix(data=NA, ncol=3,
                        nrow = nrow(scenario))
prop_rej[[2]] <- matrix(data=NA, ncol=4,
                        nrow = nrow(scenario))

# Estimated sample size for each scenario

estimatedss_allp <- matrix(data=NA, nrow=nrow(scenario),
                           ncol=nsim)

estimatedsamplesize_allp <- matrix(data=NA, ncol=1,
                                   nrow = nrow(scenario))

for (k in 1:nrow(scenario)){
  
  set.seed(seed = seed)
  
  mean0 <- scenario[k,1]
  
  mean1 <- scenario[k,2]
  
  mean2 <- scenario[k,3]
  
  for (j in 1:nsim){
    
    # Generate total population for each scenario and each simulation
    
    plac1st <- rnorm(n = popord, mean = mean0, sd = sigma)
    
    dur1st <- rnorm(n = popord, mean = mean1, sd = sigma)
    
    dur2st <- rnorm(n = popord, mean = mean2, sd = sigma)
    
    estimatedss_allp[k,j] <- length(plac1st)+length(dur1st)+length(dur2st)
    
    # Estimated response rates on the observed data
    
    p_plac11 <- mean(plac1st)
    
    p_dur11 <- mean(dur1st)
    
    p_dur12 <- mean(dur2st)
    
    # Z-statistics for the two treatment durations
    
    z1 <- (p_dur11-p_plac11)/(sqrt(2/length(plac1st))*sigma)
    
    z2 <- (p_dur12-p_plac11)/(sqrt(2/length(plac1st))*sigma)
    
    # Count when null hypothesis are rejected at the end of the stage
    
    
    if(z1 >= u1 &  z2 < u1){ 
      
      rej_h01[[1]][k,j] = TRUE
      
      
    }
    
    if( (z1 >= v1 &  
         z2 >= u1)
        
    ){
      
      rej_h01[[1]][k,j] = TRUE
      
      rej_h01andh02[[1]][k,j] = TRUE
      
    }
    
    if(z2 >= v1 &  
       z1 >= u1
       
    ){
      
      rej_h01andh02[[1]][k,j] = TRUE
      
      rej_h02[[1]][k,j] = TRUE
      
    }
    
    
    if( z2 >= u1 &  z1 < u1){
      
      rej_h02[[1]][k,j] = TRUE
      
      
    }
    
    
    if( z1 < l1 &  z2 < l1){
      
      bothineff_all[k,j] = TRUE
      
      
    }
    
    
    if(z1 < u1  & z1 >= l1  & z2 >= l1 &
       z2 < u1 ){ # 
      
      # Estimated response rates
      
      plac21 <- rnorm(n = ceiling(popord*(ratio[2]-ratio[1])), mean = mean0, sd = sigma)
      
      dur21 <- rnorm(n = ceiling(popord*(ratio[2]-ratio[1])), mean = mean1, sd = sigma)
      
      dur22 <- rnorm(n = ceiling(popord*(ratio[2]-ratio[1])), mean = mean2, sd = sigma)
      
      estimatedss_allp[k,j] <-  estimatedss_allp[k,j]+length(plac21)*3
      
      plac21 <- c(plac21, plac1st)
      
      dur21 <- c(dur21, dur1st)
      
      dur22 <- c(dur22, dur2st)
      
      p_dur21 <- mean(dur21)
      
      p_plac21 <- mean(plac21)
      
      p_dur22 <- mean(dur22)
      
      # Z-statistics for the two treatment durations
      
      z21 <- (p_dur21-p_plac21)/(sqrt(2/length(plac21))*sigma)
      
      z22 <- (p_dur22-p_plac21)/(sqrt(2/length(plac21))*sigma)
      
      
      if (z21 >= v2 & z22 >= u2){
        
        rej_h01[[2]][k,j] = TRUE
        
        rej_h01andh02[[2]][k,j] = TRUE
        
        
      }
      
      if (z21 >= u2 & z22 < u2){
        
        rej_h01[[2]][k,j] = TRUE
        
        
      }
      if (z21 < u2 & z22 < u2){
        
        
        bothineff_all[k,j] = TRUE
        
        
        
      }
      
      if (z22 >= v2 & z21 >= u2){
        
        rej_h02[[2]][k,j] = TRUE
        
        rej_h01andh02[[2]][k,j] = TRUE
        
        
      }
      
      if (z22 >= u2 & z21 < u2){
        
        rej_h02[[2]][k,j] = TRUE
        
        
      }
      
      
      
    }
    
    if((z2 >= u1 & z1 < v1  & (z1 >= l1 )) ||
       (z1 >= l1 &  z1 < u1 & z2 < l1)){
      
      # Estimated response rates
      
      plac21 <- rnorm(n = ceiling(popord*(ratio[2]-ratio[1])), mean = mean0, sd = sigma)
      
      dur21 <- rnorm(n = ceiling(popord*(ratio[2]-ratio[1])), mean = mean1, sd = sigma)
      
      estimatedss_allp[k,j] <-  estimatedss_allp[k,j]+length(plac21)*2
      
      plac21 <- c(plac21, plac1st)
      
      dur21 <- c(dur21, dur1st)
      
      p_dur21 <- mean(dur21)
      
      p_plac21 <- mean(plac21)
      
      # Z-statistics for the two treatment durations
      
      z21 <- (p_dur21-p_plac21)/(sqrt(2/length(plac21))*sigma)
      
      
      if( z21 >= u2 & (z1 >= l1 &  z1 < u1 & z2 < l1)){
        
        rej_h01[[2]][k,j] = TRUE
        
        
      }
      
      
      if( z21 < u2 & z2 < l1){
        
        bothineff_all[k,j] = TRUE
        
        
      }
      
      
      if(z21 >= v2 & (z2 >= u1 & z1 < v1 & (z1 >= l1 )) ){
        
        rej_h01[[2]][k,j] = TRUE
        
        rej_h01andh02[[2]][k,j] = TRUE
        
        
      }
      
      
      
    }
    
    if((z1 >= u1 & z2 < v1 & (z2 >= l1 )) ||
       (z2 >= l1 &  z2 < u1 & z1 < l1)){ 
      
      # Estimated response rates
      
      plac21 <- rnorm(n = ceiling(popord*(ratio[2]-ratio[1])), mean = mean0, sd = sigma)
      
      dur22 <- rnorm(n = ceiling(popord*(ratio[2]-ratio[1])), mean = mean2, sd = sigma)
      
      estimatedss_allp[k,j] <-  estimatedss_allp[k,j]+length(plac21)*2
      
      plac21 <- c(plac21, plac1st)
      
      dur22 <- c(dur22, dur2st)
      
      p_dur22 <- mean(dur22)
      
      p_plac21 <- mean(plac21)
      
      # Z-statistics for the two treatment durations
      
      z22 <- (p_dur22-p_plac21)/(sqrt(2/length(plac21))*sigma)
      
      
      if( z22 >= u2 & (z2 >= l1 &  z2 < u1 & z1 < l1)){
        
        rej_h02[[2]][k,j] = TRUE
        
        
      }
      
      
      if( z22 < u2 & z1 < l1){
        
        bothineff_all[k,j] = TRUE
        
        
      }
      
      
      if(z22 >= v2 & (z1 >= u1 & z2 < v1 
                      & 
                      (z2 >= l1 )) ){
        
        rej_h02[[2]][k,j] = TRUE
        
        rej_h01andh02[[2]][k,j] = TRUE
        
        
      }
      
      
    }
  }
  
  
  # Average Probability of rejecting the null hypotheses
  
  prop_rej[[1]][k,1] <- sum(rej_h01andh02[[1]][k,])/nsim
  
  prop_rej[[1]][k,2] <- sum(rej_h01[[1]][k,])/nsim
  
  prop_rej[[1]][k,3] <- sum(rej_h02[[1]][k,])/nsim
  
  prop_rej[[2]][k,1] <- sum(rej_h01andh02[[2]][k,])/nsim
  
  prop_rej[[2]][k,2] <- sum(rej_h01[[2]][k,])/nsim
  
  prop_rej[[2]][k,3] <- sum(rej_h02[[2]][k,])/nsim
  
  prop_rej[[2]][k,4] <- sum(bothineff_all[k,])/nsim
  
  estimatedsamplesize_allp[k] <- mean(estimatedss_allp[k,]) 
}

colnames(prop_rej[[1]]) <- c(
  "Allprom_rejH01andH02_1stage", 
  "Allprom_rejH01_1stage",
  "Allprom_rejH02_1stage")

colnames(prop_rej[[2]]) <- c(
  "Allprom_rejH01andH02_2stage", 
  "Allprom_rejH01_2stage",
  "Allprom_rejH02_2stage",
  "Allprom_rejneitherH01andH02")

colnames(estimatedsamplesize_allp) <- "ESS All promising"

summary_PU <- cbind(popord,ceiling(popord*ratio[stage]*arms),scenario,
                    prop_rej, estimatedsamplesize_allp,
                    u1, u2, l1,
                    v1, v2)

colnames(summary_PU)[1:2] <- c("patients_perarm1stage ALLProm", "total sample ALLProm")

summary_PU <- summary_PU %>% mutate(
  
  Allprom_rejH01 = Allprom_rejH01_1stage+Allprom_rejH01_2stage,
  
  Allprom_rejH02 = Allprom_rejH02_1stage+Allprom_rejH02_2stage,
  
  Allprom_rejH01andH02 = Allprom_rejH01andH02_1stage+Allprom_rejH01andH02_2stage,
  
  Allprom_FWER = Allprom_rejH01+Allprom_rejH02-Allprom_rejH01andH02,
  
  sumprob = Allprom_FWER+Allprom_rejneitherH01andH02,
  
  u1, 
  u2, 
  l1,
  v1, 
  v2
  
) %>% select(
  
  `patients_perarm1stage ALLProm`,
  `total sample ALLProm`,
  `ESS All promising`,
  thetaH,
  thetaL,
  Allprom_rejH01,
  Allprom_rejH02,
  Allprom_rejH01andH02,
  Allprom_FWER,
  Allprom_rejneitherH01andH02,
  sumprob,
  u1,
  u2,
  l1,
  v1,
  v2
  
)

##### Bayesian design ####

##### Find thresholds and sample size ####

# Parameters of interest from the paper: The Tiotropium in asthmatic adolescents symptomatic
# despite inhaled corticosteroids: A randomised dose-ranging study

diff <- 120
sd <- 340
alpha <- 0.025
beta <- 0.2

# Libraries 

library(dplyr)
library(mvtnorm)

#prior hyperparameters for model 

tau <- 1/(sd^2)
tau01 <- 10^-6
tau00 <- 0.00039
mu01 <- 602 
mu00 <- 489

# grid of values for meand and taud

taudv <- 10^-6
meandv <- 0

trteffect <- diff

# value/s of mu^(0) for which we want to control the FWER

null <- seq(mu00-30,mu00+35,by=5) #mu00

# thresholds

eta1 <- 0.9906 

#sample size

n11 <- 82

##### theoretical  simulations #######

# grid of values for sample sizes and critical bounds

grid_n11eta <- expand.grid(n11, eta1 )
colnames(grid_n11eta) <- c("n11","eta1" )

# for OBF bounds: eta2 = pnorm(qnorm(eta)/sqrt(2))
# for TRIAN bounds: ratio <- c(1,2) eta2 = pnorm((sqrt(ratio)*qnorm(eta1)/(1+ratio/max(ratio)))[1]*(1 + ratio/max(ratio))/sqrt(ratio))[2]
#                                   eps1 = pnorm((sqrt(ratio)*qnorm(eta1)/(1+ratio/max(ratio)))[1]*-((1-3*ratio/max(ratio))/sqrt(ratio))[1])
# For Pocock bounds: eta2=eta1, eps1 = -eta1

grid_n11eta <- grid_n11eta %>% mutate(
  n21=n11,n01=n11, 
  eps1=pnorm((qnorm(eta1)/(1+1/2))*-((1-3/2))), 
  eta2= pnorm((qnorm(eta1)/(1+1/2))*(1 + 1)/sqrt(2))
) %>% select(
  n11,n21,n01,eta1, eps1, eta2
  )

# Grid of values for different values of meand and taud and null configurations

grid <- cbind(expand.grid(taudv,meandv),null)

# Apply the function integ_function to each row of the grid

solfin <- apply(grid,1,integ_function, grid_n11eta=grid_n11eta)

# Results matrix

solfinmatrix <- do.call(rbind.data.frame, solfin)

tol <- 0

# Filter the results that meet the FWER and power requirements

powerrejboth <- 1-beta

filtersolcontrol <- solfinmatrix %>% filter(
  abs(rejatleastone-tol)<=alpha &  
    abs(rejectpowernull1-tol) <=alpha,
  abs(rejectpowernull2-tol) <=alpha &  
    rejpower >= powerrejboth
)

##### simulations ###

library(dplyr)
library(rjags)

# choose which model to use 

modelname <- "bayesianModel.bug"

#prior hyperparameters for model 

tau <- 1/(340^2)
tau01 <- 10^-6
tau00 <- 10^-6
mu01 <- 602 
mu00 <- 489
taud <- 10^-6
delta00 <- 0
trteffect <- diff <- 120

# sample size per arm per stage

n1 <- 102
n2 <- 102
n0 <- 102

# Thresholds 

eta1 <- 0.99347
eta2 <- 0.990363
eps1 <- 0.795988

# PARAMETERS for simulations

seed <- 64736
nsim <- 10^4

# Simulations scenarios

grid <- c(0,20,40,60,80,120)
grid2 <- expand.grid(grid,grid,grid)+mu00

scenario <- unique(data.frame(grid2)) 

colnames(scenario) <- c("mu1s", "mu2s", "mu0s")

scenario <- scenario %>% dplyr::filter(
  (mu0s ==489 & ((mu1s==609 | mu2s == 609) | (mu1s==489 & mu2s == 489)))

)

##### simulations with MCMC #######

mcmcrun <- function(modelname, data){
  
  if(modelname=="bayesianModel.bug"){
    
    inits <- list(list("mu1"=rnorm(1, 0, 1),
                       "delta"=rnorm(1, 0, 1),
                       "mu0"=rnorm(1, 0, 1)))
    
    
    parameters <- c("theta1", "theta2",
                    "mu1",
                    "mu2",
                    "mu0",
                    "delta")
  }
  
  setwd(".")
  
  model <- jags.model(modelname,
                      data=data, inits=inits, n.chains=1)
  
  update(model, n.iter=5000)
  
  samples <- coda.samples(model, 
                          variable.names=parameters,
                          n.iter=20000)
  
  ## P(theta|y)>0
  
  ptheta1eff <- mean(samples[[1]][,"theta1"]>0)
  
  ptheta2eff <- mean(samples[[1]][,"theta2"]>0)
  
  output <- c(ptheta1eff,ptheta2eff,
              samples)
  
  names(output) <- c("stop_efficacy1", "stop_efficacy2",
                     "sample")
  return(output)
}

stopefficacy12 <- stopefficacy1not2 <- stopfutility12 <- stopefficacy1not2 <- stopefficacy2not1 <- list()

stopefficacy12[[1]] <- stopefficacy1not2[[1]] <- stopefficacy2not1[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                                                                                  ncol=nsim)

stopfutility12[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                              ncol=nsim)

stopefficacy12[[2]] <- stopefficacy1not2[[2]] <- stopefficacy2not1[[2]] <-matrix(data=FALSE, nrow=nrow(scenario),
                                                                                 ncol=nsim)

stopfutility12[[2]] <- matrix(data=FALSE, nrow=nrow(scenario),
                              ncol=nsim)
prop_ord <- list()

prop_ord[[1]] <- prop_ord[[2]]  <- matrix(data=NA, ncol=4,
                                          nrow = nrow(scenario))

prejatleastone  <- matrix(data=FALSE, nrow=nrow(scenario),
                          ncol=nsim)

estimatedss <- matrix(data=NA, nrow=nrow(scenario),
                      ncol=nsim)

estimatedsamplesize <- matrix(data=NA, ncol=1,
                              nrow = nrow(scenario))

for (k in 1:nrow(scenario)){
  
  set.seed(seed = seed)
  
  mu1s <- scenario[k,1]
  mu2s <- scenario[k,2]
  mu0s <- scenario[k,3]
  
  for (j in 1:nsim){
    
    data1 <- list(
      "y1"=rnorm(n1, mu1s, sqrt((1/tau))),
      "y2"=rnorm(n2, mu2s, sqrt((1/tau))),
      "y0"=rnorm(n0, mu0s, sqrt((1/tau))),
      "mu00"=mu00,
      "mu01"=mu01,
      "n1"=n1,
      "n2"=n2,
      "n0"=n0,
      "tau01"=tau01,
      "tau00"=tau00,
      "tau"=tau,
      "taud"=taud,
      "delta00"=delta00
      
    )
    
    out1 <- mcmcrun(modelname =modelname, 
                    data = data1)
    
    estimatedss[k,j] <- n1+n2+n0
      
      if(out1["stop_efficacy1"] >= eta1 & out1["stop_efficacy2"] < eps1){
        
        stopefficacy1not2[[1]][k,j] <- TRUE
      }
      if(out1["stop_efficacy2"] >= eta1 & out1["stop_efficacy1"] < eps1){
        
        stopefficacy2not1[[1]][k,j] <- TRUE
      }
      if(out1["stop_efficacy1"] >= eta1 & out1["stop_efficacy2"] >= eta1){
        
        stopefficacy12[[1]][k,j] <- TRUE
      }
      
      if(out1["stop_efficacy1"] < eps1 & out1["stop_efficacy2"] < eps1){
        
        stopfutility12[[1]][k,j] <- TRUE
      }
      
      # stop arm 2 for futility and continue with arm 1 and control
      
      if(out1["stop_efficacy1"] >= eps1 &
         out1["stop_efficacy1"] < eta1 &
         (out1["stop_efficacy2"]< eps1 || out1["stop_efficacy2"]>= eta1)){
        
        data2 <- list(
          "y1"=c(data1$y1, rnorm(n1, mu1s, sqrt((1/tau)))),
          "y2"=c(data1$y2),
          "y0"=c(data1$y0, rnorm(n0, mu0s, sqrt((1/tau)))),
          "mu00"=mu00,
          "mu01"=mu01,
          "n1"=2*n1,
          "n2"=n2,
          "n0"=2*n0,
          "tau01"=tau01,
          "tau00"=tau00,
          "tau"=tau,
          "taud"=taud,
          "delta00"=delta00
          
        )
        
        out2 <- mcmcrun(modelname =modelname, data = data2)
        
        estimatedss[k,j] <- estimatedss[k,j]+ n1+n0
        
        if(out2["stop_efficacy1"]>= eta2 & 
           out1["stop_efficacy2"]< eps1){
          
          stopefficacy1not2[[2]][k,j] <- TRUE
        }
        
        if(out2["stop_efficacy1"]>= eta2 & 
           out1["stop_efficacy2"]>= eta1){
          
          stopefficacy12[[2]][k,j] <- TRUE
        }
        
        if(out2["stop_efficacy1"]< eta2 & 
           out1["stop_efficacy2"]< eps1){
          
          stopfutility12[[2]][k,j] <- TRUE
        }
        
        if(out2["stop_efficacy1"]< eta2 & 
           out1["stop_efficacy2"]>= eta1){
          
          stopefficacy2not1[[2]][k,j] <- TRUE
        }
        
        
      }
      
      # stop arm 1 for futility and continue with arm 2 and control
      
      if(out1["stop_efficacy2"] >= eps1 &
         out1["stop_efficacy2"] < eta1 &
         (out1["stop_efficacy1"]< eps1 || out1["stop_efficacy1"]>= eta1)){
        
        data2 <- list(
          "y1"=c(data1$y1),
          "y2"=c(data1$y2, rnorm(n2, mu2s, sqrt((1/tau)))),
          "y0"=c(data1$y0, rnorm(n0, mu0s, sqrt((1/tau)))),
          "mu00"=mu00,
          "mu01"=mu01,
          "n1"=n1,
          "n2"=2*n2,
          "n0"=2*n0,
          "tau01"=tau01,
          "tau00"=tau00,
          "tau"=tau,
          "taud"=taud,
          "delta00"=delta00
          
        )
        
        out2 <- mcmcrun(modelname =modelname, data = data2)
        
        estimatedss[k,j] <- estimatedss[k,j]+ n2+n0
        
        if(out2["stop_efficacy2"]>= eta2 & 
           out1["stop_efficacy1"]< eps1){
          
          stopefficacy2not1[[2]][k,j] <- TRUE
        }
        
        if(out2["stop_efficacy2"]>= eta2 & 
           out1["stop_efficacy1"] >= eta1){
          
          stopefficacy12[[2]][k,j] <- TRUE
        }
        
        if(out2["stop_efficacy2"]< eta2 & 
           out1["stop_efficacy1"]< eps1){
          
          stopfutility12[[2]][k,j] <- TRUE
        }
        
        if(out2["stop_efficacy2"]< eta2 &
           out1["stop_efficacy1"] >= eta1){
          
          stopefficacy1not2[[2]][k,j] <- TRUE
        }
        
        
      }
      
      # continue with both arms
      
      if(out1["stop_efficacy2"] >= eps1 &
         out1["stop_efficacy2"] < eta1 &
         out1["stop_efficacy1"] >= eps1 &
         out1["stop_efficacy1"] < eta1 ){
        
        data2 <- list(
          "y1"=c(data1$y1, rnorm(n1, mu1s, sqrt((1/tau)))),
          "y2"=c(data1$y2, rnorm(n2, mu2s, sqrt((1/tau)))),
          "y0"=c(data1$y0, rnorm(n0, mu0s, sqrt((1/tau)))),
          "mu00"=mu00,
          "mu01"=mu01,
          "n1"=2*n1,
          "n2"=2*n2,
          "n0"=2*n0,
          "tau01"=tau01,
          "tau00"=tau00,
          "tau"=tau,
          "taud"=taud,
          "delta00"=delta00
          
        )
        
        out2 <- mcmcrun(modelname =modelname, data = data2)
        
        estimatedss[k,j] <- estimatedss[k,j]+ n1+n2+n0
        
        if(out2["stop_efficacy2"]>= eta2 & 
           out2["stop_efficacy1"]>= eta2){ 
          
          stopefficacy12[[2]][k,j] <- TRUE
        }
        
        if(out2["stop_efficacy1"]>= eta2 & 
           out2["stop_efficacy2"]< eta2
        ){
          
          stopefficacy1not2[[2]][k,j] <- TRUE
        }
        
        if(out2["stop_efficacy2"]>= eta2 & 
           out2["stop_efficacy1"]< eta2 
        ){
          
          stopefficacy2not1[[2]][k,j] <- TRUE
        }
        
        if(out2["stop_efficacy1"]< eta2 & 
           out2["stop_efficacy2"]< eta2
        ){
          
          stopfutility12[[2]][k,j] <- TRUE
        }
        
        
      }
      
  }
    
    prop_ord[[1]][k,1] <- sum(stopfutility12[[1]][k,])/nsim
    
    prop_ord[[1]][k,2] <- sum(stopefficacy12[[1]][k,])/nsim
    prop_ord[[1]][k,3] <- sum(stopefficacy1not2[[1]][k,])/nsim
    prop_ord[[1]][k,4] <- sum(stopefficacy2not1[[1]][k,])/nsim
    
    
    prop_ord[[2]][k,1] <- sum(stopfutility12[[2]][k,])/nsim
    
    prop_ord[[2]][k,2] <- sum(stopefficacy12[[2]][k,])/nsim
    prop_ord[[2]][k,3] <- sum(stopefficacy1not2[[2]][k,])/nsim
    prop_ord[[2]][k,4] <- sum(stopefficacy2not1[[2]][k,])/nsim
  
  
   estimatedsamplesize[k] <- mean(estimatedss[k,])
}

  colnames(prop_ord[[1]]) <- c(
    "Stop_futility12_1stage",
    "Stop_efficacy12_1stage",
    "Stop_efficacy1not2_1stage",
    "Stop_efficacy2not1_1stage")
  
  colnames(prop_ord[[2]]) <- c(
    "Stop_futility12_2stage",
    
    "Stop_efficacy12_2stage",
    "Stop_efficacy1not2_2stage",
    "Stop_efficacy2not1_2stage")
  
  summary_bayesian <- cbind(scenario, prop_ord, n1,n2, n0, estimatedsamplesize,
                           eta1,eps1,eta2, tau01,tau00, 
                           taud,delta00,mu00,mu01) %>% dplyr::mutate(
                             
                             rejefficacy_atleastone = (Stop_efficacy1not2_1stage+Stop_efficacy1not2_2stage)+
                               (Stop_efficacy2not1_1stage+Stop_efficacy2not1_2stage)+
                               (Stop_efficacy12_1stage+Stop_efficacy12_2stage),
                             
                             Pr_stop_futility12 = (Stop_futility12_1stage+Stop_futility12_2stage),
                             
                             Pr_stop_futility12_1stage = (Stop_futility12_1stage),
                             
                             Pr_stop_efficacy1 = (Stop_efficacy12_1stage+Stop_efficacy12_2stage)+
                               (Stop_efficacy1not2_1stage+Stop_efficacy1not2_2stage),
                             
                             Pr_stop_efficacy2 = (Stop_efficacy12_1stage+Stop_efficacy12_2stage)+
                               (Stop_efficacy2not1_1stage+Stop_efficacy2not1_2stage),
                             
                             Pr_stop_efficacy12 = (Stop_efficacy12_1stage+Stop_efficacy12_2stage),
                             
                             Pr_stop_efficacy1not2 = (Stop_efficacy1not2_1stage+Stop_efficacy1not2_2stage),
                             
                             Pr_stop_efficacy2not1 = (Stop_efficacy2not1_1stage+Stop_efficacy2not1_2stage),
                             
                             Pr_stop_efficacy1not2_1stage = Stop_efficacy1not2_1stage,
                             
                             Pr_stop_efficacy2not1_1stage = Stop_efficacy2not1_1stage,
                             
                             Pr_stop_efficacy12_1stage = Stop_efficacy12_1stage,
                             
                             Pr_stop_efficacy12_2stage = Stop_efficacy12_2stage,
                             
                             ssarm1stage = n1,
                             
                             ssarm2stage = n2,
                             
                             sscontrolstage = n0,
                             
                             totss = (n1+n2)*2+n0*2,
                             
                             sumprob = rejefficacy_atleastone + Pr_stop_futility12,
                             
                             theta1 = mu1s-mu0s,
                             
                             theta2 = mu2s-mu0s
                             
                           ) %>% dplyr::select(
                             mu1s,
                             mu2s,
                             mu0s,
                             theta1,
                             theta2,
                             ssarm1stage,
                             ssarm2stage,
                             sscontrolstage,
                             totss,
                             estimatedsamplesize,
                             rejefficacy_atleastone,
                             Pr_stop_futility12,
                             Pr_stop_futility12_1stage,
                             Pr_stop_efficacy1,
                             Pr_stop_efficacy2,
                             Pr_stop_efficacy12,
                             Pr_stop_efficacy12_1stage,
                             Pr_stop_efficacy12_2stage,
                             Pr_stop_efficacy1not2_1stage,
                             Pr_stop_efficacy1not2,
                             Pr_stop_efficacy2not1,
                             Pr_stop_efficacy2not1_1stage, 
                             sumprob,
                             eta1,
                             eps1,
                             eta2,
                             tau01,
                             tau00,
                             taud,
                             delta00,
                             mu00,
                             mu01)



  
