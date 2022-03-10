####### Functions MAMS(m) ####
##### Function to find critical bounds for 3-arm or 4-arm J-stage MAMS(m) ##############

#' @u1: parameter for the grid search
#' @interimN: number of interim analyses
#' @alpha: alpha-level of the hypothesis test
#' @cov: covariance matrix
#' @first: covariance matrix at the first stage
#' @r: allocation ratio sample size in the first and second stage
#' @prec: precision for the boundaries: numbers after the comma
#' @arms: number of arms (3)
#' @ushape: upper boundary shape: "pocock", "obf", "triangular"
#' @lshape: lower boundary shape: "pocock", "obf", "triangular"
#' @lfix: fixed value for lower bound

bounds_3armJstagem <- function(u1, interimN, arms,r,
                               alpha, cov, first, ushape,lshape,lfix, prec){
  
  library(mvtnorm)
  library(gtools)
  
  alphat <- rep(0,length = interimN)
  
  alphatcaseup<- NULL
  
  alphatcasedown<- NULL
  
  alphatcasefin<- NULL
  
  for(j in 1:length(u1)){
    
    for (i in 1:interimN){
      
      
      L <- array(0, dim = i*(arms-1))
      U <- array(0, dim = i*(arms-1))
      
      if (ushape == "obf") {
        u <- u1[j] * 1/sqrt(r)
      }
      else if (ushape == "pocock") {
        u <- rep(u1[j], i)
      }
      else if (ushape == "triangular") {
        u <- u1[j] * (1 + r/max(r))/sqrt(r)
      }
      
      if (lshape == "obf") {
        l <- c(-u1[j] * 1/sqrt(r))[1:(i-1)]
      }
      else if (lshape == "pocock") {
        l <- c(rep(-u1[j], i-1 ))
      }
      else if (lshape == "triangular") {
        if (ushape == "triangular") {
          l <- c(-u1[j] * ((1 -3* r/max(r))/sqrt(r)))[1:(i-1)]
        }
      }
      else if (lshape == "fixed") {
        l <- c(rep(lfix, i-1 ))
      }
      
      
      if(i ==1){
        set.seed(123456)
        
        alphat[i] <- 1-pmvnorm(lower = rep(-Inf,times = arms-1),
                               upper = rep(u[1],times = arms-1),
                               sigma = first)[1]
        
      }
      if(i>1){
        
        st <- seq(1, i*(arms-1), by=(arms-1))
        
        st2 <- seq(2, i*(arms-1), by=(arms-1))
        
        L[st] <- c(l, u[i])
        
        U[st[1:(i-1)]] <- u[1:(i-1)]
        
        L[-st] <- rep(-Inf, times = arms-2)
        
        o <- 2 
        it <- 1
        
        while(o<max(st2)){
          
          U[o:(it*(arms-1))] <- u[1]
          
          it <- it+1
          o <- st2[it]
        }
        
        U[which(U==0)] <- Inf
        
        L <- c(L,rep(-Inf, times = ncol(cov)-length(L)))
        
        U <- c(U,rep(Inf, times = ncol(cov)-length(U)))
        
        
        index <- seq(2, ncol(cov), by=arms-1)
        
        Lneg <- L
        
        Lneg[index] <- Lneg[index-1]
        
        Uneg <- U
        
        Uneg[index] <- Uneg[index-1]
        
        int <- seq(0, ncol(cov), by=(arms-1))
        
        permL <- matrix(0, ncol = ncol(cov), nrow = (arms-1))
        
        permU <- matrix(0, ncol = ncol(cov), nrow = (arms-1))
        
        permLneg <- matrix(0, ncol = ncol(cov), nrow = (arms-1)*(arms-2)/2)
        
        permUneg <- matrix(0, ncol = ncol(cov), nrow = (arms-1)*(arms-2)/2)
        
        s <- 1
        
        for(p in 2:length(int)){
          if(length(unique(L[s:int[p]]))==1){
            permL[,s:int[p]] <- permutations(arms-1,arms-1,L[s:int[p]], set=FALSE)[1:nrow(permL),]
          }
          
          if(length(unique(U[s:int[p]]))==1){
            permU[,s:int[p]] <- permutations(arms-1,arms-1,U[s:int[p]], set=FALSE)[1:nrow(permU),]
          }
          
          if(length(unique(Lneg[s:int[p]]))==1){
            permLneg[,s:int[p]]  <-permutations(arms-1,arms-1,Lneg[s:int[p]], set=FALSE)[1:nrow(permLneg),]
          }
          
          
          if(length(unique(Uneg[s:int[p]]))==1){
            permUneg[,s:int[p]]  <-permutations(arms-1,arms-1,Uneg[s:int[p]], set=FALSE)[1:nrow(permUneg),]
          }
          
          if(length(unique(L[s:int[p]]))!=1){
            permL[,s:int[p]]  <- unique(permutations(arms-1,arms-1,L[s:int[p]], set=FALSE))
            
          }
          
          if(length(unique(U[s:int[p]]))!=1){
            permU[,s:int[p]]  <- unique(permutations(arms-1,arms-1,U[s:int[p]], set=FALSE))
            
          }
          
          if(length(unique(Lneg[s:int[p]]))!=1){
            permLneg[,s:int[p]]  <- unique(permutations(arms-1,arms-1,Lneg[s:int[p]], set=FALSE))
          }
          
          if(length(unique(Uneg[s:int[p]]))!=1){
            permUneg[,s:int[p]]  <- unique(permutations(arms-1,arms-1,Uneg[s:int[p]], set=FALSE))
          }
          s <- int[p]+1
          
        }
        
        permLneg <- unique(permLneg)
        
        permL <- unique(permL)
        
        index2 <- seq(3, ncol(cov), by=arms-1)
        
        Lfin <- Lneg
        
        Lfin[index2] <- Lfin[index2-1]
        
        Ufin <- Uneg
        
        Ufin[index2] <- Ufin[index2-1]
        
        for(c in 1:nrow(permL)){
          
          set.seed(123456)
          alphatcaseup[c] <- pmvnorm(lower = as.vector(permL[c,]),
                                     upper = as.vector(permU[c,]),
                                     sigma = cov)[1]
        }
        for(c in 1:nrow(permLneg)){
          
          set.seed(123456)
          alphatcasedown[c] <- pmvnorm(lower = as.vector(permLneg[c,]),
                                       upper = as.vector(permUneg[c,]),
                                       sigma = cov)[1]
        }
        if(arms>3){
          set.seed(123456)
          alphatcasefin <- pmvnorm(lower = as.vector(Lfin),
                                   upper = as.vector(Ufin),
                                   sigma = cov)[1]
        }
        
        alphat[i] <- sum(alphatcaseup)-sum(alphatcasedown)+sum(alphatcasefin)
        
        
      }
      
      
    }
    
    if(sum(alphat) <= alpha){
      
      if(i==1){
        
        upperbound <- u[1]
        lowerbound <- u[1]
        totalpha <- sum(alphat)
        alow <- u1[j]
        finalalpha <- alphat
        
      }
      else{
        upperbound <- u
        lowerbound <- l
        totalpha <- sum(alphat)
        alow <- u1[j]
        finalalpha <- alphat
      }
      
    }
    
    if(sum(alphat) > alpha){
      
      aup <- u1[j]
      break
      
      
    }
    
  }
  
  results <- list(upperbound,alow, aup, lowerbound,totalpha)
  
  names(results) <- c("upperbound", "alow", "aup", "lowerbound", "totalpha")
  return(results)
  
}

##### Function to find sample size for 3-arm 2-stage or 3-arm 3-stage design ###############

#' @theta: vector of clinically relevant difference 
#' @sigma: standard deviation of the population
#' @interimN: number of interim analyses
#' @alpha: alpha-level of the hypothesis test
#' @beta: beta-level of the test
#' @r: allocation ratio sample size in the first and second stage
#' @rho: allocation ratio sample size in each treatment divided by sample size in the control
#' @prec: precision for the boundaries: numbers after the comma
#' @arms: number of arms (3)
#' @ushape: upper boundary shape : "pocock", "obf", "triangular"
#' @lshape: lower boundary shape : "pocock", "obf", "triangular"
#' @lfix: fixed value for lower bound

boundaries_3armJstagem <- function(theta,
                                   sigma, 
                                   interimN, 
                                   alpha, 
                                   beta, 
                                   r, 
                                   rho,
                                   prec,
                                   arms,
                                   ushape,
                                   lshape,
                                   lfix){ 
  
  library(mvtnorm)
  library(cubature)
  
  #construction of the covariance matrix
  
  R <- rho[1:(length(rho)-1)]/(rho[1:(length(rho)-1)]+rho[length(rho)])
  
  bottom<-matrix(r,length(r),length(r))
  top<-matrix(rep(r,rep(length(r),length(r))),length(r),length(r))
  top[upper.tri(top)]<-t(top)[upper.tri(top)]
  bottom[upper.tri(bottom)]<-t(bottom)[upper.tri(bottom)] 
  
  cov <-sqrt(top/bottom)
  
  ma <- matrix(0, nrow = length(R), ncol = length(R))
  
  for (i in 1:length(R)){
    
    for(j in 2:length(R)){
      
      ma[i,j] <- sqrt(R[i]*R[j])
    }
  }
  
  diag(ma) <- 1
  
  ma[lower.tri(ma)]  <- t(ma)[lower.tri(ma)]
  
  block <- list()
  
  dim <- 1
  
  for (i in 1:nrow(cov)){
    
    for(j in 1:ncol(cov)){
      
      block[[dim]] <- cov[i,j]*ma
      
      dim <- dim+1
    }
    
  }
  
  covb <- NULL
  
  for (i in 1:length(block)){
    
    covb <- cbind(covb, block[[i]])
  }
  
  int <- seq(0, ncol(covb), by=((arms-1)*length(r)))
  
  covbfin <- NULL
  s <- 1
  for(j in 2:length(int)){
    covbfin <- rbind(covbfin, covb[,s:int[j]])
    s <- int[j]+1
  }
  
  diag(covbfin) <- 1
  
  #covariance matrix at the first stage
  
  firststagematrix <- covbfin[1:(arms-1),1:(arms-1)]
  
  diag(firststagematrix) <- 1
  
  ## Research of the critical bounds
  
  u1 <- seq(from = 7, to = 0, by = -1)
  
  first <-bounds_3armJstagem(u1, interimN, arms,r, 
                             alpha, covbfin, firststagematrix,
                             ushape, lshape,lfix,prec)
  
  low <- first$alow
  
  up <- first$aup
  
  p <- rep(1, times = prec)
  
  for (p in 1:length(p)){
    
    callf <- bounds_3armJstagem(seq(from = low, 
                                    to = up, by = -1/(10^p)), 
                                interimN= interimN,
                                arms,
                                r,
                                alpha,
                                covbfin,
                                firststagematrix,
                                ushape,
                                lshape,
                                lfix,
                                prec)
    
    low <- callf$alow
    
    up <- callf$aup
  }
  
  
  upperbound <- callf$upperbound
  lowerbound <- callf$lowerbound
  a <- callf$alow
  totalalpha <- callf$totalpha
  
  # Search of the sample size to reach the desired power
  
  maxsample <- NULL
  
  betat <- NULL
  
  betatcaseup <- NULL
  
  betatcasedown <- NULL
  
  betatcasefin <- NULL
  
  pop <- seq(from = 1, to = 1500, by = 1)
  
  if (length(theta)!=(arms-1)){
    print("Error: the length of theta should be equal to arms-1")
  }
  
  rarm <- c(rep(r, each=arms-1))
  
  for (n in 1:length(pop)){
    
    mean2 <- (theta)*(sqrt(rarm*R*pop[n]))/(sigma)
    
    mean1 <- mean2[1:2] #at the first stage
    
    for (i in 1:interimN){
      
      L <- array(0, dim = i*(arms-1))
      U <- array(0, dim = i*(arms-1))
      
      if (ushape == "obf") {
        u <- a * 1/sqrt(r)
      }
      else if (ushape == "pocock") {
        u <- rep(a, i)
      }
      else if (ushape == "triangular") {
        u <- a * (1 + r/max(r))/sqrt(r)
      }
      
      if (lshape == "obf") {
        l <- c(-a * 1/sqrt(r))[1:(i-1)]
      }
      else if (lshape == "pocock") {
        l <- c(rep(-a, i-1 ))
      }
      else if (lshape == "triangular") {
        if (ushape == "triangular") {
          l <- c(-a * ((1 -3* r/max(r))/sqrt(r)))[1:(i-1)]
        }
      }
      else if (lshape == "fixed") {
        l <- c(rep(lfix, i-1 ))
      }
      
      
      # For 3-arm 2-stage or 3-arm 3-stage MAMS non-binding constraints
      
      
      if(i ==1){
        set.seed(123456)
        
        betat[i] <- pmvnorm(lower = rep(upperbound[1],times = arms-1),
                            upper = rep(Inf,times = arms-1),
                            sigma = firststagematrix,
                            mean = mean1)[1]
        
        
      }
      else{
        
        st <- seq(1, i*(arms-1), by=(arms-1))
        
        st2 <- seq(2, i*(arms-1), by=(arms-1))
        
        L[st] <- c(l, u[i])
        
        U[st[1:(i-1)]] <- u[1:(i-1)]
        
        L[-st] <- rep(-Inf, times = arms-2)
        
        o <- 2 
        it <- 1
        
        while(o<max(st2)){
          
          U[o:(it*(arms-1))] <- u[1]
          
          it <- it+1
          o <- st2[it]
        }
        
        U[which(U==0)] <- Inf
        
        L <- c(L,rep(-Inf, times = ncol(covbfin)-length(L)))
        
        U <- c(U,rep(Inf, times = ncol(covbfin)-length(U)))
        
        index <- seq(2, ncol(covbfin), by=arms-1)
        
        Lneg <- L
        
        Lneg[index] <- Lneg[index-1]
        
        Uneg <- U
        
        Uneg[index] <- Uneg[index-1]
        
        y <- seq(2, length(L), by =2)
        
        L1 <- Lneg
        
        U1 <- Uneg
        
        L2 <- Lneg
        
        U2 <- Uneg
        
        L1[y] <- c(u[1], rep(-Inf, (length(L)/2)-1))
        
        U1[y] <- rep(Inf, (length(U)/2))
        
        L2[y-1] <- c(u[1], rep(-Inf, (length(L)/2)-1))
        
        U2[y-1] <- rep(Inf, (length(U)/2))
        
        betat[i] <- pmvnorm(lower = as.vector(Lneg)[1:(2*i)],
                            upper = as.vector(Uneg)[1:(2*i)],
                            sigma = covbfin[1:(2*i),1:(2*i)],
                            mean = mean2[1:(2*i)])[1]+
          pmvnorm(lower = as.vector(L1)[1:(2*i)],
                  upper = as.vector(U1)[1:(2*i)],
                  sigma = covbfin[1:(2*i),1:(2*i)],
                  mean = mean2[1:(2*i)])[1]+
          pmvnorm(lower = as.vector(L2)[1:(2*i)],
                  upper = as.vector(U2)[1:(2*i)],
                  sigma = covbfin[1:(2*i),1:(2*i)],
                  mean = mean2[1:(2*i)])[1]
        if(i>2){
          
          betat[i] <- betat[i] +
            pmvnorm(lower = c(l[1], l[1], l[2], u[2], u[3], -Inf),
                    upper = c(u[1],u[1], u[2], Inf, Inf, Inf),
                    sigma = covbfin[1:(2*i),1:(2*i)],
                    mean = mean2[1:(2*i)])[1]+
            pmvnorm(lower = c(l[1], l[1],  u[2],l[2], -Inf, u[3]),
                    upper = c(u[1],u[1],  Inf,u[2], Inf, Inf),
                    sigma = covbfin[1:(2*i),1:(2*i)],
                    mean = mean2[1:(2*i)])[1]
          
        }
        
      }
      
      
      
      
      
      
    }
    
    
    
    
    
    if(sum(betat) > 1-beta){
      
      maxsample <- pop[n]
      break
      
    }
  }
  
  results <- list(totalalpha,
                  lowerbound,upperbound,maxsample )
  
  names(results) <- c("sum_alphat", 
                      "lowerbounds", "upperbounds", "sample size per arm per stage")
  
  return(results)
}




##### Function to simulate 3-arm J-stage MAMS(m) #########

simulAll <- function(scenario, nsim, stage, popord, ratio,rho, sigma, seed, func){
  
  # save rejection first and/or second hypothesis in each stage and
  # save rejection for the power computation when theta = (theta1, theta2) and
  # theta1 >= theta2 > delta0 > 0
  
  powerII <- matrix(data=FALSE, nrow=nrow(scenario),
                    ncol=nsim)
  
  bothineff <- matrix(data=FALSE, nrow=nrow(scenario),
                      ncol=nsim)
  
  rej_ord_h01andh02 <- matrix(data=FALSE, nrow=nrow(scenario),
                              ncol=nsim)
  
  rej_ord_h01noh02 <- matrix(data=FALSE, nrow=nrow(scenario),
                             ncol=nsim)
  
  rej_ord_h02noh01 <- matrix(data=FALSE, nrow=nrow(scenario),
                             ncol=nsim)
  
  prop_ord <- matrix(data=NA, ncol=5,
                     nrow = nrow(scenario))
  
  # Estimated sample size for each scenario. --> The average of the sample size used for each simulation
  
  estimatedss <- matrix(data=NA, nrow=nrow(scenario),
                        ncol=nsim)
  
  estimatedsamplesize <- matrix(data=NA, ncol=1,
                                nrow = nrow(scenario))
  
  for (k in 1:nrow(scenario)){
    
    set.seed(seed = seed)
    
    mean0 <- scenario[k,1]
    
    mean1 <- scenario[k,2]
    
    mean2 <- scenario[k,3]
    
    for (j in 1:nsim){
      
      continue <- TRUE
      
      plac1st <- NULL
      
      dur1st <- NULL
      
      dur2st <- NULL
      
      cell2 <- FALSE
      
      cell4 <- FALSE
      
      cell3 <- FALSE
      
      cell6 <- FALSE
      
      # Simulation of MAMS(m)
      
      for(s in 1:stage){
        
        if (continue==TRUE){
          
          if(s==1|| (cell5==TRUE)){
            
            continue <- FALSE
            
            
            cell5 <- FALSE
            
            if(s==1){
              
              plac1 <- rnorm(n = popord*(ratio[s]*rho[3]), mean = mean0, sd = sigma)
              
              dur1 <- rnorm(n = popord*(ratio[s]*rho[1]), mean = mean1, sd = sigma)
              
              dur2 <- rnorm(n = popord*(ratio[s]*rho[2]), mean = mean2, sd = sigma)
            }
            
            if(s>1){
              
              plac1 <- rnorm(n = ceiling(popord*(ratio[s]-ratio[s-1])*rho[3]), mean = mean0, sd = sigma)
              
              dur1 <- rnorm(n = ceiling(popord*(ratio[s]-ratio[s-1])*rho[1]), mean = mean1, sd = sigma)
              
              dur2 <- rnorm(n = ceiling(popord*(ratio[s]-ratio[s-1])*rho[2]), mean = mean2, sd = sigma)
            }
            
            # Estimated response rates on the observed data
            
            plac1st <- c(plac1, plac1st)
            
            dur1st <- c(dur1, dur1st)
            
            dur2st <- c(dur2, dur2st)
            
            estimatedss[k,j] <- length(plac1st)+length(dur1st)+length(dur2st)
            
            p_plac11 <- mean(plac1st)
            
            p_dur11 <- mean(dur1st)
            
            p_dur12 <- mean(dur2st)
            
            # Z-statistics for the two treatment durations
            
            z1 <- (p_dur11-p_plac11)/(sqrt((1/length(plac1st)+1/length(dur1st)))*sigma)
            
            z2 <- (p_dur12-p_plac11)/(sqrt((1/length(plac1st)+1/length(dur2st)))*sigma)
            
            # The null hypothesis is early rejected in the 1st stage
            
            if(s < stage){
              
              
              if( z1 >= func$upperbounds[s] &  z2 >= func$upperbounds[s] & continue == FALSE){
                
                rej_ord_h01andh02[k,j] = TRUE
                
                powerII[k,j] = TRUE
                
                break 
                
              }
              
              
              if( z1 >= func$upperbounds[s] &  z2 < func$lowerbounds[s] & continue == FALSE){
                
                rej_ord_h01noh02[k,j] = TRUE
                
                break
                
              }
              
              if( z2 >= func$upperbounds[s] &  z1 < func$lowerbounds[s] & continue == FALSE){
                
                rej_ord_h02noh01[k,j] = TRUE
                
                break
                
              }
              
              if( z1 <= func$lowerbounds[s] &  z2 <= func$lowerbounds[s] & continue == FALSE){
                
                bothineff[k,j] = TRUE
                
                break
                
              }
              
              # Continue to second stage
              # P(z21 >= u2, z22 >= u2) + P(z21 >= u2, z22 < u2): cell 5 of the decision matrix
              
              if(z1 < func$upperbounds[s] & z1 > func$lowerbounds[s] &
                 z2 < func$upperbounds[s] & z2 > func$lowerbounds[s] & continue == FALSE){
                
                cell5 <- TRUE
                
                continue <- TRUE
              }
              
              
              if(z1 < func$upperbounds[s] & z1 > func$lowerbounds[s]&
                 z2 <= func$lowerbounds[s] & continue == FALSE){
                
                cell2 <- TRUE
                
                continue <- TRUE
              }
              
              if(z1 < func$upperbounds[s] & z1 > func$lowerbounds[s]&
                 z2 >= func$upperbounds[s] & continue == FALSE){
                
                cell3 <- TRUE
                
                continue <- TRUE
              }
              
              
              if(z2 < func$upperbounds[s] & z2 > func$lowerbounds[s] &
                 z1 <= func$lowerbounds[s] & continue == FALSE){
                
                
                cell4 <- TRUE
                
                continue <- TRUE
              }
              
              if(z2 < func$upperbounds[s] & z2 > func$lowerbounds[s] &
                 z1 >= func$upperbounds[s] & continue == FALSE){
                
                
                cell6 <- TRUE
                
                continue <- TRUE
              }
              
              
            }
            
            if(s==stage & continue == FALSE){
              
              if (z1 >= func$upperbounds[s] & z2 >= func$upperbounds[s]){
                
                rej_ord_h01andh02[k,j] = TRUE
                
                powerII[k,j] = TRUE
                
                break
                
                
              }
              
              if (z1 >= func$upperbounds[s] & z2 < func$upperbounds[s]){
                
                rej_ord_h01noh02[k,j] = TRUE
                
                break
                
                
              }
              
              if (z2 >= func$upperbounds[s] & z1 < func$upperbounds[s]){
                
                rej_ord_h02noh01[k,j] = TRUE
                
                break
                
              }
              
              if( z1 < func$upperbounds[s] &  z2 < func$upperbounds[s]){
                
                bothineff[k,j] = TRUE
                
                break
                
              }
            }
            
          }
          
          
          if(cell2==TRUE){
            
            s <- s+1
            
            cell2 <- FALSE
            
            continue <- FALSE
            
            plac21 <- rnorm(n = ceiling(popord*(ratio[s]-ratio[s-1])*rho[3]), mean = mean0, sd = sigma)
            
            dur21 <- rnorm(n = ceiling(popord*(ratio[s]-ratio[s-1])*rho[1]), mean = mean1, sd = sigma)
            
            plac1st <- c(plac21, plac1st)
            
            dur1st <- c(dur21, dur1st)
            
            estimatedss[k,j] <- estimatedss[k,j] +length(plac21)+length(dur21)
            
            # Estimated response rates
            
            p_plac21 <- mean(plac1st)
            
            p_dur21 <- mean(dur1st)
            
            # Z-statistics for the two treatment durations
            
            z1 <- (p_dur21-p_plac21)/(sqrt((1/length(plac1st)+1/length(dur1st)))*sigma)
            
            # Rejection of the null hypothesis at the second stage
            
            if(s < stage){  
              
              if (z1 >= func$upperbounds[s] & continue == FALSE){
                
                rej_ord_h01noh02[k,j] = TRUE
                
                break
                
              }
              if( z1 < func$lowerbounds[s] & continue == FALSE){
                
                bothineff[k,j] = TRUE
                
                break
                
                
              }
              if(z1 < func$upperbounds[s] & z1 > func$lowerbounds[s] & continue == FALSE){
                
                cell2 <- TRUE
                
                continue <- TRUE
              }
              
            }
            
            if (s==stage & continue == FALSE){
              
              if (z1 >= func$upperbounds[s]){
                
                
                rej_ord_h01noh02[k,j] = TRUE
                
                break
                
              }
              if( z1 < func$upperbounds[s]){
                
                bothineff[k,j] = TRUE
                
                break
                
                
              }
              
            }
          }
          
          if(cell3==TRUE){
            
            s <- s+1
            
            cell3 <- FALSE
            
            continue <- FALSE
            
            plac21 <- rnorm(n = ceiling(popord*(ratio[s]-ratio[s-1])*rho[3]), mean = mean0, sd = sigma)
            
            dur21 <- rnorm(n = ceiling(popord*(ratio[s]-ratio[s-1])*rho[1]), mean = mean1, sd = sigma)
            
            plac1st <- c(plac21, plac1st)
            
            dur1st <- c(dur21, dur1st)
            
            estimatedss[k,j] <- estimatedss[k,j] +length(plac21)+length(dur21)
            
            # Estimated response rates
            
            p_plac21 <- mean(plac1st)
            
            p_dur21 <- mean(dur1st)
            
            # Z-statistics for the two treatment durations
            
            z1 <- (p_dur21-p_plac21)/(sqrt((1/length(plac1st)+1/length(dur1st)))*sigma)
            
            # Rejection of the null hypothesis at the second stage
            
            if(s < stage){  
              
              if (z1 >= func$upperbounds[s] & continue == FALSE){
                
                rej_ord_h01andh02[k,j] = TRUE
                
                break
                
              }
              if( z1 < func$lowerbounds[s] & continue == FALSE){
                
                rej_ord_h02noh01[k,j] = TRUE
                
                break
                
                
              }
              if(z1 < func$upperbounds[s] & z1 > func$lowerbounds[s] & continue == FALSE){
                
                cell3 <- TRUE
                
                continue <- TRUE
              }
              
            }
            
            if (s==stage & continue == FALSE){
              
              if (z1 >= func$upperbounds[s]){
                
                
                rej_ord_h01andh02[k,j] = TRUE
                
                break
                
              }
              if( z1 < func$upperbounds[s]){
                
                rej_ord_h02noh01[k,j] = TRUE
                
                break
                
                
              }
              
            }
          }
          
          if(cell4==TRUE){
            
            s <- s+1
            
            cell4 <- FALSE
            
            continue <- FALSE
            
            plac21 <- rnorm(n = ceiling(popord*(ratio[s]-ratio[s-1])*rho[3]), mean = mean0, sd = sigma)
            
            dur22 <- rnorm(n = ceiling(popord*(ratio[s]-ratio[s-1])*rho[2]), mean = mean2, sd = sigma)
            
            plac1st <- c(plac21, plac1st)
            
            dur2st <- c(dur22, dur2st)
            
            estimatedss[k,j] <- estimatedss[k,j] +length(plac21)+length(dur22)
            
            # Estimated response rates
            
            p_plac21 <- mean(plac1st)
            
            p_dur22 <- mean(dur2st)
            
            # Z-statistics for the two treatment durations
            
            z2 <- (p_dur22-p_plac21)/(sqrt((1/length(plac1st)+1/length(dur2st)))*sigma)
            
            # Rejection of the null hypothesis at the second stage
            
            if(s < stage){  
              
              if (z2 >= func$upperbounds[s] & continue == FALSE){
                
                rej_ord_h02noh01[k,j] = TRUE
                
                break
                
              }
              if( z2 < func$lowerbounds[s] & continue == FALSE){
                
                bothineff[k,j] = TRUE
                
                break
                
                
              }
              if(z2 < func$upperbounds[s] & z2 > func$lowerbounds[s] & continue == FALSE){
                
                cell4 <- TRUE
                
                continue <- TRUE
              }
              
            }
            
            if (s==stage & continue == FALSE){
              
              if (z2 >= func$upperbounds[s]){
                
                
                rej_ord_h02noh01[k,j] = TRUE
                
                break
                
              }
              if( z2 < func$upperbounds[s]){
                
                bothineff[k,j] = TRUE
                
                break
                
                
              }
              
            }
          }
          
          if(cell6==TRUE){
            
            s <- s+1
            
            cell6 <- FALSE
            
            continue <- FALSE
            
            plac21 <- rnorm(n = ceiling(popord*(ratio[s]-ratio[s-1])*rho[3]), mean = mean0, sd = sigma)
            
            dur22 <- rnorm(n = ceiling(popord*(ratio[s]-ratio[s-1])*rho[2]), mean = mean2, sd = sigma)
            
            plac1st <- c(plac21, plac1st)
            
            dur2st <- c(dur22, dur2st)
            
            estimatedss[k,j] <- estimatedss[k,j] +length(plac21)+length(dur22)
            
            # Estimated response rates
            
            p_plac21 <- mean(plac1st)
            
            p_dur22 <- mean(dur2st)
            
            # Z-statistics for the two treatment durations
            
            z2 <- (p_dur22-p_plac21)/(sqrt((1/length(plac1st)+1/length(dur2st)))*sigma)
            
            # Rejection of the null hypothesis at the second stage
            
            if(s < stage){  
              
              if (z2 >= func$upperbounds[s] & continue == FALSE){
                
                rej_ord_h01andh02[k,j] = TRUE
                
                break
                
              }
              if( z2 < func$lowerbounds[s] & continue == FALSE){
                
                rej_ord_h01noh02[k,j] = TRUE
                
                break
                
                
              }
              if(z2 < func$upperbounds[s] & z2 > func$lowerbounds[s] & continue == FALSE){
                
                cell6 <- TRUE
                
                continue <- TRUE
              }
              
            }
            
            if (s==stage & continue == FALSE){
              
              if (z2 >= func$upperbounds[s]){
                
                
                rej_ord_h01andh02[k,j] = TRUE
                
                break
                
              }
              if( z2 < func$upperbounds[s]){
                
                rej_ord_h01noh02[k,j] = TRUE
                
                break
                
                
              }
              
            }
          }
          
          
        }
      }
      
    }
    
    # Average probabilities of rej the hypotheses
    
    prop_ord[k,1] <- sum(rej_ord_h01andh02[k,])/nsim
    
    prop_ord[k,2] <- sum(rej_ord_h01noh02[k,])/nsim
    
    prop_ord[k,3] <- sum(powerII[k,])/nsim
    
    prop_ord[k,4] <- sum(bothineff[k,])/nsim
    
    prop_ord[k,5] <- sum(rej_ord_h02noh01[k,])/nsim
    
    estimatedsamplesize[k] <- mean(estimatedss[k,]) 
    
    
  }
  
  colnames(prop_ord) <- c(     "Allprom_rejH01andH02", 
                               "Allprom_rejH01noH02",
                               "power",
                               "Allprom_rejneitherH01andH02",
                               "Allprom_rejH02noH01")
  
  colnames(estimatedsamplesize) <- "ESS Ordered Allprom"
  
  ssarm1stage = popord*rho[1]
  
  ssarm2stage = popord*rho[2]
  
  sscontrolstage = popord
  
  totss = ceiling(ratio[stage]*(ssarm1stage+ssarm2stage+sscontrolstage))
  
  summary_Allprom <- cbind(ssarm1stage,
                           ssarm2stage,
                           sscontrolstage,
                           totss,
                           scenario,
                           prop_ord, 
                           estimatedsamplesize)
  
  colnames(summary_Allprom)[1:4]<- c("patients_arm1stage1 Allprom", 
                                     "patients_arm2stage1 Allprom",
                                     "patients_controlarmstage1 Allprom",
                                     "total sample Allprom")
  
  library(dplyr)
  
  summary_Allprom <- summary_Allprom %>% mutate(
    sumprob = Allprom_rejH01andH02+Allprom_rejH01noH02+Allprom_rejneitherH01andH02+Allprom_rejH02noH01,
    typeIprob = Allprom_rejH01andH02+Allprom_rejH01noH02+Allprom_rejH02noH01
  ) %>% select(
    -power
  )
  
  return(summary_Allprom)
  
}











####### Functions ORD ####

##### Function to find critical bounds ##########

#' @u1: parameter for the grid search
#' @interimN: number of interim analyses
#' @alpha: alpha-level of the hypothesis test
#' @cov: covariance matrix
#' @first: covariance matrix at the first stage
#' @r: allocation ratio sample size in the first and second stage
#' @prec: precision for the boundaries: numbers after the comma
#' @arms: number of arms (3)
#' @ushape1 - @ushape2: upper boundary shape for first and second treatment: "pocock", "obf", "triangular"
#' @lshape1 - @lshape2: lower boundary shape for first and second treatment: "pocock", "obf", "triangular"

bounds_3armsORD_differentbounds <- function(u1, interimN, arms,r,
                                            alpha, cov, first, ushape1,ushape2,lshape1,lshape2,prec){
  
  library(mvtnorm)
  library(gtools)
  library(cubature)
  
  alphat <- rep(0,length = interimN)
  
  for(t in 1:length(u1)){
    
    for (i in 1:interimN){
      
      if (ushape1 == "obf") {
        u11 <- (u1[t] * 1/sqrt(r))[1:i]
      }
      else if (ushape1 == "pocock") {
        u11 <- rep(u1[t], i)
      }
      else if (ushape1 == "triangular") {
        u11 <- (u1[t] * (1 + r/max(r))/sqrt(r))[1:i]
      }
      
      if (lshape1 == "obf") {
        l11 <- c(-u1[t] * 1/sqrt(r))[1:(i-1)]
      }
      else if (lshape1 == "pocock") {
        l11 <- c(rep(-u1[t], i-1))
      }
      else if (lshape1 == "triangular") {
        if (ushape1 == "triangular") {
          l11 <- c(-u1[t] * ((1 -3* r/max(r))/sqrt(r)))[1:(i-1)]
        }
      }
      
      if (ushape2 == "obf") {
        u2 <- (u1[t] * 1/sqrt(r))[1:i]
      }
      else if (ushape2 == "pocock") {
        u2 <- rep(u1[t], i)
      }
      else if (ushape2 == "triangular") {
        u2 <- (u1[t] * (1 + r/max(r))/sqrt(r))[1:i]
      }
      
      if (lshape2 == "obf") {
        l2 <- c(-u1[t] * 1/sqrt(r))[1:(i-1)]
      }
      else if (lshape2 == "pocock") {
        l2 <- c(rep(-u1[t], i-1))
      }
      else if (lshape2 == "triangular") {
        if (ushape2 == "triangular") {
          l2 <- c(-u1[t] * ((1 -3* r/max(r))/sqrt(r)))[1:(i-1)]
        }
      }
      
      if(i ==1){
        set.seed(123456)
        
        alphat[i] <- pmvnorm(lower = u11[1],
                             upper = Inf,
                             sigma = first[1])[1]
        
      }
      
      
      if(i>1){
        
        alphat[i] <- pmvnorm(lower = c(-Inf,u2[1], u11[2],-Inf),
                             upper = c(u11[1],Inf,Inf,Inf),
                             sigma = cov)[1]+
          pmvnorm(lower = c(l11[1],l2[1], u11[2],-Inf),
                  upper = c(u11[1],u2[1],Inf,Inf),
                  sigma = cov)[1]+
          pmvnorm(lower = c(l11[1],-Inf, u11[2],-Inf),
                  upper = c(u11[1],l2[1],Inf,Inf),
                  sigma = cov)[1]
        
        
      }
      
    }
    
    
    if(sum(alphat) <= alpha){
      
      upperbound1 <- u11
      lowerbound1 <- l11
      upperbound2 <- u2
      lowerbound2 <- l2
      totalpha <- sum(alphat)
      alow <- u1[t]
      finalalpha <- alphat
      
    }
    
    if(sum(alphat) > alpha){
      
      aup <- u1[t]
      break
      
      
    }
    
  }
  
  results <- list(upperbound1, upperbound2, alow, aup, lowerbound1,lowerbound2,totalpha)
  
  names(results) <- c("upperbound1", "upperbound2", "alow", "aup", "lowerbound1", "lowerbound2","totalpha")
  return(results)
  
}

##### Function to find sample size for 3-arm 2-stage ORD ##########

#' @theta: vector of clinically relevant difference 
#' @sigma: standard deviation of the population
#' @interimN: number of stages
#' @alpha: alpha-level of the hypothesis test
#' @beta: beta-level of the test
#' @r: allocation ratio sample size in the first and second stage
#' @rho: allocation ratio sample size in each treatment divided by sample size in the control
#' @prec: precision for the boundaries: numbers after the comma
#' @arms: number of arms (3)
#' @ushape1 - @ushape2: upper boundary shape for first and second treatment: "pocock", "obf", "triangular"
#' @lshape1 - @lshape2: lower boundary shape for first and second treatment: "pocock", "obf", "triangular"
#' @power: "reject all" or "reject at least one"

boundariesSample_3armsORD_differentbounds <- function(theta,
                                                      sigma, 
                                                      interimN,
                                                      alpha, 
                                                      beta,
                                                      r, 
                                                      rho,
                                                      prec,
                                                      arms,
                                                      ushape1,
                                                      ushape2,
                                                      lshape1,
                                                      lshape2,
                                                      power
){ 
  
  library(mvtnorm)
  library(cubature)
  
  #construction of the covariance matrix
  
  R <- rho[1:(length(rho)-1)]/(rho[1:(length(rho)-1)]+rho[length(rho)])
  
  bottom<-matrix(r,length(r),length(r))
  top<-matrix(rep(r,rep(length(r),length(r))),length(r),length(r))
  top[upper.tri(top)]<-t(top)[upper.tri(top)]
  bottom[upper.tri(bottom)]<-t(bottom)[upper.tri(bottom)] 
  
  cov <-sqrt(top/bottom)
  
  ma <- matrix(0, nrow = length(R), ncol = length(R))
  
  for (i in 1:length(R)){
    
    for(j in 2:length(R)){
      
      ma[i,j] <- sqrt(R[i]*R[j])
    }
  }
  
  diag(ma) <- 1
  
  ma[lower.tri(ma)]  <- t(ma)[lower.tri(ma)]
  
  block <- list()
  
  dim <- 1
  
  for (i in 1:nrow(cov)){
    
    for(j in 1:ncol(cov)){
      
      block[[dim]] <- cov[i,j]*ma
      
      dim <- dim+1
    }
    
  }
  
  covb <- NULL
  
  for (i in 1:length(block)){
    
    covb <- cbind(covb, block[[i]])
  }
  
  int <- seq(0, ncol(covb), by=((arms-1)*length(r)))
  
  covbfin <- NULL
  s <- 1
  for(j in 2:length(int)){
    covbfin <- rbind(covbfin, covb[,s:int[j]])
    s <- int[j]+1
  }
  
  diag(covbfin) <- 1
  
  #covariance matrix at the first stage
  
  firststagematrix <- covbfin[1:(arms-1),1:(arms-1)]
  
  diag(firststagematrix) <- 1
  
  if(interimN == 1){
    
    u1 <- qnorm(1-alpha)
    
    pop <- seq(1,2000, by = 1)
    
    maxsample <- NULL
    
    betat <- NULL
    
    n <- 1
    
    bound <- matrix(NA, nrow = 2, ncol = arms-1)
    
    bound[1,] <- rep(u1, arms-1)
    
    bound[2,] <- rep(Inf, arms-1)
    
    rarm <- c(rep(r, each=arms-1))
    
    for (n in 1:length(pop)){
      
      prob <- 0
      
      mean1 <- ((theta)*(sqrt(rarm*R*pop[n]))/(sigma))[1:2]
      
      if(power == "reject all"){
        
        prob <- pmvnorm(lower = bound[1,],
                        upper = bound[2,],
                        sigma = covbfin[1:(arms-1), 1:(arms-1)],
                        mean = mean1)[1]
      }
      
      if(power == "reject at least one"){
        
        prob <- pmvnorm(lower = u1,
                        upper = Inf,
                        sigma = covbfin[1, 1],
                        mean = mean1[1])[1]
      }
      
      if(prob > 1-beta){
        
        maxsample <- pop[n]
        break
        
      }
      
      
    }
    results <- list(u1,-u1,maxsample)
    
    names(results) <- c("upperbounds1", "lowerbounds1",
                        "sample size per arm per stage"
    )
    
    
  }
  
  else{
    
    ## Research of the critical bounds
    
    u1 <- seq(from = 7, to = 0, by = -1)
    
    first <- bounds_3armsORD_differentbounds(u1, interimN, arms,r, 
                                             alpha, covbfin, firststagematrix,
                                             ushape1, ushape2,lshape1, lshape2,prec)
    
    low <- first$alow
    
    up <- first$aup
    
    p <- rep(1, times = prec)
    
    for (p in 1:length(p)){
      
      callf <- bounds_3armsORD_differentbounds(seq(from = low, 
                                                   to = up, by = -1/(10^p)), 
                                               interimN= interimN,
                                               arms,
                                               r,
                                               alpha,
                                               covbfin,
                                               firststagematrix,
                                               ushape1, 
                                               ushape2,
                                               lshape1,
                                               lshape2,
                                               prec)
      
      low <- callf$alow
      
      up <- callf$aup
    }
    
    
    upperbound1 <- callf$upperbound1
    upperbound2 <- callf$upperbound2
    lowerbound1 <- callf$lowerbound1
    lowerbound2 <- callf$lowerbound2
    a <- callf$alow
    totalalpha <- callf$totalpha
    
    
    maxsample <- NULL
    
    betat <- NULL
    
    pop <- seq(from = 1, to = 500, by = 1)
    
    if (length(theta)!=(arms-1)){
      print("Error: the length of theta should be equal to arms-1")
    }
    
    rarm <- c(rep(r, each=arms-1))
    
    for (n in 1:length(pop)){
      
      mean2 <- (theta)*(sqrt(rarm*R*pop[n]))/(sigma)
      
      mean1 <- mean2[1:2] #at the first stage
      
      for (i in 1:interimN){
        
        u11 <- upperbound1[1:i]
        l11 <- lowerbound1[1:(i-1)]
        
        u2 <- upperbound2[1:i]
        l2 <- lowerbound2[1:(i-1)]
        
        if(power == "reject all"){
          
          if(i ==1){
            set.seed(123456)
            
            betat[i] <- pmvnorm(lower = c(u11[1],u2[1]),
                                upper = c(Inf, Inf),
                                sigma = firststagematrix,
                                mean = mean1)[1]
            
          }
          if(i>1){
            
            set.seed(123456)
            betat[i] <- pmvnorm(lower = c(l11[1],u2[1], u11[2],u2[2]),
                                upper = c(u11[1],Inf,Inf,Inf),
                                sigma = covbfin,
                                mean = mean2)[1]+
              pmvnorm(lower = c(-Inf,u2[1], u11[2],u2[2]),
                      upper = c(l11[1],Inf,Inf,Inf),
                      sigma = covbfin,
                      mean = mean2)[1]+
              pmvnorm(lower = c(u11[1],l2[1], -Inf,u2[2]),
                      upper = c(Inf,u2[1],Inf,Inf),
                      sigma = covbfin,
                      mean = mean2)[1]+
              pmvnorm(lower = c(l11[1],l2[1], u11[2],u2[2]),
                      upper = c(u11[1],u2[1],Inf,Inf),
                      sigma = covbfin,
                      mean = mean2)[1]
            
          }
        }
        if(power == "reject at least one"){
          
          if(i ==1){
            set.seed(123456)
            
            betat[i] <- pmvnorm(lower = u11[1],
                                upper = Inf,
                                sigma = firststagematrix[1,1],
                                mean = mean1[1])[1]
            
          }
          
          
          if(i>1){
            
            set.seed(123456)
            betat[i] <- pmvnorm(lower = c(l11[1],-Inf, u2[2],-Inf),
                                upper = c(u11[1],Inf,Inf,Inf),
                                sigma = covbfin,
                                mean = mean2)[1]+
              pmvnorm(lower = c(-Inf,u2[1], u11[2],-Inf),
                      upper = c(l11[1],Inf,Inf,Inf),
                      sigma = covbfin,
                      mean = mean2)[1]
            
            
          }
          
          
          
        }
        
      }
      
      if(sum(betat) > 1-beta){
        
        maxsample <- pop[n]
        break
        
      }
      
    }
    
    results <- list(totalalpha,
                    lowerbound1,upperbound1,lowerbound2,upperbound2,maxsample )
    
    names(results) <- c("sum_alphat", 
                        "lowerbounds1", "upperbounds1", "lowerbounds2", "upperbounds2","sample size per arm per stage")
  }
  
  return(results)
}






##### Function to simulate 3-arm J-stage ORD ###########

simul <- function(scenario, nsim, stage, popord, ratio, rho, sigma, seed, fun ){
  
  # save rejection first and/or second hypothesis in each stage and
  # save rejection for the power computation when theta = (theta1, theta2) and
  # theta1 >= theta2 >= delta0 > 0
  
  bothineff <- matrix(data=FALSE, nrow=nrow(scenario),
                      ncol=nsim)
  
  rej_ord_h01andh02 <- matrix(data=FALSE, nrow=nrow(scenario),
                              ncol=nsim)
  
  rej_ord_h01noh02 <- matrix(data=FALSE, nrow=nrow(scenario),
                             ncol=nsim)
  
  rej_ord_h02noh01 <- matrix(data=FALSE, nrow=nrow(scenario),
                             ncol=nsim)
  
  prop_ord <- matrix(data=NA, ncol=4,
                     nrow = nrow(scenario))
  
  # Estimated sample size for each scenario. 
  
  estimatedss <- matrix(data=NA, nrow=nrow(scenario),
                        ncol=nsim)
  
  estimatedsamplesize <- matrix(data=NA, ncol=1,
                                nrow = nrow(scenario))
  
  for (k in 1:nrow(scenario)){
    
    set.seed(seed = seed)
    
    mean0 <- scenario[k,1]
    
    mean1 <- scenario[k,2]
    
    mean2 <- scenario[k,3]
    
    for (j in 1:nsim){
      
      continue <- TRUE
      
      plac1st <- NULL
      
      dur1st <- NULL
      
      dur2st <- NULL
      
      cell4 <- FALSE
      
      cell8 <- FALSE
      
      for(s in 1:stage){
        
        if (continue==TRUE){
          
          if(s==1|| (cell5==TRUE || cell2==TRUE)){
            
            continue <- FALSE
            
            cell2 <- FALSE
            
            cell5 <- FALSE
            
            if(s==1){
              
              plac1 <- rnorm(n = popord*(ratio[s]*rho[3]), mean = mean0, sd = sigma)
              
              dur1 <- rnorm(n = popord*(ratio[s]*rho[1]), mean = mean1, sd = sigma)
              
              dur2 <- rnorm(n = popord*(ratio[s]*rho[2]), mean = mean2, sd = sigma)
            }
            
            if(s>1){
              
              plac1 <- rnorm(n = ceiling(popord*(ratio[s]-ratio[s-1])*rho[3]), mean = mean0, sd = sigma)
              
              dur1 <- rnorm(n = ceiling(popord*(ratio[s]-ratio[s-1])*rho[1]), mean = mean1, sd = sigma)
              
              dur2 <- rnorm(n = ceiling(popord*(ratio[s]-ratio[s-1])*rho[2]), mean = mean2, sd = sigma)
            }
            
            # Estimated response rates on the observed data
            
            plac1st <- c(plac1, plac1st)
            
            dur1st <- c(dur1, dur1st)
            
            dur2st <- c(dur2, dur2st)
            
            estimatedss[k,j] <- length(plac1st)+length(dur1st)+length(dur2st)
            
            p_plac11 <- mean(plac1st)
            
            p_dur11 <- mean(dur1st)
            
            p_dur12 <- mean(dur2st)
            
            # Z-statistics for the two treatment durations
            
            z1 <- (p_dur11-p_plac11)/(sqrt((1/length(plac1st)+1/length(dur1st)))*sigma)
            
            z2 <- (p_dur12-p_plac11)/(sqrt((1/length(plac1st)+1/length(dur2st)))*sigma)
            
            # The null hypothesis is early rejected in the 1st stage
            
            if(s < stage){
              
              
              if( z1 >= fun$upperbounds1[s] &  z2 >= fun$upperbounds2[s] & continue == FALSE){
                
                rej_ord_h01andh02[k,j] = TRUE
                
                break 
                
              }
              
              
              if( z1 >= fun$upperbounds1[s] &  z2 < fun$lowerbounds2[s] & continue == FALSE){
                
                rej_ord_h01noh02[k,j] = TRUE
                
                break
                
              }
              
              if( z1 <= fun$lowerbounds1[s] &  z2 < fun$upperbounds2[s] & continue == FALSE){
                
                bothineff[k,j] = TRUE
                
                break
                
              }
              
              # Continue to second stage
              # P(z21 >= u2, z22 >= u2) + P(z21 >= u2, z22 < u2): cell 5 of the decision matrix
              
              if(z1 < fun$upperbounds1[s] & z1 > fun$lowerbounds1[s] &
                 z2 < fun$upperbounds2[s] & z2 > fun$lowerbounds2[s] & continue == FALSE){
                
                cell5 <- TRUE
                
                continue <- TRUE
              }
              
              # P(z21 >= u2, z22 >= u2) + P(z21 >= u2, z22 < u2): cell 2 of the decision matrix
              
              if(z1 < fun$upperbounds1[s] &
                 z2 >=fun$upperbounds2[s] & continue == FALSE){
                
                cell2 <- TRUE
                
                continue <- TRUE
              }
              
              # P(z21 >= u2): cell 4 of the decision matrix
              
              if(z1 >= fun$upperbounds1[s] &
                 z2 < fun$upperbounds2[s] & z2 > fun$lowerbounds2[s] & continue == FALSE){
                
                
                cell4 <- TRUE
                
                continue <- TRUE
              }
              
              
              # P(z21 >= u2): cell 8 of the decision matrix
              
              if(z2 <= fun$lowerbounds2[s] &
                 z1 < fun$upperbounds1[s] & z1 > fun$lowerbounds1[s] & continue == FALSE){
                
                cell8 <- TRUE
                
                continue <- TRUE
                
              }
            }
            
            if(s==stage & continue == FALSE){
              
              if (z1 >= fun$upperbounds1[s] & z2 >= fun$upperbounds2[s]){
                
                rej_ord_h01andh02[k,j] = TRUE
                
                break
                
                
              }
              
              if (z1 >= fun$upperbounds1[s] & z2 < fun$upperbounds2[s]){
                
                rej_ord_h01noh02[k,j] = TRUE
                
                break
                
                
              }
              
              if (z2 >= fun$upperbounds2[s] & z1 < fun$upperbounds1[s]){
                
                rej_ord_h02noh01[k,j] = TRUE
                bothineff[k,j] = TRUE
                
                break
                
              }
              
              if( z1 < fun$upperbounds1[s] &  z2 < fun$upperbounds2[s]){
                
                bothineff[k,j] = TRUE
                
                break
                
              }
            }
            
          }
          
          if(cell4==TRUE){
            
            s <- s+1
            
            cell4 <- FALSE
            
            continue <- FALSE
            
            plac21 <- rnorm(n = ceiling(popord*(ratio[s]-ratio[s-1])*rho[3]), mean = mean0, sd = sigma)
            
            dur22 <- rnorm(n = ceiling(popord*(ratio[s]-ratio[s-1])*rho[2]), mean = mean2, sd = sigma)
            
            # Estimated response rates on the observed data
            
            plac1st <- c(plac21, plac1st)
            
            dur2st <- c(dur22, dur2st)
            
            estimatedss[k,j] <- estimatedss[k,j]  + length(plac21)+length(dur22)
            
            p_plac21 <- mean(plac1st)
            
            p_dur22 <- mean(dur2st)
            
            # Z-statistics for the two treatment durations
            
            z2 <- (p_dur22-p_plac21)/(sqrt((1/length(plac1st)+1/length(dur2st)))*sigma)
            
            if(s < stage){
              
              
              if (z2 >= fun$upperbounds2[s] & continue == FALSE){
                
                rej_ord_h01andh02[k,j] = TRUE
                
                break
                
              }
              
              if (z2<= fun$lowerbounds2[s] & continue == FALSE){
                
                rej_ord_h01noh02[k,j] = TRUE
                
                break
                
              }
              
              if(z2 < fun$upperbounds2[s] & z2 > fun$lowerbounds2[s] & continue == FALSE){
                
                cell4 <- TRUE
                
                continue <- TRUE
                
              }
              
            }
            
            if (s==stage & continue == FALSE){
              
              if (z2 >= fun$upperbounds2[s]){
                
                
                rej_ord_h01andh02[k,j] = TRUE
                
                break
                
              }
              if( z2 < fun$upperbounds2[s]){
                
                rej_ord_h01noh02[k,j] = TRUE
                
                break
                
                
              }
              
            }
          }
          
          if(cell8==TRUE){
            
            s <- s+1
            
            cell8 <- FALSE
            
            continue <- FALSE
            
            plac21 <- rnorm(n = ceiling(popord*(ratio[s]-ratio[s-1])*rho[3]), mean = mean0, sd = sigma)
            
            dur21 <- rnorm(n = ceiling(popord*(ratio[s]-ratio[s-1])*rho[1]), mean = mean1, sd = sigma)
            
            plac1st <- c(plac21, plac1st)
            
            dur1st <- c(dur21, dur1st)
            
            estimatedss[k,j] <- estimatedss[k,j] +length(plac21)+length(dur21)
            
            # Estimated response rates
            
            p_plac21 <- mean(plac1st)
            
            p_dur21 <- mean(dur1st)
            
            # Z-statistics for the two treatment durations
            
            z1 <- (p_dur21-p_plac21)/(sqrt((1/length(plac1st)+1/length(dur1st)))*sigma)
            
            # Rejection of the null hypothesis at the second stage
            
            if(s < stage){  
              
              if (z1 >= fun$upperbounds1[s] & continue == FALSE){
                
                rej_ord_h01noh02[k,j] = TRUE
                
                break
                
              }
              if( z1 < fun$lowerbounds1[s] & continue == FALSE){
                
                bothineff[k,j] = TRUE
                
                break
                
                
              }
              if(z1 < fun$upperbounds1[s] & z1 > fun$lowerbounds1[s] & continue == FALSE){
                
                cell8 <- TRUE
                
                continue <- TRUE
              }
              
            }
            
            if (s==stage & continue == FALSE){
              
              if (z1 >= fun$upperbounds1[s]){
                
                
                rej_ord_h01noh02[k,j] = TRUE
                
                break
                
              }
              if( z1 < fun$upperbounds1[s]){
                
                bothineff[k,j] = TRUE
                
                break
                
                
              }
              
            }
          }
          
          
        }
      }
    }
    
    # Average probabilities of rej the hypotheses
    
    prop_ord[k,1] <- sum(rej_ord_h01andh02[k,])/nsim
    
    prop_ord[k,2] <- sum(rej_ord_h01noh02[k,])/nsim
    
    prop_ord[k,3] <- sum(bothineff[k,])/nsim
    
    prop_ord[k,4] <- sum(rej_ord_h02noh01[k,])/nsim
    
    estimatedsamplesize[k] <- mean(estimatedss[k,]) 
    
  }
  
  colnames(prop_ord) <- c(     "ORD_rejH01andH02", 
                               "ORD_rejH01noH02",
                               "ORD_rejneitherH01andH02",
                               "ORD_rejH02noH01")
  
  colnames(estimatedsamplesize) <- "ESS Ordered ORD"
  
  ssarm1stage = popord*rho[1]
  
  ssarm2stage = popord*rho[2]
  
  sscontrolstage = popord
  
  totss = ceiling(ratio[stage]*(ssarm1stage+ssarm2stage+sscontrolstage))
  
  summary_ORD  <- cbind(ssarm1stage,
                        ssarm2stage,
                        sscontrolstage,
                        totss,
                        scenario,
                        prop_ord, 
                        estimatedsamplesize)
  
  colnames( summary_ORD )[1:4]<- c("patients_arm1stage1 ORD", 
                                   "patients_arm2stage1 ORD",
                                   "patients_controlarmstage1 ORD",
                                   "total sample ORD")
  
  library(dplyr)
  
  summary_ORD <- summary_ORD %>% mutate(
    sumprob = ORD_rejH01andH02+ORD_rejH01noH02+ORD_rejneitherH01andH02,
    typeIprob = ORD_rejH01andH02+ORD_rejH01noH02
  ) %>% select(
    -ORD_rejH02noH01
  )
  
  return(summary_ORD)
  
}












####### Functions Urach & Posch ####
##### Function to find critical bounds & SS ####

#' @theta: clinically relevant difference 
#' @sigma: standard deviation of the population
#' @interimN: number of interim analyses
#' @alpha: alpha-level of the hypothesis test
#' @beta: beta-level of the test
#' @r: ratio of population in the stages (vector)
#' @prec: precision term for the boundaries, number of values after the comma
#' @ushape: "pocock", "obf", "triangular"
#' @lshape: "pocock", "obf", "triangular"

boundpowerUrach3arm2stage <- function(theta, 
                                      sigma, 
                                      interimN, 
                                      alpha, 
                                      beta, 
                                      r,
                                      prec,
                                      ushape, 
                                      lshape){ 
  
  library(mvtnorm)
  library(MAMS)
  
  betat <- NULL
  
  maxsample <- NULL
  
  covbfinB <- matrix(c(1, 1/2, sqrt(1/2),1/2*sqrt(1/2), 
                       1/2, 1, 1/2*sqrt(1/2), sqrt(1/2),
                       sqrt(1/2),  1/2*sqrt(1/2), 1, 1/2,
                       1/2*sqrt(1/2), sqrt(1/2), 1/2,1),
                     nrow = 4,
                     ncol = 4,
                     byrow = TRUE)
  
  firststagematrixB <- covbfinB[1:2,1:2]
  
  global <- mams(K=2,J=interimN, 
                 ushape = ushape, 
                 lshape = lshape, 
                 nstart = 1,
                 delta = theta[1],
                 r = c(1,2),
                 r0 = c(1,2),
                 p = NULL,
                 p0 = NULL,
                 delta0 = 0,
                 sd = sigma,
                 power = 1-beta,
                 alpha=alpha)
  
  u1 <- round(global$u[1],prec)
  u2 <- round(global$u[2],prec)
  l1 <- round(global$l[1],prec)
  
  single <- mams(K=1,J=interimN, 
                 ushape = ushape, 
                 lshape = lshape, 
                 nstart = 1,
                 delta = theta[1],
                 r = c(1,2),
                 r0 = c(1,2),
                 p = NULL,
                 p0 = NULL,
                 delta0 = 0,
                 sd = sigma,
                 power = 1-beta,
                 alpha=alpha)
  
  v1 <- round(single$u[1],prec)
  v2 <- round(single$u[2],prec)
  vl1 <- round(single$l[1],prec)
  
  pop <- seq(from = 1, to = 1500, by = 1)
  
  for (n in 1:length(pop)){
    
    mean1 <- rep(theta*(sqrt(2*r*rep(pop[n], length(r))))/(2*sigma), each = 2)
    
    for (i in 1:interimN){
      
      if(i ==1){
        set.seed(123456)
        
        betat[i] <- pmvnorm(lower = c(u1,u1),
                            upper = c(Inf,Inf),
                            sigma = firststagematrixB,
                            mean = mean1[1:2])[1]+
          pmvnorm(lower = c(v1,u1),
                  upper = c(u1,Inf),
                  sigma = firststagematrixB,
                  mean = mean1[1:2])[1]+
          pmvnorm(lower = c(u1,v1),
                  upper = c(Inf,u1),
                  sigma = firststagematrixB,
                  mean = mean1[1:2])[1]
        
        
      }
      else{
        
        set.seed(123456)
        
        betat[i] <- pmvnorm(lower = c(l1, l1,u2, u2),
                            upper = c(u1,u1, Inf, Inf),
                            sigma = covbfinB,
                            mean = mean1)[1]+
          pmvnorm(lower = c(l1, l1,v2, u2),
                  upper = c(u1,u1, u2, Inf),
                  sigma = covbfinB,
                  mean = mean1)[1]+
          pmvnorm(lower = c(l1, l1,u2, v2),
                  upper = c(u1,u1, Inf, u2),
                  sigma = covbfinB,
                  mean = mean1)[1]+
          pmvnorm(lower = c(u1, l1,-Inf, v2),
                  upper = c(Inf,v1, Inf, Inf),
                  sigma = covbfinB,
                  mean = mean1)[1]+
          pmvnorm(lower = c(l1,u1,v2,-Inf),
                  upper = c(v1,Inf, Inf, Inf),
                  sigma = covbfinB,
                  mean = mean1)[1]
        
      }
    }
    
    if(sum(betat) > 1-beta){
      maxsample <- pop[n]
      break
      
    }
  }
  
  results <- list(
    u1,u2,l1,v1,v2,vl1, maxsample)
  
  names(results) <- c(
    "u1","u2","l1","v1","v2","vl1","sample size per arm per stage")
  
  return(results)
}



####### Functions for Bayesian design ####
##### Function to find critical bounds & SS ####

# Function to compute the FWER under the global and partial nulls and the power for
# 3-arm 2-stage Bayesian design

#' @x: vector contianing sample sizes per arm per stage and critical bounds
#' @mu01: mean of the prior distribution for mu^(1)
#' @mu00: mean of the prior distribution for mu^(0)
#' @tau: precision fro patients' outcomes
#' @tau00: prior precision on distribution for mu^(0)
#' @tau01: prior precision on distribution for mu^(1)
#' @meand: mean of the prior distribution for delta^(1)
#' @taud: prior precision on distribution for delta^(1)
#' @alpha: alpha-level for FWER control
#' @beta: beta-level for power requirement
#' @null: global null configuration
#' @power: power configuration
#' @power2: partial null configuration with theta^(2)=0
#' @power3: partial null configuration with theta^(1)=0

sim <- function(x, mu01, mu00, tau,
                tau00, meand, tau01,taud,
                alpha, beta, null,power,power2,power3){
  
  
  ################# null #####################
  
  tau02 <- 1/(1/taud+1/tau01)
  
  mu02 <- mu01-meand
  
  var_delta <- 1/taud
  
  ##### First stage 
  
  omega <- matrix(c(1/tau01, 1/tau01,
                    1/tau01, 1/tau01+var_delta), nrow=2, byrow = TRUE)
  
  sigma <- matrix(c(1/tau, 0,
                    0, 1/tau), nrow=2, byrow = TRUE)
  
  mu <- c(mu01, mu02)
  
  invsigma <- solve(sigma)
  
  invomega <- solve(omega)
  
  A2 <- invomega+c(x[1],x[2])*invsigma
  
  invA2 <- solve(A2)
  
  coeff1y1 <- invA2[1,1]*(invsigma*x[1])[1,1] 
  
  coeff1y2 <- invA2[1,2]*(invsigma*x[2])[2,2] 
  
  coeffy0 <- x[3]*tau/(tau00+x[3]*tau) 
  
  const0 <- mu00*tau00/(tau00+x[3]*tau)
  
  const1 <- -(invA2%*%invomega%*%(mu))[1,1]+const0 
  
  const2 <- -(invA2%*%invomega%*%(mu))[2,1]+const0 
  
  coeff2y1 <- invA2[2,1]*(invsigma*x[1])[1,1] 
  
  coeff2y2 <- invA2[2,2]*(invsigma*x[2])[2,2] 
  
  post_varmu0 <- 1/(tau00+x[3]*tau)
  
  post_sdtheta1 <- sqrt(invA2[1,1]+post_varmu0) 
  
  post_sdtheta2 <- sqrt(invA2[2,2]+post_varmu0) 
  
  # means and standard deviation for multivariate normal 
  
  mean1 <- sum(c(coeff1y1,coeff1y2,-coeffy0)*null)
  
  sd1 <- sqrt(sum(c(coeff1y1,coeff1y2)^2)/(tau*x[1])+coeffy0^2/(tau*x[3]))
  
  mean2 <-  sum(c(coeff2y1,coeff2y2,-coeffy0)*null)
  
  sd2 <- sqrt(sum(c(coeff2y1,coeff2y2)^2)/(tau*x[1])+coeffy0^2/(tau*x[3]))
  
  # covariance matrix 
  
  cov12 <- (coeff1y1*coeff2y1+coeff1y2*coeff2y2)/(tau*x[1])+(coeffy0^2)/(tau*x[3])
  
  cov_matr12 <- matrix(c(sd1^2, cov12, 
                         cov12, sd2^2), byrow = TRUE, nrow =2, ncol=2)
  
  # thresholds
  
  thr1u <- const1+qnorm(x[4])*post_sdtheta1
  thr1l <- const1+qnorm(x[5])*post_sdtheta1
  
  thr2u <- const2+qnorm(x[4])*post_sdtheta2
  thr2l <- const2+qnorm(x[5])*post_sdtheta2
  
  # integrals
  
  int1null <-  1-pnorm(thr1u, sd=sd1,mean =mean1)
  
  int2null <-  1-pnorm(thr2u, sd=sd2,mean =mean2)
  
  int12null <- pmvnorm(lower = c(thr1u, thr2u),
                       upper = c(Inf, Inf),
                       sigma = cov_matr12,
                       mean = c(mean1, mean2),
                       algorithm = GenzBretz(abseps = 1*10^-20))[1]
  
  ####### Continue to second stage with one arm
  
  #continue only with arm 1 (coeff1-..)
  
  A21_2st <-invomega+c(2*x[1],x[2])*invsigma
  
  invA21_2st <- solve(A21_2st)
  
  coeff1y1_2st <- invA21_2st[1,1]*(invsigma*2*x[1])[1,1] 
  
  coeff1y2_2st <- invA21_2st[1,2]*(invsigma*x[2])[2,2] 
  
  coeffy0_2st <- 2*x[3]*tau/(tau00+2*x[3]*tau) 
  
  const0_2st <- mu00*tau00/(tau00+2*x[3]*tau)
  
  const1_2st <- -(invA21_2st%*%invomega%*%(mu))[1,1]+const0_2st
  
  post_varmu0_2st <- 1/(tau00+2*x[3]*tau) #posterior variance of mu0
  
  post_sdtheta1_2st <- sqrt(invA21_2st[1,1]+post_varmu0_2st) #posterior standard deviation of theta1
  
  # means and standard deviation for multivariate normal 
  
  mean11_2st <- sum(c(coeff1y1_2st,coeff1y2_2st,-coeffy0_2st)*null)
  
  sd11_2st <- sqrt(coeff1y1_2st^2/(tau*2*x[1])+coeff1y2_2st^2/(tau*x[2])+coeffy0_2st^2/(tau*2*x[3]))
  
  # covariance matrix 
  
  cov13 <- (coeff1y1*coeff1y1_2st/2+coeff1y2*coeff1y2_2st)/(tau*x[1])+(coeffy0*coeffy0_2st/2)/(tau*x[3])
  
  cov23 <- (coeff2y1*coeff1y1_2st/2+coeff2y2*coeff1y2_2st)/(tau*x[1])+(coeffy0*coeffy0_2st/2)/(tau*x[3])
  
  cov_matr112 <- matrix(c(sd1^2, cov12, cov13,
                          cov12, sd2^2, cov23,
                          cov13, cov23, sd11_2st^2), byrow = TRUE, nrow =3, ncol=3)
  
  # thresholds
  
  thr1u1_2st <- const1+qnorm(x[4])*post_sdtheta1
  thr1l1_2st <- const1+qnorm(x[5])*post_sdtheta1
  
  thr2u1_2st <- const2+qnorm(x[5])*post_sdtheta2
  thr2l1_2st <- -Inf
  
  thr3l1_2st <- const1_2st+qnorm(x[6])*post_sdtheta1_2st
  thr3u1_2st <- Inf
  
  # integrals
  
  int_1_2st <- pmvnorm(lower = c(thr1l1_2st,thr2l1_2st,thr3l1_2st),
                       upper = c(thr1u1_2st, thr2u1_2st,thr3u1_2st),
                       sigma = cov_matr112,
                       mean = c(mean1, mean2, mean11_2st),
                       algorithm = GenzBretz(abseps = 1*10^-20))[1]
  
  #continue only with arm 2 (coeff2-..)
  
  A22_2st <- invomega+c(x[1],2*x[2])*invsigma
  
  invA22_2st <- solve(A22_2st)
  
  const2_2st <- -(invA22_2st%*%invomega%*%(mu))[2,1]+const0_2st 
  
  coeff2y1_2st <- invA22_2st[2,1]*(invsigma*x[1])[1,1] 
  
  coeff2y2_2st <- invA22_2st[2,2]*(invsigma*2*x[2])[2,2] 
  
  post_varmu0_2st <- 1/(tau00+2*x[3]*tau) #posterior variance of mu0
  
  post_sdtheta2_2st <- sqrt(invA22_2st[2,2]+post_varmu0_2st) #posterior standard deviation of theta2
  
  # means and standard deviation for multivariate normal 
  
  mean22_2st <- sum(c(coeff2y1_2st,coeff2y2_2st,-coeffy0_2st)*null)
  
  sd22_2st <- sqrt(coeff2y1_2st^2/(tau*x[1])+coeff2y2_2st^2/(tau*2*x[2])+coeffy0_2st^2/(tau*2*x[3]))
  
  # covariance matrix 
  
  cov13_2 <- (coeff1y1*coeff2y1_2st+coeff1y2*coeff2y2_2st/2)/(tau*x[1])+(coeffy0*coeffy0_2st/2)/(tau*x[3])
  
  cov23_2 <- (coeff2y1*coeff2y1_2st+coeff2y2*coeff2y2_2st/2)/(tau*x[1])+(coeffy0*coeffy0_2st/2)/(tau*x[3])
  
  cov_matr212 <- matrix(c(sd1^2, cov12, cov13_2,
                          cov12, sd2^2, cov23_2,
                          cov13_2, cov23_2, sd22_2st^2), byrow = TRUE, nrow =3, ncol=3)
  
  # thresholds
  
  thr1u2_2st <- const1+qnorm(x[5])*post_sdtheta1
  thr1l2_2st <- -Inf
  
  thr2u2_2st <- const2+qnorm(x[4])*post_sdtheta2
  thr2l2_2st <- const2+qnorm(x[5])*post_sdtheta2
  
  thr3l2_2st <- const2_2st+qnorm(x[6])*post_sdtheta2_2st
  thr3u2_2st <- Inf
  
  # integrals
  
  int_2_2st <- pmvnorm(lower = c(thr1l2_2st,thr2l2_2st,thr3l2_2st),
                       upper = c(thr1u2_2st, thr2u2_2st,thr3u2_2st),
                       sigma = cov_matr212,
                       mean = c(mean1, mean2, mean22_2st),
                       algorithm = GenzBretz(abseps = 1*10^-20))[1]
  
  
  # continue with 2 arms
  
  A212_2st <- invomega+2*c(x[1],x[2])*invsigma
  
  invA212_2st <- solve(A212_2st)
  
  coeff1y1_12_2st <- invA212_2st[1,1]*(invsigma*2*x[1])[1,1] 
  
  coeff1y2_12_2st <- invA212_2st[1,2]*(invsigma*2*x[2])[2,2] 
  
  coeffy0_12_2st <- 2*x[3]*tau/(tau00+2*x[3]*tau) 
  
  const1_12_2st <- -(invA212_2st%*%invomega%*%(mu))[1,1] +const0_2st 
  
  const2_12_2st <- -(invA212_2st%*%invomega%*%(mu))[2,1] +const0_2st 
  
  coeff2y1_12_2st <- invA212_2st[2,1]*(invsigma*2*x[1])[1,1] 
  
  coeff2y2_12_2st <- invA212_2st[2,2]*(invsigma*2*x[2])[2,2] 
  
  post_varmu0_12_2st <- 1/(tau00+2*x[3]*tau) #posterior variance of mu0
  
  post_sdtheta1_12_2st <- sqrt(invA212_2st[1,1]+post_varmu0_12_2st) #posterior standard deviation of theta1
  
  post_sdtheta2_12_2st <- sqrt(invA212_2st[2,2]+post_varmu0_12_2st) #posterior standard deviation of theta2
  
  # means and standard deviation for multivariate normal 
  
  mean1_12_2st <- sum(c(coeff1y1_12_2st,coeff1y2_12_2st,-coeffy0_12_2st)*null)
  
  sd1_12_2st <- sqrt(sum(c(coeff1y1_12_2st,coeff1y2_12_2st)^2)/(tau*2*x[1])+coeffy0_12_2st^2/(tau*2*x[3]))
  
  mean2_12_2st <-  sum(c(coeff2y1_12_2st,coeff2y2_12_2st,-coeffy0_12_2st)*null)
  
  sd2_12_2st <- sqrt(sum(c(coeff2y1_12_2st,coeff2y2_12_2st)^2)/(tau*2*x[1])+coeffy0_12_2st^2/(tau*2*x[3]))
  
  # covariance matrix 
  
  cov13_3 <- (coeff1y1*coeff1y1_12_2st/2+coeff1y2*coeff1y2_12_2st/2)/(tau*x[1])+(coeffy0*coeffy0_12_2st/2)/(tau*x[3])
  
  cov23_3 <- (coeff2y1*coeff1y1_12_2st/2+coeff2y2*coeff1y2_12_2st/2)/(tau*x[1])+(coeffy0*coeffy0_12_2st/2)/(tau*x[3])
  
  cov13_4 <- (coeff1y1*coeff2y1_12_2st/2+coeff1y2*coeff2y2_12_2st/2)/(tau*x[1])+(coeffy0*coeffy0_12_2st/2)/(tau*x[3])
  
  cov23_4 <- (coeff2y1*coeff2y1_12_2st/2+coeff2y2*coeff2y2_12_2st/2)/(tau*x[1])+(coeffy0*coeffy0_12_2st/2)/(tau*x[3])
  
  cov34_4 <- (coeff2y1_12_2st*coeff1y1_12_2st/2+coeff2y2_12_2st*coeff1y2_12_2st/2)/(tau*x[1])+(coeffy0_12_2st*coeffy0_12_2st/2)/(tau*x[3])
  
  cov_matr1212 <- matrix(c(sd1^2, cov12, cov13_3,cov13_4,
                           cov12, sd2^2, cov23_3, cov23_4,
                           cov13_3, cov23_3,sd1_12_2st^2,cov34_4, 
                           cov13_4, cov23_4, cov34_4,sd2_12_2st^2), 
                         byrow = TRUE, nrow =4, ncol=4)
  # thresholds
  
  thr3l12_2st <- const1_12_2st+qnorm(x[6])*post_sdtheta1_12_2st
  thr3u12_2st <- Inf
  
  thr4l12_2st <- const2_12_2st+qnorm(x[6])*post_sdtheta2_12_2st
  thr4u12_2st <- Inf
  
  int12_1_2st <- pmvnorm(lower = c(thr1l,thr2l,thr3l12_2st,-Inf),
                         upper = c(thr1u,thr2u,thr3u12_2st,Inf),
                         sigma = cov_matr1212,
                         mean = c(mean1, mean2, mean1_12_2st, mean2_12_2st),
                         algorithm = GenzBretz(abseps = 1*10^-20))[1]
  
  int12_2_2st <-pmvnorm(lower = c(thr1l,thr2l,-Inf,thr4l12_2st),
                        upper = c(thr1u,thr2u,Inf,thr4u12_2st),
                        sigma = cov_matr1212,
                        mean = c(mean1, mean2, mean1_12_2st, mean2_12_2st),
                        algorithm = GenzBretz(abseps = 1*10^-20))[1]
  
  int12_12_2st <-  pmvnorm(lower = c(thr1l,thr2l,thr3l12_2st,thr4l12_2st),
                           upper = c(thr1u,thr2u,thr3u12_2st,thr4u12_2st),
                           sigma = cov_matr1212,
                           mean = c(mean1, mean2, mean1_12_2st, mean2_12_2st),
                           algorithm = GenzBretz(abseps = 1*10^-20))[1]
  
  rejatleastone=  int1null+int2null-int12null+int_1_2st+int_2_2st+int12_1_2st+int12_2_2st-int12_12_2st
  
  ################ power ###############
  
  mean1pow <- sum(c(coeff1y1,coeff1y2,-coeffy0)*power)
  
  mean2pow <- sum(c(coeff2y1,coeff2y2,-coeffy0)*power)
  
  mean1_12_2stpow <- sum(c(coeff1y1_12_2st,coeff1y2_12_2st,-coeffy0_12_2st)*power)
  
  mean2_12_2stpow <- sum(c(coeff2y1_12_2st,coeff2y2_12_2st,-coeffy0_12_2st)*power)
  
  mean11_2stpow <- sum(c(coeff1y1_2st,coeff1y2_2st,-coeffy0_2st)*power)
  
  mean22_2stpow <- sum(c(coeff2y1_2st,coeff2y2_2st,-coeffy0_2st)*power)
  
  int12pow <- pmvnorm(lower = c(thr1u, thr2u),
                      upper = c(Inf, Inf),
                      sigma = cov_matr12,
                      mean = c(mean1pow, mean2pow),
                      algorithm = GenzBretz(abseps = 1*10^-20))[1]
  
  int12pow_12_2st <-  pmvnorm(lower = c(thr1l,thr2l,thr3l12_2st,thr4l12_2st),
                              upper = c(thr1u,thr2u,thr3u12_2st,thr4u12_2st),
                              sigma = cov_matr1212,
                              mean = c(mean1pow, mean2pow, mean1_12_2stpow, mean2_12_2stpow),
                              algorithm = GenzBretz(abseps = 1*10^-20))[1]
  
  intpow_1_2st <- pmvnorm(lower = c(thr1l1_2st,const2+qnorm(x[4])*post_sdtheta2,thr3l1_2st),
                          upper = c(thr1u1_2st, Inf,thr3u1_2st),
                          sigma = cov_matr112,
                          mean = c(mean1pow, mean2pow, mean11_2stpow),
                          algorithm = GenzBretz(abseps = 1*10^-20))[1]
  
  intpow_2_2st <- pmvnorm(lower = c(const1+qnorm(x[4])*post_sdtheta1,thr2l2_2st,thr3l2_2st),
                          upper = c(Inf, thr2u2_2st,thr3u2_2st),
                          sigma = cov_matr212,
                          mean = c(mean1pow, mean2pow, mean22_2stpow),
                          algorithm = GenzBretz(abseps = 1*10^-20))[1]
  
  rejpower = int12pow+intpow_1_2st+intpow_2_2st+int12pow_12_2st
  
  
  ################ reject both under power2 ###############
  
  mean1pow2 <- sum(c(coeff1y1,coeff1y2,-coeffy0)*power2)
  
  mean2pow2 <- sum(c(coeff2y1,coeff2y2,-coeffy0)*power2)
  
  mean1_12_2stpow2 <- sum(c(coeff1y1_12_2st,coeff1y2_12_2st,-coeffy0_12_2st)*power2)
  
  mean2_12_2stpow2 <- sum(c(coeff2y1_12_2st,coeff2y2_12_2st,-coeffy0_12_2st)*power2)
  
  mean11_2stpow2 <- sum(c(coeff1y1_2st,coeff1y2_2st,-coeffy0_2st)*power2)
  
  mean22_2stpow2 <- sum(c(coeff2y1_2st,coeff2y2_2st,-coeffy0_2st)*power2)
  
  int12pow2 <- pmvnorm(lower = c(thr1u, thr2u),
                       upper = c(Inf, Inf),
                       sigma = cov_matr12,
                       mean = c(mean1pow2, mean2pow2),
                       algorithm = GenzBretz(abseps = 1*10^-20))[1]
  
  int12pow_12_2st2 <-  pmvnorm(lower = c(thr1l,thr2l,thr3l12_2st,thr4l12_2st),
                               upper = c(thr1u,thr2u,thr3u12_2st,thr4u12_2st),
                               sigma = cov_matr1212,
                               mean = c(mean1pow2, mean2pow2, mean1_12_2stpow2, mean2_12_2stpow2),
                               algorithm = GenzBretz(abseps = 1*10^-20))[1]
  
  intpow_1_2st2 <- pmvnorm(lower = c(thr1l1_2st,const2+qnorm(x[4])*post_sdtheta2,thr3l1_2st),
                           upper = c(thr1u1_2st, Inf,thr3u1_2st),
                           sigma = cov_matr112,
                           mean = c(mean1pow2, mean2pow2, mean11_2stpow2),
                           algorithm = GenzBretz(abseps = 1*10^-20))[1]
  
  intpow_2_2st2 <- pmvnorm(lower = c(const1+qnorm(x[4])*post_sdtheta1,thr2l2_2st,thr3l2_2st),
                           upper = c(Inf, thr2u2_2st,thr3u2_2st),
                           sigma = cov_matr212,
                           mean = c(mean1pow2, mean2pow2, mean22_2stpow2),
                           algorithm = GenzBretz(abseps = 1*10^-20))[1]
  
  rejpowernull1 = int12pow2+intpow_1_2st2+intpow_2_2st2+int12pow_12_2st2
  
  
  ################ reject both under power3 ###############
  
  mean1pow3 <- sum(c(coeff1y1,coeff1y2,-coeffy0)*power3)
  
  mean2pow3 <- sum(c(coeff2y1,coeff2y2,-coeffy0)*power3)
  
  mean1_12_2stpow3 <- sum(c(coeff1y1_12_2st,coeff1y2_12_2st,-coeffy0_12_2st)*power3)
  
  mean2_12_2stpow3 <- sum(c(coeff2y1_12_2st,coeff2y2_12_2st,-coeffy0_12_2st)*power3)
  
  mean11_2stpow3 <- sum(c(coeff1y1_2st,coeff1y2_2st,-coeffy0_2st)*power3)
  
  mean22_2stpow3 <- sum(c(coeff2y1_2st,coeff2y2_2st,-coeffy0_2st)*power3)
  
  int12pow3 <- pmvnorm(lower = c(thr1u, thr2u),
                       upper = c(Inf, Inf),
                       sigma = cov_matr12,
                       mean = c(mean1pow3, mean2pow3),
                       algorithm = GenzBretz(abseps = 1*10^-20))[1]
  
  int12pow_12_2st3 <-  pmvnorm(lower = c(thr1l,thr2l,thr3l12_2st,thr4l12_2st),
                               upper = c(thr1u,thr2u,thr3u12_2st,thr4u12_2st),
                               sigma = cov_matr1212,
                               mean = c(mean1pow3, mean2pow3, mean1_12_2stpow3, mean2_12_2stpow3),
                               algorithm = GenzBretz(abseps = 1*10^-20))[1]
  
  intpow_1_2st3 <- pmvnorm(lower = c(thr1l1_2st,const2+qnorm(x[4])*post_sdtheta2,thr3l1_2st),
                           upper = c(thr1u1_2st, Inf,thr3u1_2st),
                           sigma = cov_matr112,
                           mean = c(mean1pow3, mean2pow3, mean11_2stpow3),
                           algorithm = GenzBretz(abseps = 1*10^-20))[1]
  
  intpow_2_2st3 <- pmvnorm(lower = c(const1+qnorm(x[4])*post_sdtheta1,thr2l2_2st,thr3l2_2st),
                           upper = c(Inf, thr2u2_2st,thr3u2_2st),
                           sigma = cov_matr212,
                           mean = c(mean1pow3, mean2pow3, mean22_2stpow3),
                           algorithm = GenzBretz(abseps = 1*10^-20))[1]
  
  rejpowernull2 = int12pow3+intpow_1_2st3+intpow_2_2st3+int12pow_12_2st3
  
  res <- list(x[1], x[4],x[5],x[6], rejatleastone, rejpower, 
              rejpowernull1, rejpowernull2)
  
  names(res) <- c("n11","eta1","l1","eta2", "rejatleastone", "rejpower", 
                  "rejectpowernull1","rejectpowernull2")
  return(res)
  
  
}

##### Function to compute integrals for each value of the grid x and grid_n11eta ####

integ_function <- function(x, grid_n11eta){
  
  taud <- x[1]
  
  meand <- x[2]
  
  null <- c(x[3],x[3],x[3])
  
  power <- c(trteffect,trteffect,0)+null
  
  power2 <- c(trteffect,0,0)+null
  
  power3 <- c(0,trteffect,0)+null
  
  st = Sys.time()
  
  out <- apply(grid_n11eta, 1, sim,
               tau=tau, tau00=tau00, taud=taud, mu00=mu00, mu01=mu01,
               tau01=tau01, alpha=alpha, beta=beta,
               meand=meand, null=null,power=power, power2=power2,power3=power3)
  out
  
  outmatrix <- data.frame(matrix(unlist(out), ncol = 8, byrow=TRUE))
  
  
  colnames(outmatrix) <- c("n11","eta1","eps1","eta2", 
                           "rejatleastone", "rejpower", 
                           "rejectpowernull1","rejectpowernull2")
  
  outmatrix <- outmatrix  %>%  mutate(
    SS = 4*n11+2*n11,
    diff_FWER_alpha= abs(rejatleastone-alpha)) 
  
  end = Sys.time()
  print(end-st)
  
  truemu00 <- x[3]
  
  res <- cbind(outmatrix,meand, taud, truemu00,
               mu01, tau01, mu00, tau00, tau)
  
  return(res)
  
}

