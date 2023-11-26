SEIRD8V <- function( S0, E0, I0, II0,
                     RE0, RI0,  RR0, D0, V0, alpha1, beta1,  gamma1, gamma2,  eta1, 
                     rho1,  phi1,  zeta1, zeta2, kappa1, n1, mchpt1){
  # Here alpha1 must be one dimension higher than chpt1
  ##### add ho1-3 in vector form too! Others are not
  
  X1a <- mchpt1$X1a   #alpha
  X1b <- mchpt1$X1b   #beta
  X1g <- mchpt1$X1g   #gamma
  X1p1 <- mchpt1$X1p1   #phi1
  X1m <- mchpt1$X1m  # mu
  # X1mI <- mchpt1$X1mI   #muI
  
  
  #beta1v <- X1b*betaI1
  
  alpha11v <- X1a%*% alpha1
  #alpha11v <- X1a[,1]*alpha1[1] + X1a[,2]*alpha1[2]
  rho1v <- X1m[,1]*rho1[1] + X1m[,2]*rho1[2]
  #rho1Iv <- X1mI*rho1I
  gamma1v <- X1g%*%gamma1 
  phi1v <- X1p1%*%phi1
  #gamma1v <- X1g[,1]*gamma1[1] + X1g[,2]*gamma1[2]
  
  
  Out2a <- SEIRDV8cpp(S0, E0, I0, II0,
                      RE0, RI0,  RR0, D0, V0, alpha11v, beta1, gamma1v, 
                      gamma2, eta1,
                      rho1v, phi1v, zeta1, zeta2, kappa1, n1 )
  
  Out2 <- data.frame( S = Out2a$S,
                      E = Out2a$E,
                      I = Out2a$I,
                      II = Out2a$II,
                      
                      RE = Out2a$RE,
                      RI = Out2a$RI,
                     
                      RR = Out2a$RR,
                      
                      D =  Out2a$D,
                      V = Out2a$V,
                      
                      Rep0 = Out2a$R0)
  return( Out2 )
}


LikeSEIRD8V2 <- function( data3, S0, E0, I0, II0, 
                          RE0, RI0,  RR0, D0, V0, alpha1, beta1,  gamma1, 
                          gamma2, eta1, 
                          rho1,  phi1, zeta1, zeta2, kappa1, n1, mchpt1){
  # Call the solver
  Out2 <- SEIRD8V(  S0, E0, I0, II0, 
                    RE0, RI0,  RR0, D0, V0, alpha1, beta1,  gamma1, gamma2,
                    eta1, 
                    rho1,  phi1, zeta1, zeta2, kappa1, n1, mchpt1=mchpt1 )
  
  # Evaluate the likelihood at each point.
  ##### write the total recovered and total reinfection here
  
  Out2$RE[Out2$RE==0] <- 1.0e-15
  Out2$II[Out2$II==0] <- 1.0e-15
  Out2$I[Out2$I==0] <- 1.0e-15
  Out2$RI[Out2$RI==0] <- 1.0e-15
  
  Out2$RR[Out2$RR==0] <- 1.0e-15
  # 
  # 
  Out2$V[Out2$V==0] <- 1.0e-15
  # 
  
  Lik2I <- sum( dpois( data3$AdjInfect, Out2$I, log = TRUE ) )  #infected
  Lik2II <- sum( dpois( data3$CReinfection, Out2$II, log = TRUE ) )  #reinfection
  Lik2RI <- sum( dpois( data3$CRecovered, Out2$RI, log = TRUE ) )  #recovered
  Lik2RR <- sum( dpois( data3$CReRecovered, Out2$RR, log = TRUE ) )  #Re_recovered
  Lik2D <- sum( dpois( data3$CDeaths, Out2$D, log = TRUE ) )
  Lik2V <- sum( dpois( data3$CVaccine, Out2$V, log = TRUE ) )
  
  res2 <- Lik2I +  Lik2II + Lik2RI + Lik2RR + Lik2D + Lik2V
  
  return( res2 )
}



PostSEIRD8V1 <- function( data3, S0, E0, I0, II0,
                          RE0, RI0,  RR0, D0, V0, alpha1, beta1,  gamma1, 
                          gamma2, eta1, 
                          rho1,  phi1,  zeta1, zeta2, kappa1, n1, mchpt1,
                          prior1am=prior1am, prior1b=prior1b, prior1g=prior1g,
                          prior1g2=prior1g2, prior1e=prior1e, prior1m=prior1m,
                          prior1p1=prior1p1, prior1z1=prior1z1, 
                          prior1z2=prior1z2, prior1k=prior1k ){
  #
  if( min(gamma1 )  > 0 & 
      min( alpha1) > 0  ){
    Post2 <- sum(dexp( eta1, prior1e, log = TRUE )) 
    Post2 <- Post2 + sum( dexp( rho1, prior1m, log = TRUE ) )
    Post2 <- Post2 + sum( dexp( beta1, prior1b, log = TRUE ) )
    Post2 <- Post2 + sum( dexp( gamma1, prior1g, log = TRUE ) )
    Post2 <- Post2 + sum( dexp( gamma2, prior1g2, log = TRUE ) )
    Post2 <- Post2 + sum( dexp( phi1, prior1p1, log = TRUE) )
    Post2 <- Post2 + sum( dexp( zeta1, prior1z1, log = TRUE ) )
    Post2 <- Post2 + sum( dexp( zeta2, prior1z2, log = TRUE ) )
    Post2 <- Post2 + sum( dexp( alpha1, prior1am, log = TRUE) )
    Post2 <- Post2 + sum( dbeta( kappa1, prior1k, 1, log = TRUE ) )
    

    Lik2 <- LikeSEIRD8V2( data3, S0, E0, I0, II0,  
                          RE0, RI0,  RR0, D0, V0, alpha1, beta1,  gamma1,
                          gamma2, eta1, rho1,  phi1,  zeta1, zeta2, kappa1, n1, mchpt1=mchpt1)
    res2 <- Post2 + Lik2
  }else{
    res2 <- -Inf
  }
  return( res2 )
}



# Build the matrix for change points.
MatrixBuild1 <- function( chpt1, n1 ){
  l2 <- length( chpt1 )
  res2 <- matrix( 0, nrow = n1, ncol = l2 )
  for( i in 1:(l2-1)){
    res2[(chpt1[i]:(chpt1[i+1]-1)),i] <- 1 
  }
  res2[ chpt1[l2]:n1, l2 ] <- 1
  return(res2)
}



MCMCSEIRD8VPred <- function( data3, S0, E0, I0, II0,
                             RE0, RI0,  RR0, D0, V0, alpha1, beta1,  gamma1,  
                             gamma2, eta1, rho1,  phi1, zeta1, zeta2, kappa1, n1, chpt1, rchpt1, gchpt1,  
                             chpt2, nMCMC1, alpha1Step, beta1Step, 
                             gamma1Step, gamma2Step, eta1Step, rho1Step, rho1IStep,
                             phi1Step, zeta1Step=zeta1Step, zeta2Step, kappa1Step, 
                             prior1am=prior1am, prior1b=prior1b, prior1g=prior1g,
                             prior1g2=prior1g2, prior1e=prior1e, prior1m=prior1m,
                             prior1p1=prior1p1, prior1z1=prior1z1, 
                             prior1z2=prior1z2, prior1k=prior1k){
  
  X1a <- MatrixBuild1( chpt1, n1 )
  X1g <- MatrixBuild1( gchpt1, n1 )
  X1m <- MatrixBuild1( rchpt1, n1 )
  X1p1 <- MatrixBuild1( chpt2, n1 )
  
  X1g2 <- rep(0,n1)
  X1z1 <- rep(0, n1)
  X1z2 <- rep(0, n1)
  X1e <- rep(0, n1)
  X1b <- rep(0, n1)
  #X1b[ImpI1] <- 1
  # X1mI <- rep(0, n1)
  # X1mI[ImpRI1 ] <- 1
  # 
  
  
  #mchpt1 <- list( X1a = X1a, X1b = X1b, X1g = X1g, X1m = X1m,  X1mI= X1mI)
  mchpt1 <- list( X1a = X1a, X1g = X1g, X1m = X1m, X1p1=X1p1)
  
  Post2 <- PostSEIRD8V1( data3, S0, E0, I0, II0, 
                         RE0, RI0,  RR0, D0, V0, alpha1, beta1,  gamma1,
                         gamma2, eta1, 
                         rho1, phi1,  zeta1, zeta2, kappa1, n1, mchpt1,
                         prior1am=prior1am, prior1b=prior1b, prior1g=prior1g,
                         prior1g2=prior1g2, prior1e=prior1e, prior1m=prior1m,
                         prior1p1=prior1p1, prior1z1=prior1z1, 
                         prior1z2=prior1z2, prior1k=prior1k)
  
  alpha1Out <- matrix( 0, nrow = nsamp1, ncol = length( alpha1 ) )
  beta1Out <- rep( 0, nsamp1)
  gamma1Out <- matrix( 0, nrow = nsamp1, ncol = length( gamma1 ) )
  eta1Out <- rep( 0, nsamp1 )
  zeta1Out <- rep( 0, nsamp1 )
  gamma2Out <- rep(0, nsamp1)
  zeta2Out <- rep( 0, nsamp1 )
  phi1Out <- matrix(0, nrow = nsamp1, ncol= length(phi1) )
  rho1Out <- matrix( 0, nrow = nsamp1, ncol = length( rho1 ) )
  kappa1Out <- rep( 0, nsamp1 )
  
  
  Post2Out <- rep(0, nsamp1 )
  IPred1 <- matrix( 0, nrow = nsamp1, ncol = n1 ) 
  IIPred1 <- matrix( 0, nrow = nsamp1, ncol = n1 )
  
  REPred1 <- matrix( 0, nrow = nsamp1, ncol = n1 )
  RIPred1 <- matrix( 0, nrow = nsamp1, ncol = n1 )
  RRPred1 <- matrix( 0, nrow = nsamp1, ncol = n1 )
  
  DPred1 <- matrix( 0, nrow = nsamp1, ncol = n1 )
  VPred1 <- matrix( 0, nrow = nsamp1, ncol = n1 )
  
  
  # MCMC (Metropolis Hasting sampler) for estimating the parameters
  
  for( i in 1:nsamp1 ){
    for( j in 1:length( alpha1 ) ){
      alpha1t <- alpha1
      #  Draw a candidate for random walk chain
      alpha1t[j] <- alpha1[j] + alpha1Step[j]*rnorm( 1, 0, 1 )
      if( alpha1t[j] > 0 ){
        Post2t <- PostSEIRD8V1( data3, S0, E0, I0, II0, 
                                RE0, RI0,  RR0, D0, V0, alpha1t, beta1,  gamma1, 
                                gamma2, eta1, rho1,  phi1, zeta1, zeta2, kappa1, n1, mchpt1,
                                prior1am=prior1am, prior1b=prior1b, prior1g=prior1g,
                                prior1g2=prior1g2, prior1e=prior1e, prior1m=prior1m,
                                prior1p1=prior1p1, prior1z1=prior1z1, 
                                prior1z2=prior1z2, prior1k=prior1k)
        diff1 <- Post2t - Post2
        U1 <- log( runif( 1, 0, 1) )
        # to accept candidate draw with probability diff1
        if( diff1 > U1 ){
          Post2 <- Post2t
          alpha1 <- alpha1t 
        }
      }
    }
    
    beta1t <- beta1 + beta1Step*rnorm( 1, 0, 1 )
    if( beta1t > 0 ){
      Post2t <- PostSEIRD8V1(  data3, S0, E0, I0, II0, 
                               RE0, RI0,  RR0, D0, V0, alpha1, beta1t,  gamma1, 
                               gamma2, eta1, 
                               rho1,  phi1,  zeta1, zeta2, kappa1, n1, mchpt1,
                               prior1am=prior1am, prior1b=prior1b, prior1g=prior1g,
                               prior1g2=prior1g2, prior1e=prior1e, prior1m=prior1m,
                               prior1p1=prior1p1, prior1z1=prior1z1, 
                               prior1z2=prior1z2, prior1k=prior1k)
      diff1 <- Post2t - Post2
      U1 <- log( runif( 1, 0, 1) )
      if( diff1 > U1 ){
        Post2 <- Post2t
        beta1 <- beta1t 
      }
    }
    
    # betaI1t <- betaI1 + betaI1Step*rnorm( 1, 0, 1 )
    # if( betaI1t > 0  ){
    #   Post2t <- PostSEIRD8V1(  data3, S0, E0, I0, II0, 
    #                            RE0, RI0, RR0, D0, V0, alpha1, beta1, betaI1t, gamma1, 
    #                             eta1, 
    #                            rho1,  phi1,  n1, mchpt1,
    #                            prior1am, prior1b, prior1g,  prior1e, prior1m,
    #                            prior1p1 = prior1p1) 
    #   diff1 <- Post2t - Post2
    #   U1 <- log( runif( 1, 0, 1) )
    #   if( diff1 > U1 ){
    #     Post2 <- Post2t
    #     betaI1 <- betaI1t 
    #   }
    # }
    
    for( j in 1:length( gamma1 ) ){
      gamma1t <- gamma1
      gamma1t[j] <- gamma1[j] + gamma1Step[j]*rnorm( 1, 0, 1 ) 
      if( gamma1t[j] > 0) {
        Post2t <- PostSEIRD8V1(  data3, S0, E0, I0, II0, 
                                 RE0, RI0,  RR0, D0, V0, alpha1, beta1,  gamma1t, 
                                 gamma2, eta1, 
                                 rho1,  phi1,  zeta1, zeta2, kappa1, n1, mchpt1,
                                 prior1am=prior1am, prior1b=prior1b, prior1g=prior1g,
                                 prior1g2=prior1g2, prior1e=prior1e, prior1m=prior1m,
                                 prior1p1=prior1p1, prior1z1=prior1z1, 
                                 prior1z2=prior1z2, prior1k=prior1k)
        diff1 <- Post2t - Post2
        U1 <- log( runif( 1, 0, 1) )
        if( diff1 > U1 ){
          Post2 <- Post2t
          gamma1 <- gamma1t 
        }
      }
    }
    
    eta1t <- eta1 + eta1Step*rnorm( 1, 0, 1 )
    if( eta1t > 0 ){
      Post2t <- PostSEIRD8V1(  data3, S0, E0, I0, II0, 
                               RE0, RI0,  RR0, D0, V0, alpha1, beta1,  gamma1, 
                               gamma2, eta1t, 
                               rho1,  phi1, zeta1, zeta2, kappa1, n1, mchpt1,
                               prior1am=prior1am, prior1b=prior1b, prior1g=prior1g,
                               prior1g2=prior1g2, prior1e=prior1e, prior1m=prior1m,
                               prior1p1=prior1p1, prior1z1=prior1z1, 
                               prior1z2=prior1z2, prior1k=prior1k)
      diff1 <- Post2t - Post2
      U1 <- log( runif( 1, 0, 1) )
      if( diff1 > U1 ){
        Post2 <- Post2t
        eta1 <- eta1t 
      }
    }
    
    
    zeta1t <- zeta1 + zeta1Step*rnorm( 1, 0, 1 )
    if( zeta1t > 0 ){
      Post2t <- PostSEIRD8V1(  data3, S0, E0, I0, II0, 
                               RE0, RI0,  RR0, D0, V0, alpha1, beta1,  gamma1, 
                               gamma2, eta1, 
                               rho1,  phi1, zeta1t, zeta2, kappa1, n1, mchpt1,
                               prior1am=prior1am, prior1b=prior1b, prior1g=prior1g,
                               prior1g2=prior1g2, prior1e=prior1e, prior1m=prior1m,
                               prior1p1=prior1p1, prior1z1=prior1z1, 
                               prior1z2=prior1z2, prior1k=prior1k)
      diff1 <- Post2t - Post2
      U1 <- log( runif( 1, 0, 1) )
      if( diff1 > U1 ){
        Post2 <- Post2t
        zeta1 <- zeta1t 
      }
    }
    
    zeta2t <- zeta2 + zeta2Step*rnorm( 1, 0, 1 )
    if( zeta2t > 0 ){
      Post2t <- PostSEIRD8V1(  data3, S0, E0, I0, II0, 
                               RE0, RI0,  RR0, D0, V0, alpha1, beta1,  gamma1, 
                               gamma2, eta1, 
                               rho1,  phi1, zeta1, zeta2t, kappa1, n1, mchpt1,
                               prior1am=prior1am, prior1b=prior1b, prior1g=prior1g,
                               prior1g2=prior1g2, prior1e=prior1e, prior1m=prior1m,
                               prior1p1=prior1p1, prior1z1=prior1z1, 
                               prior1z2=prior1z2, prior1k=prior1k)
      diff1 <- Post2t - Post2
      U1 <- log( runif( 1, 0, 1) )
      if( diff1 > U1 ){
        Post2 <- Post2t
        zeta2 <- zeta2t 
      }
    }
    
    for(j in 2:length(rho1)){
      rho1t <- rho1
      rho1t[j] <- rho1[j] + rho1Step[j]*rnorm( 1, 0, 1 )
      if( rho1t[j] > 0 ){
        Post2t <- PostSEIRD8V1(  data3, S0, E0, I0, II0, 
                                 RE0, RI0,  RR0, D0, V0, alpha1, beta1,  gamma1, 
                                 gamma2, eta1, 
                                 rho1t,  phi1, zeta1, zeta2, kappa1, n1, mchpt1,
                                 prior1am=prior1am, prior1b=prior1b, prior1g=prior1g,
                                 prior1g2=prior1g2, prior1e=prior1e, prior1m=prior1m,
                                 prior1p1=prior1p1, prior1z1=prior1z1, 
                                 prior1z2=prior1z2, prior1k=prior1k )
        diff1 <- Post2t - Post2
        U1 <- log( runif( 1, 0, 1) )
        if( diff1 > U1 ){
          Post2 <- Post2t
          rho1 <- rho1t 
        }
      }
    }
    
    gamma2t <- gamma2
    gamma2t <- gamma2 + gamma2Step*rnorm( 1, 0, 1 )
    if( gamma2t > 0 ){
      Post2t <- PostSEIRD8V1(  data3, S0, E0, I0, II0, 
                               RE0, RI0,  RR0, D0, V0, alpha1, beta1,  gamma1, 
                               gamma2t, eta1, 
                               rho1t,  phi1, zeta1, zeta2, kappa1, n1, mchpt1,
                               prior1am=prior1am, prior1b=prior1b, prior1g=prior1g,
                               prior1g2=prior1g2, prior1e=prior1e, prior1m=prior1m,
                               prior1p1=prior1p1, prior1z1=prior1z1, 
                               prior1z2=prior1z2, prior1k=prior1k )
      diff1 <- Post2t - Post2
      U1 <- log( runif( 1, 0, 1) )
      if( diff1 > U1 ){
        Post2 <- Post2t
        gamma2 <- gamma2t
      }
    }
 
    
    kappa1t <- kappa1
    kappa1t <- kappa1 + kappa1Step*rnorm( 1, 0, 1 )
    if( kappa1t > 0 ){
      Post2t <- PostSEIRD8V1(  data3, S0, E0, I0, II0, 
                               RE0, RI0,  RR0, D0, V0, alpha1, beta1,  gamma1, 
                               gamma2, eta1, 
                               rho1t,  phi1, zeta1, zeta2, kappa1t, n1, mchpt1,
                               prior1am=prior1am, prior1b=prior1b, prior1g=prior1g,
                               prior1g2=prior1g2, prior1e=prior1e, prior1m=prior1m,
                               prior1p1=prior1p1, prior1z1=prior1z1, 
                               prior1z2=prior1z2, prior1k=prior1k )
      diff1 <- Post2t - Post2
      U1 <- log( runif( 1, 0, 1) )
      if( diff1 > U1 ){
        Post2 <- Post2t
        kappa1 <- kappa1t
      }
    }
    
    
    
    for( j in 1:length( phi1 ) ){
      phi1t <- phi1 
      phi1t[j] <- phi1[j] + phi1Step[j]*rnorm( 1, 0, 1 ) 
      if( phi1t[j] > 0) {
        Post2t <- PostSEIRD8V1(  data3, S0, E0, I0, II0, 
                                 RE0, RI0,  RR0, D0, V0, alpha1, beta1,  gamma1, 
                                 gamma2, eta1, 
                                 rho1,  phi1t,  zeta1, zeta2, kappa1, n1, mchpt1,
                                 prior1am=prior1am, prior1b=prior1b, prior1g=prior1g,
                                 prior1g2=prior1g2, prior1e=prior1e, prior1m=prior1m,
                                 prior1p1=prior1p1, prior1z1=prior1z1, 
                                 prior1z2=prior1z2, prior1k=prior1k)
        diff1 <- Post2t - Post2
        U1 <- log( runif( 1, 0, 1) )
        if( diff1 > U1 ){   #acceptance criteria
          Post2 <- Post2t
          phi1 <- phi1t 
        }
      }
    }
    
    
    
    Fit2Q <- SEIRD8V( S0, E0, I0, II0, 
                      RE0, RI0,  RR0, D0, V0, alpha1, beta1, gamma1,  
                      gamma2, eta1, rho1,  phi1, zeta1, zeta2, kappa1, n1, mchpt1 )
    
    
    # Fit2Q$II[Fit2Q$II==0] <- 2.672700e-06
    # Fit2Q$I[Fit2Q$I==0] <- 2.672700e-05
    # Fit2Q$RI[Fit2Q$RI==0] <- 2.672700e-05
    # Fit2Q$RE[Fit2Q$RE==0] <- 2.672700e-05
    # #Fit2Q$D[Fit2Q$D==0] <- 2.672700e-05
    # Fit2Q$V[Fit2Q$V==0] <- 2.672700e-04
    
    IPred1[i,] <- ifelse( Fit2Q$I < 1000, rpois( n1, Fit2Q$I), round(rnorm(n1, Fit2Q$I, sqrt(Fit2Q$I)),0))
    IIPred1[i,] <- ifelse( Fit2Q$II < 1000, rpois( n1, Fit2Q$II), round(rnorm(n1, Fit2Q$II, sqrt(Fit2Q$II)),0))
    
    REPred1[i,] <- ifelse( Fit2Q$RE < 1000, rpois( n1, Fit2Q$RE), round(rnorm(n1, Fit2Q$RE, sqrt(Fit2Q$RE)),0))
    RIPred1[i,] <- ifelse( Fit2Q$RI < 1000, rpois( n1, Fit2Q$RI), round(rnorm(n1, Fit2Q$RI, sqrt(Fit2Q$RI)),0))
    RRPred1[i,] <- ifelse( Fit2Q$RR < 1000, rpois( n1, Fit2Q$RR), round(rnorm(n1, Fit2Q$RR, sqrt(Fit2Q$RR)),0))
    
    
    DPred1[i,] <- ifelse( Fit2Q$D < 1000, rpois( n1, Fit2Q$D), round(rnorm(n1, Fit2Q$D, sqrt(Fit2Q$D)),0))
    VPred1[i,] <- ifelse( Fit2Q$V < 1000, rpois( n1, Fit2Q$V), round(rnorm(n1, Fit2Q$V, sqrt(Fit2Q$V)),0))
    #Out1IN[i,] <- diff( Out1I[i,] )
    #Out1Max1[i] <- which( Fit2Q$I == max(Fit1Q$I) )
    
    
    alpha1Out[i,] <- alpha1
    beta1Out[i] <- beta1
    kappa1Out[i] <- kappa1
    gamma1Out[i,] <- gamma1
    eta1Out[i] <- eta1
    gamma2Out[i] <- gamma2
    zeta1Out[i] <- zeta1
    zeta2Out[i] <- zeta2
    rho1Out[i,] <- rho1
    phi1Out[i, ] <- phi1
    
    
    
    Post2Out[i] <- Post2
    
  }
  res2 <- list( alpha1 = alpha1Out,
                beta1 = beta1Out,
                #betaI1 = betaI1Out,
                gamma1 = gamma1Out,
                gamma2 = gamma2Out,
                eta1 = eta1Out,
                rho1 = rho1Out,
                kappa1 = kappa1Out,
                phi1 = phi1Out,
                zeta1 = zeta1Out,
                zeta2 = zeta2Out,
                Post2 = Post2Out,
                
                IPred1 = IPred1,
                IIPred1 = IIPred1,
                
                REPred1 = REPred1,
                RIPred1 = RIPred1,
                RRPred1 = RRPred1,
                
                DPred1 = DPred1,
                VPred1 = VPred1,
                
                Rep0 =    Rep0)
  
  return( res2 )
}




