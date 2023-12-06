################################################################################
########### This analysis no S2 compartment with some interventions--- Model 2 in the paper

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# This code produces the results in model 2 (SEIRDVI_I2R_R model) 
#Note: This code is similar in structure to model 1 code, however, the mode and results are entirely different.
#Also, the start values and initial values of the parameters are different, so be careful not to mess up the analysis here.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#Load the required library
library( Rcpp )

###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read in the functions we will use %%%%%%%%%%%%%%%%%%%%%%%%%%%%
sourceCpp( "SVEIRD.cpp" )

source("SVEIRD.R")

#Load the data
Cdata3 <- read.csv("DataForAnalysis.csv", header = TRUE, sep = ",", na.strings = "NA")
dim(Cdata3)

data3 <- Cdata3[1:550,]  # this data is until June 30 2022

n1Q <- nrow(data3)
Q1 <- data3


# # Adjust Infected to remove recovered and remove deaths
Q1$AdjInfect <-  cumsum(Q1$Number.of.new.infections) - Q1$CRecovered - Q1$CDeaths 


colnames(Q1)
data3 <- Q1
#View(data3)
n1 <- n1Q
# Starting Values
# Data for Qatar
S0 <- 2695122    #total population as of 2022 
E0 <- 5 #exposed

I0 <- 1 #Infected
II0 <- 0 #Reinfected
#Recovered
RE0 <- 0 
RI0 <- 0
#Re-recovered
RR0 <- 0

D0 <- 0  #death

V0 <- 0 ##vaccination



#%%%%%%%%%%%%% priors initial values %%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prior1am = 1 #prior for transmission rate (alpha)
prior1b = 1 #prior for Infectious rate (beta)
prior1g = 1 # prior for first recovery rate (gamma1)
prior1g2 = 1 #prior for second recovery rate (gamma2)
prior1e = 1 # prior for death rate (eta)
prior1m = 1  # prior for vaccination rate (mu=rho in the paper)
prior1p1 = 1  #prior for reinfection rate (phi)
prior1z1 = 1 #prior natural immunity waning rate (zeta1) 
prior1z2 = 1 #prior for vaccination protection waning rate (zeta2)
prior1k = 1 #prior for vaccination efficacy (kappa, 0 < kappa <= 1)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Parameter initial values and change points
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta1 <-  0.06987875
eta1 <- 0.0002108282
gamma2 <- 1/7
zeta1 <- 0.00116255   #90-3
zeta2 <- 0.00116255
kappa1 <- 0.94

########################################################### Steps
beta1Step <- 0.0001
gamma2Step <- 0.0001
eta1Step <- 0.0000001
zeta1Step <- 0.000001
zeta2Step <- 0.000001
kappa1Step <- 0.0001

rchpt1 <- c( 1, 390)  # Time at which vaccine1 was deployed
rho1 <- c(0,  0.004597066)  #1/100 persons
rho1Step <- c(0,0.0001)
chpt2 <- c( 117, 350, 450, 500)
phi1 <- c(0.005377395, 0.005377395, 0.005377395, 0.005377395)
phi1Step <- c( 1.00e-06, 1.00e-6, 1.00e-06, 1.00e-6)

chpt1 <- c(  20, 40, 80, 95, 160)
chpt1 <- c(chpt1, 170, 190, 250, 325, 360, 448, 500)

alpha1 <- c(2.861767e-07, 8.415194e-08, 9.087136e-08, 7.590418e-08, 1.099209e-07)
alpha1 <- c(alpha1, 8.856122e-08, 8.577855e-08, 9.725080e-08, 1.122183e-07, 1.321569e-07, 2.576241e-07, 5.008670e-07)
alpha1Step <- c( 1.00e-11, 1.00e-11, 1.00e-11, 1.00e-11, 1.00e-11)
alpha1Step <- c( alpha1Step, 1.00e-11, 1.00e-11, 1.00e-11, 1.00e-11,  1.00e-11, 1.00e-11,  1.00e-11)

gchpt1 <- c(6, 60, 350, 428)
gamma1 <- c(0.0909695, 0.1356373, 0.1475781, 0.1133254)
gamma1Step <- rep(1.00e-4, length(gamma1) )

chpt2 <- c( 117, 350, 450, 500)
phi1 <- c(0.005377395, 0.011865828, 0.022132963, 0.007401364)
phi1Step <- c( 1.00e-06, 1.00e-6, 1.00e-06, 1.00e-6)


####### Read in the start step values from a file...

filename8_SVEIRD <- "StartStepValuesMCSVEIRD.RData"
#load(filename8_SVEIRD) # Uncomment this if you want to load my startvalues.



tic <- Sys.time()
nsamp1 = 50000
nMCMC1 <- nsamp1
Q1MCMC <- MCMCSEIRD8VPred(data3, S0, E0, I0, II0,
                          RE0, RI0, RR0, D0, V0, alpha1, beta1, gamma1,
                          gamma2, eta1, rho1, phi1, zeta1, zeta2, kappa1, n1, chpt1, rchpt1, gchpt1,  
                          chpt2, nMCMC1, alpha1Step, beta1Step,
                          gamma1Step=gamma1Step, gamma2Step=gamma2Step, eta1Step=eta1Step, 
                          rho1Step=rho1Step, phi1Step= phi1Step, zeta1Step=zeta1Step,
                          zeta2Step=zeta2Step, kappa1Step=kappa1Step,
                          prior1am=prior1am, prior1b=prior1b, prior1g=prior1g, 
                          prior1g2=prior1g2,
                          prior1e=prior1e, prior1m=prior1m,
                          prior1p1= prior1p1, prior1z1=prior1z1, prior1z2=prior1z2, prior1k=prior1k)
Sys.time() - tic





#######################################################################
#
#  Trace plots to ensure convergence 
#
#######################################################################
for( i in 1:length(alpha1) ){
  plot( Q1MCMC$alpha1[,i], type = "l" )
}
plot( Q1MCMC$beta1, type = "l" )
for( i in 1:length( gamma1 ) ){
  plot( Q1MCMC$gamma1[,i], type = "l" )
}
plot( Q1MCMC$gamma2, type = "l" )
plot( Q1MCMC$kappa1, type = "l" )
plot( Q1MCMC$eta1, type = "l" )
plot( Q1MCMC$rho1[,1], type = "l" )
plot( Q1MCMC$rho1[,2], type = "l" )
#plot( Q1MCMC$rho1I, type = "l" )

for( i in 1:length( phi1 ) ){
  plot( Q1MCMC$phi1[,i], type = "l" )
}
plot( Q1MCMC$zeta1, type = "l" )
plot( Q1MCMC$zeta2, type = "l" )



# Create files so you can read these in 
alpha1Step <- apply( Q1MCMC$alpha1, 2, sd )/4
gamma1Step <- apply( Q1MCMC$gamma1, 2, sd )/4
gamma2Step <- sd(Q1MCMC$gamma2)/4
beta1Step <- sd( Q1MCMC$beta1 )/4
eta1Step <- sd( Q1MCMC$eta1 )/4
zeta1Step <- sd( Q1MCMC$zeta1 )/4
zeta2Step <- sd( Q1MCMC$zeta2 )/4
rho1Step <- apply( Q1MCMC$rho1, 2, sd )/4
kappa1Step <- sd( Q1MCMC$kappa1)/4
phi1Step <- apply( Q1MCMC$phi1, 2, sd )/4




#########  This will create a R database and save it.  When you load it everything will be there.

save(alpha1, beta1, gamma1, gamma2, eta1, rho1, phi1, zeta1, zeta2, kappa1,
     n1Q, chpt1, rchpt1, gchpt1, chpt2, alpha1Step,
     beta1Step, gamma1Step, gamma2Step, eta1Step, rho1Step, 
     phi1Step, zeta1Step, zeta2Step, kappa1Step, file = filename8_SVEIRD)


########  You can add a load statement above and it will pull all these in so you don't need to start all over.
#########  Why not just use a save function?  This will create a R database and save it.  When you load it everything
#########  will be there.

Q1MCMC20_SV <- Q1MCMC
save( Q1MCMC20_SV, file = "Q1MCMC20_SV.Rdat" )

#################################################################################
########  You can add a load statement above and it will pull all these in so you don't need to start all over.

ISVPred1m <- apply( Q1MCMC20_SV$IPred1, 2, quantile, c(0.5,0.25,0.975) )
IISVPred1m <- apply( Q1MCMC20_SV$IIPred1, 2, quantile, c(0.5,0.25,0.975) )

RESVPred1m <- apply( Q1MCMC20_SV$REPred1, 2, quantile, c(0.5,0.25,0.975) )
RISVPred1m <- apply( Q1MCMC20_SV$RIPred1, 2, quantile, c(0.5,0.25,0.975) )
RRSVPred1m <- apply( Q1MCMC20_SV$RRPred1, 2, quantile, c(0.5,0.25,0.975) )

DSVPred1m <- apply( Q1MCMC20_SV$DPred1, 2, quantile, c(0.5,0.25,0.975) )

VSVPred1m <- apply( Q1MCMC20_SV$VPred1, 2, quantile, c(0.5,0.25,0.975) )

###Extracting the "median" of all states (Infected, reinfected, deaths, vaccinated, recovered and re-recovered)
meanSV1 <- apply(Q1MCMC20_SV$IPred1, 2, function(x) quantile(x, 0.5)) #Infected 

meanSV2 <- apply(Q1MCMC20_SV$IIPred1, 2, function(x) quantile(x, 0.5)) #Reinfected
meanSV3 <- apply(Q1MCMC20_SV$RIPred1, 2, function(x) quantile(x, 0.5)) #Recovered
meanSV4 <- apply(Q1MCMC20_SV$RRPred1, 2, function(x) quantile(x, 0.5)) #Re-recovered
meanSV5 <- apply(Q1MCMC20_SV$DPred1, 2, function(x) quantile(x, 0.5)) #Deaths
meanSV6 <- apply(Q1MCMC20_SV$VPred1, 2, function(x) quantile(x, 0.5))# Vaccinated


#####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

################### Checking the model fits #########################
# Basic Model Fit Statistics
ISV1Error <- sum( (Q1$AdjInfect - meanSV1)^2 )
IISV1Error <- sum( (Q1$CReinfection - meanSV2)^2 )

RISV1Error <- sum( (Q1$CRecovered - meanSV3)^2 )
RRSV1Error <- sum( (Q1$CReRecover - meanSV4)^2 )
DSV1Error <- sum( (Q1$CDeaths - meanSV5)^2 )
VSV1Error <- sum( (Q1$CVaccine - meanSV6)^2 )
TotErrorSV1 <- (n6-1)*(sd( Q1$AdjInfect)^2 + sd(Q1$CReinfection)^2 +
                         sd( Q1$CRecovered )^2
                       + sd( Q1$CReRecover)^2 + sd( Q1$CDeaths )^2 + 
                         sd( Q1$CVaccine )^2 )
TotErrorSV1

# Pseudo R2

PsuedoR2.SV1 <- 1 - (ISV1Error + IISV1Error + RISV1Error + RRSV1Error + DSV1Error + VSV1Error)/TotErrorSV1


PsuedoR2.SV1

# Other analysis
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha1 <- Q1MCMC20_SV$alpha1
beta1 <- Q1MCMC20_SV$beta1
betaI1 <- Q1MCMC20_SV$betaI1
gamma1 <- Q1MCMC20_SV$gamma1
gamma2 <- Q1MCMC20_SV$gamma2
phi1 <- Q1MCMC20_SV$phi1
zeta1 <- Q1MCMC20_SV$zeta1
zeta2 <- Q1MCMC20_SV$zeta2
rho1 <- Q1MCMC20_SV$rho1
kappa1 <- Q1MCMC20_SV$kappa1
eta1 <-Q1MCMC20_SV$eta1
#
#####%%%%%%%%%  Estimating the means, standard deviation and quantiles of the parameters %%%%%%%%%%%%%%
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  Mean median, sd, and quantile for alpha1 
mean_alpha1 <- apply(alpha1, 2, mean);mean_alpha1
sd_alpha1 <- apply(alpha1, 2, sd);sd_alpha1
median_alpha1 <- apply(alpha1, 2, median);median_alpha1
quantile_alpha1 <- apply(alpha1, 2, quantile, c(0.025,0.5,0.975));quantile_alpha1

###  Mean median, sd, and quantile for beta1 
mean(beta1)
sd(beta1)
median(beta1)
quantile(beta1, c(0.025,0.5,0.975))

###  Mean median, sd, and quantile for Kappa
mean(kappa1)
median(kappa1)
sd(kappa1)
quantile(kappa1, c(0.025,0.5,0.975))

###  Mean median, sd, and quantile for eta
mean(eta1)
sd(eta1)
median(eta1)
quantile(eta1, c(0.025,0.5,0.975))

###  Mean median, sd, and quantile for phi

mean_phi1 <- apply(phi1, 2, mean);mean_phi1
median_phi1 <- apply(phi1, 2, median);median_phi1
quantile_phi1 <- apply(phi1, 2, quantile, c(0.025,0.5,0.975));quantile_phi1

###  Mean median, sd, and quantile for rho
mean(rho1)
sd(rho1)
median(rho1)
quantile(rho1, c(0.025,0.5,0.975))

###  Mean median, sd, and quantile for zeta1
mean(zeta1)
sd(zeta1)
median(zeta1)
quantile(zeta1, c(0.025,0.5,0.975))

###  Mean median, sd, and quantile for zeta2
mean(zeta2)
sd(zeta2)
median(zeta2)
quantile(zeta2, c(0.025,0.5,0.975))

###  Mean median, sd, and quantile for gamma2
mean(gamma2)
sd(gamma2)
median(gamma2)
quantile(gamma2, c(0.025,0.5,0.975))

# ###  Mean median, sd, and quantile for gamma1
mean_gamma1 <- apply(gamma1, 2, mean);mean_gamma1
sd_gamma1 <- apply(gamma1, 2, sd);sd_gamma1
median_gamma1 <- apply(gamma1, 2, median);median_gamma1
quantile_gamma1 <- apply(gamma1, 2, quantile, c(0.025,0.5,0.975));quantile_gamma1

############################  RESULTS #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#         See the output file for the results
####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

############# Plotting the Posterior Predictive bands 

pdf("infectedSV.pdf")
n6 <- n1
plot( 1:n6, Q1$AdjInfect, type = "l" , col = "red", lwd = 3, lty = 2, ylim = c(min(Q1$AdjInfect), max(Q1$AdjInfect)+400), ylab = " Actively Infected",
      xlab = "Days since 28Feb2020", main="Infected- Model 2")
lines(1:n6, IPredSV[2,], lty = 1, lwd=2, col = "darkred") #median
# Add confidence bands
polygon(c(1:n6, rev(1:n6)), c(IPredSV[1,], rev(IPredSV[3,])), col = rgb(0, 0, 0, alpha = 0.3))

lines( 1:n6, IPredSV[1,], lty = 2, lwd=1 ) #lower
lines( 1:n6, IPredSV[3,], lty = 2, lwd=1 ) #upper

# Add a legend
legend("topright", legend = c("Actual", "Fitted"), lwd = c(3,2), lty = c(2, 1), col = c("red", "darkred"))

dev.off()



pdf("RecoverSV.pdf")          
plot( 1:n6, Q1$CRecovered, type = "l", lty = 2, col = "green4", lwd = 3, ylab = " Cumulative Recovered",
      xlab = "Days since 6Mar2020", main="Recovered - Model 2")
lines(1:n6, RIPredSV[2,], lty=1, lwd = 2, col = "darkgreen",) #median
# Add confidence bands
polygon(c(1:n6, rev(1:n6)), c(RIPredSV[1,], rev(RIPredSV[3,])), col = rgb(0, 0, 0, alpha = 0.3), border = NULL)

lines( 1:n6, RIPredSV[1,], lty = 2) #lower
lines( 1:n6, RIPredSV[3,], lty = 2) # uppe
# Create a legend
legend("topleft", legend = c("Actual", "Fitted"), lty = c(2, 1), lwd = c(3, 2), col = c("green4", "darkgreen"))

dev.off()


#####################
pdf("deathSV.pdf")
plot( 1:n6, Q1$CDeaths, type = "l", lty = 2, col = "grey0", lwd = 3, ylim = c(0, 650),
      xlab = "Days since 11Mar2020",pch = 1, ylab = "Cumulative Deaths", main = "Deaths - Model 2")
lines(1:n6, DPredSV[2,], lty=1, lwd = 2, col = "black") # median
# Add confidence bands
polygon(c(1:n6, rev(1:n6)), c(DPredSV[1,], rev(DPredSV[3,])), col = rgb(0, 0, 0, alpha = 0.5), border = NULL)

lines( 1:n6, DPredSV[1,], lty = 2) #lower
lines( 1:n6, DPredSV[3,],lty = 2) #upper

legend("topleft", legend = c("Actual", "Fitted"), lty = c(2, 1), lwd = c(3, 2), col=c("grey0", "black"))

dev.off()



###################
pdf("ReinfectSV.pdf")
plot( 1:n6, Q1$CReinfection, type = "l", lty = 2, lwd = 3, col = "brown4", ylim= c(0, 45), ylab = " Cumulative Reinfected",
      xlab = "Days since 23Jun2020", main="Number of Reinfections - model 2")
lines(1:n6, IIPredSV[2,], lty=1, lwd = 2, col="orange3") #median
polygon(c(1:n6, rev(1:n6)), c(IIPredSV[1,], rev(IIPredSV[3,])), col = rgb(0, 0, 0, alpha = 0.6), border = NULL)

lines( 1:n6, IIPredSV[1,], lty = 2) #lower
lines( 1:n6, IIPredSV[3,], lty = 2) #upper

legend("topright", legend = c("Actual", "Fitted"), lty = c(2,1), lwd = c(3, 2), col = c("brown4", "orange3"))

dev.off()

#############################
pdf("Re-recoverSV.pdf")   

plot( 1:n6, Q1$CReRecovered, type = "l", col = "darkblue", lty = 2, lwd = 3,  ylab = " Cumulative Re-Recovered",
      xlab = "Days since 30Jun2020", main="Those who recovered after reinfections - Model 2")
lines( 1:n6, RRPredSV[2,],lty = 1, lwd = 2, col = "blue1") #median
polygon(c(1:n6, rev(1:n6)), c(RRPredSV[1,], rev(RRPredSV[3,])), col = rgb(0, 0, 0, alpha = 0.3), border = NULL)
lines( 1:n6, RRPredSV[1,],lty = 2) #lower
lines( 1:n6, RRPredSV[3,],lty = 2) #upper
legend("topleft", legend = c("Actual", "Fitted"), lty = c(2,1), lwd = c(3, 2), col = c("darkblue", "blue1"))

dev.off()

##################
pdf("vaccineSV.pdf")
plot( 1:n6, Q1$CVaccine, type = "l",
      xlab = "Days since 31Dec2020", lwd = 2, lty = 2, col = "darkviolet",  ylab = " Cumulative number of those Vaccinated", main = "vaccinated - Model 2 ")
lines( 1:n6, VPredSV[2,],lty = 1, lwd = 2, col = "darkviolet") #median
polygon(c(1:n6, rev(1:n6)), c(VPredSV[1,], rev(VPredSV[3,])), col = rgb(0, 0, 0, alpha = 0.5), border = NULL)
lines( 1:n6, VPredSV[1,],lty = 2) #lower
lines( 1:n6, VPredSV[3,],lty = 2) #upper
legend("topleft", legend = c("Actual", "Fitted"), lty = c(2,1), lwd = c(2, 2), col = c("darkviolet", "darkviolet"))

dev.off()

######### The end of parameter estimation and model fitting for model 2 ##########################
