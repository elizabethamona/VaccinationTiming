######################################################
#
### This code includes the parameter estimation and model fitting to the data.
### To run this code, you need to first load the "MH55.Cpp" and "MH55.R" files
####Warning: carefully run each section of this code to reproduce the results, skipping any section might result in an error message.

##############################################################################
################################################################################

#Load the required library
library( Rcpp )

###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read in the functions we will use %%%%%%%%%%%%%%%%%%%%%%%%%%%%
sourceCpp( "MH55.cpp" ) ### C++ file

source("MH55.R") #R file

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


####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Starting Values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Data for Qatar
S0 <- 2695122  #(S_1): First Susceptible  #total population as of 2022 
SI0 <- 2695122/6  #second susceptible (S_2)
E0 <- 5 #Exposed
I0 <- 1 #infected
II0 <- 0 #Reinfected
#Recovered
RE0 <- 0
RI0 <- 0
#Re-recovered
RR0 <- 0
D0 <- 0  #death
##vaccination
V0 <- 0


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
beta1 <-  0.07072258
eta1 <- 0.000202525
gamma2 <- 1/7

###########################################################
beta1Step <- 0.0001
gamma2Step <- 0.0001
eta1Step <- 0.0000001
zeta1Step <- 0.00001
zeta2Step <- 0.000001


#rchpt1 <- c( 1, 390)
#phi1 <-  0.00002125089
#rho1 <- c(0, 0.003850046)  
###efficacy: so we chose 94% (1-0.94=0.06)
# zeta1 <- 1/90*0.06   #90-3 days
# zeta2 <- 1/90*0.04

zeta1 <- 0.0003300381   
zeta2 <- 1/90
kappa1 <- 0.94

zeta1Step <- 0.000001
zeta2Step <- 0.000001
kappa1Step <- 1.99e-03

rchpt1 <- c( 1, 390)  # Time at which vaccine1 was deployed
rho1 <- c(0, 0.01)  #1/100 persons
rho1Step <- c(0,0.0001)
#gamma change point
gchpt1 <- c(6, 60, 350, 428)
gamma1 <- c(0.04589753, 0.14013582, 0.13010453, 0.11339254)
gamma1Step <- rep(1.00e-4, length(gamma1) )


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chpt1 <- c(  20, 40, 80, 95, 160)
chpt1 <- c(chpt1, 170, 190, 250, 325, 360, 448, 500)

# These are good values to start the analysis with. Also called initial values.
alpha1 <- c(2.420616e-07, 8.199343e-08, 9.372360e-08, 7.839590e-08, 1.106528e-07)
alpha1 <- c(alpha1, 9.132672e-08, 8.826232e-08, 1.005170e-07, 1.133644e-07, 1.293249e-07, 2.646762e-07, 5.044460e-07)

alpha1Step <- c( 1.00e-11, 1.00e-11, 1.00e-11, 1.00e-11, 1.00e-11)
alpha1Step <- c( alpha1Step, 1.00e-11, 1.00e-11, 1.00e-11, 1.00e-11, 1.00e-11, 1.00e-11,  1.00e-11)


#phi1 change point. These are also initial values
chpt2 <- c( 117, 350, 450, 500)
phi1 <- c(0.0002125089, 0.00015668648, 0.0012775134, 0.0022775134)
#phi1 <- c(0.0010000000, 0.012447373, 0.0012775134, 0.012775134)
phi1Step <- c( 1.00e-06, 1.00e-6, 1.00e-6, 1.00e-6)

####### Create start step values
filename8_new <- "StartStepValuesMC8_new" # These are the values that have converged and used for the analysis
#load(filename8_new) # Loading these values allows you to reproduce my start values. But you can leave it commented if you want to generate your own startvalues. 



tic <- Sys.time()
nsamp1 = 50000 # You can change this number to suit your need. I generated 50,000 samples after I ensured my analysis has converged. Feel free to change the value to a smaller sample size or bigger, depending on your goal.
nMCMC1 <- nsamp1
Q1MCMC <- MCMCSEIRD8VPred(data3, S0, E0, I0, II0,
                          RE0, RI0, SI0, RR0, D0, V0, alpha1, beta1, gamma1,
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
# Trace plots to ensure convergence 
#
#######################################################################

for( i in 1:length(alpha1) ){
  plot( Q1MCMC$alpha1[,5], ylab="alpha",type = "l" )
}

plot( Q1MCMC$beta1, xlab="samples" , ylab="beta",type = "l" )

for( i in 1:length( gamma1 ) ){
  plot( Q1MCMC$gamma1[,4] , ylab="gamma", type = "l" )
}

# pdf("Trace Deaths.pdf")
plot( Q1MCMC$eta1, ylab="Deaths", type = "l" )
# dev.off()

# plot( Q1MCMC$rho1[,1], type = "l" )
#pdf("Trace rho.pdf")
plot( Q1MCMC$rho1[,2], ylab="rho", type = "l" )
#dev.off()

for( i in 1:length( phi1 ) ){
  #pdf("Trace phi1.pdf")
  plot( Q1MCMC$phi1[,1] , ylab="phi", type = "l" )
  #dev.off()
}
# pdf("Trace zeta1.pdf")
plot( Q1MCMC$zeta1, ylab="zeta1", type = "l" )
# dev.off()
# 
# pdf("Trace zeta2.pdf")
plot( Q1MCMC$zeta2, ylab="zeta2", type = "l" )
# dev.off()

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
     phi1Step, zeta1Step, zeta2Step, kappa1Step, file = filename8_new)



########  You can add a load statement above and it will pull all these in so you don't need to start all over.
#########  Why not just use a save function?  This will create a R database and save it.  When you load it everything
#########  will be there.

save(alpha1, beta1, gamma1, gamma2, eta1, rho1, phi1, zeta1, zeta2, kappa1,
     n1Q, chpt1, rchpt1, gchpt1, chpt2, alpha1Step,
     beta1Step, gamma1Step, gamma2Step, eta1Step, rho1Step, 
     phi1Step, zeta1Step, zeta2Step, kappa1Step, file = filename8_new)



########  You can add a load statement above and it will pull all these in so you don't need to start all over.

Q1MCMC20 <- Q1MCMC
save( Q1MCMC20, file = "Q1MCMC20.Rdat" )
load("Q1MCMC20.Rdat")
##This contains all the predicted values and the log.lik values for all the 50,000 samples


#########%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
########## generating the quantiles of the predicted value for Infected, Reinfected, Recovered, Deaths, Re-recovered and Vaccinated.

IPred1m <- apply( Q1MCMC20$IPred1, 2, quantile, c(0.5,0.25,0.975) ) #Infected
IIPred1m <- apply( Q1MCMC20$IIPred1, 2, quantile, c(0.5,0.25,0.975) )#Reinfected

RIPred1m <- apply( Q1MCMC20$RIPred1, 2, quantile, c(0.5,0.25,0.975) ) #Recovered
RRPred1m <- apply( Q1MCMC20$RRPred1, 2, quantile, c(0.5,0.25,0.975) ) #Re-recovered

DPred1m <- apply( Q1MCMC20$DPred1, 2, quantile, c(0.5,0.25,0.975) ) #Deaths

VPred1m <- apply( Q1MCMC20$VPred1, 2, quantile, c(0.5,0.25,0.975) ) #Vaccinated

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##### Checking the model fit 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Basic Model Fit Statistics
IError <- sum( (Q1$AdjInfect - IPred1m[1,])^2 );IError # error of Infected
IIError <- sum( (Q1$CReinfection - IIPred1m[1,])^2 );IIError  # Reinfected

RIError <- sum( (Q1$CRecovered - RIPred1m[1,])^2 );RIError# Recovered
RRError <- sum( (Q1$CReRecover - RRPred1m[1,])^2 );RRError #Re-recovered

DError <- sum( (Q1$CDeaths - DPred1m[1,])^2 );DError  #deaths
VError <- sum( (Q1$CVaccine - VPred1m[1,])^2 );VError #vaccinated
TotError1 <- (n1Q-1)*(sd( Q1$AdjInfect)^2 + sd(Q1$CReinfection)^2 +       # Total error
                        sd( Q1$CRecovered )^2
                      + sd( Q1$CReRecover)^2 + sd( Q1$CDeaths )^2 + 
                        sd( Q1$CVaccine )^2 )
TotError1

# Pseudo R2

PsuedoR2 <- 1 - (IError + IIError + RIError + RRError + DError + VError)/TotError1


PsuedoR2

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####################### OTHER Statistical ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Renaming the parameters 
alpha1 <- Q1MCMC20$alpha1
beta1 <- Q1MCMC20$beta1
betaI1 <- Q1MCMC20$betaI1
gamma1 <- Q1MCMC20$gamma1
gamma2 <- Q1MCMC20$gamma2
phi1 <- Q1MCMC20$phi1
zeta1 <- Q1MCMC20$zeta1
zeta2 <- Q1MCMC20$zeta2
rho1 <- Q1MCMC20$rho1
kappa1 <- Q1MCMC20$kappa1
eta1 <- Q1MCMC20$eta1

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

IPred1m <- apply( Q1MCMC20 $IPred1, 2, quantile, c(0.5,0.25,0.975) )
IIPred1m <- apply( Q1MCMC20$IIPred1, 2, quantile, c(0.5,0.25,0.975) )

REPred1m <- apply( Q1MCMC20$REPred1, 2, quantile, c(0.5,0.25,0.975) )
RIPred1m <- apply( Q1MCMC20$RIPred1, 2, quantile, c(0.5,0.25,0.975) )
RRPred1m <- apply( Q1MCMC20$RRPred1, 2, quantile, c(0.5,0.25,0.975) )

DPred1m <- apply( Q1MCMC20$DPred1, 2, quantile, c(0.5,0.25,0.975) )

VPred1m <- apply( Q1MCMC20$VPred1, 2, quantile, c(0.5,0.25,0.975) )

######################## RESULTS OF THE POSTERIOR with predictive bands

pdf("infected.pdf")
n1Q <- n1
plot( 1:n1Q, Q1$AdjInfect, type = "l",  lty = 2, col = "red", lwd = 3, ylim = c(min(Q1$AdjInfect), max(Q1$AdjInfect)+400), ylab = " Actively Infected",
      xlab = "Days since 28Feb2020", main="Infected- Model 1")
lines( 1:n1Q, IPred1m[1,], lwd = 2, col = "darkred") #median
# Add confidence bands
polygon(c(1:n1Q, rev(1:n1Q)), c(IPred1m[2,], rev(IPred1m[3,])), col = rgb(0, 0, 0, alpha = 0.5), border = NULL)

lines( 1:n1Q, IPred1m[2,],lty = 2 ) #lower
lines( 1:n1Q, IPred1m[3,], lty = 2 ) #upper
# Add a legend
legend("topright", legend = c("Actual", "Fitted"), lty = c(2, 1), lwd = c(3, 2), col= c("red", "darkred"))
dev.off()
#######################

pdf("Recover.pdf")          
plot( 1:n1Q, Q1$CRecovered, type = "l", lwd=3, lty = 2,  col = "green4", ylab = " Cumulative Recovered",
      xlab = "Days since 6Mar2020", main="Recovered - Model 1")
lines( 1:n1Q, RIPred1m[1,], lwd = 2, col = "darkgreen", lty = 1) #median
# Add confidence bands
polygon(c(1:n1Q, rev(1:n1Q)), c(RIPred1m[2,], rev(RIPred1m[3,])), col = rgb(0, 0, 0, alpha = 0.5), border = NULL)
lines( 1:n1Q, RIPred1m[2,], lty = 2) #lower
lines( 1:n1Q, RIPred1m[3,], lty = 2) #upper
# Create a legend
legend("topleft", legend = c("Actual", "Fitted"), lty = c(2, 1), lwd = c(3, 2), col = c("green4", "darkgreen"))
dev.off()

#####################
pdf("death.pdf")
plot( 1:n1Q, Q1$CDeaths, type = "l", lty = 2, 
      xlab = "Days since 11Mar2020", col = "grey0", lwd = 3, ylab = "Cumulative Deaths", main = "Deaths - Model 1")
lines( 1:n1Q, DPred1m[1,], lwd = 2, col = "black")
# Add confidence bands
polygon(c(1:n1Q, rev(1:n1Q)), c(DPred1m[2,], rev(DPred1m[3,])), col = rgb(0, 0, 0, alpha = 0.5), border = NULL)

lines( 1:n1Q, DPred1m[2,], lty = 2)
lines( 1:n1Q, DPred1m[3,],lty = 2)
legend("topleft", legend = c("Actual", "Fitted"), lty = c(2,1), lwd = c(3, 2),col=c("grey0", "black"))

dev.off()

##################
pdf("vaccine.pdf")
plot( 1:n1Q, Q1$CVaccine, type = "l", lty =2, 
      xlab = "Days since 31Dec2020", lwd = 3, col = "darkviolet", ylab = " Cumulative number of those Vaccinated", main = "vaccinated - Model 1 ")
lines( 1:n1Q, VPred1m[1,], lty=1, lwd = 2, col = "darkviolet")
polygon(c(1:n1Q, rev(1:n1Q)), c(VPred1m[2,], rev(VPred1m[3,])), col = rgb(0, 0, 0, alpha = 0.6), border = NULL)
lines( 1:n1Q, VPred1m[2,], lty = 2)
lines( 1:n1Q, VPred1m[3,], lty = 2)
legend("topleft", legend = c("Actual", "Fitted"), lty = c(2, 1), lwd = c(3, 2), col = c("darkviolet", "darkviolet"))

dev.off()

###################
pdf("Reinfect.pdf")
plot( 1:n1Q, Q1$CReinfection, type = "l", lty = 2, lwd = 3, col = "brown4", ylim = c(0, 45), ylab = " Cumulative Reinfected",
      xlab = "Days since 23Jun2020", main="Number of Reinfections - model 1")
lines( 1:n1Q, IIPred1m[1,], lty = 1, lwd = 2, col="orange3")
polygon(c(1:n1Q, rev(1:n1Q)), c(IIPred1m[2,], rev(IIPred1m[3,])), col = rgb(0, 0, 0, alpha = 0.6), border = NULL)
lines( 1:n1Q, IIPred1m[2,], lty = 2)
lines( 1:n1Q, IIPred1m[3,], lty = 2)
legend("topright", legend = c("Actual", "Fitted"), lty = c(2,1), lwd = c(3, 2), col = c("brown4", "orange3"))

dev.off()

#############################
pdf("Re-recover.pdf")   
plot( 1:n1Q, Q1$CReRecovered, type = "l", lty=2, col = "darkblue", lwd = 3, ylab = " Cumulative Re-Recovered",
      xlab = "Days since 30Jun2020", main="Those who recovered after reinfections - Model 1")
lines( 1:n1Q, RRPred1m[1,], lty=1, lwd = 2, col = "blue1")
polygon(c(1:n1Q, rev(1:n1Q)), c(RRPred1m[2,], rev(RRPred1m[3,])), col = rgb(0, 0, 0, alpha = 0.5), border = NULL)

lines( 1:n1Q, RRPred1m[2,],lty = 2)
lines( 1:n1Q, RRPred1m[3,],lty = 2)
legend("topleft", legend = c("Actual", "Fitted"),  lty = c(2, 1), lwd = c(3, 2), col = c("darkblue", "blue1"))
dev.off()
#######################################################################################################

######### The end of parameter estimation and model fitting for model 1 ####################
