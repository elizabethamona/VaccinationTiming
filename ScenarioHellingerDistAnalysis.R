################################################################################
#### This code provides the scenario Analysis for model 1 (S1EI1R1DVI1R2 model)
#There are 6 scenario analysis performed here:
##    Scenario 1: Distribution of Infected, Reinfected and Death with 94% vaccine efficacy (day vaccine was available)
#     Scenario 2: Distribution of Infected, Reinfected and Death with 100% vaccine efficacy (day vaccine was available)
#     Scenario 3: Distribution of Infected, Reinfected and Death with 94% vaccine efficacy (Early vaccination date (day 200))
#     Scenario 4: Distribution of Infected, Reinfected and Death with 94% vaccine efficacy (Late vaccination date (day 450))
#     Scenario 5: Distribution of Infected, Reinfected and Death with 100% vaccine efficacy (Early vaccination date (day 200))
#     Scenario 6: Distribution of Infected, Reinfected and Death with 100% vaccine efficacy (Late vaccination date (day 450))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####Note: Run each section and make sure you load all the results in the form of Rdat file before any analysis.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

####### Loading all results 
load("Q1MCMC20.Rdat")  #Model 1 result
load("Q1MCMC20_SV.Rdat") # Model 2 result
load("Q1MCMC20_555.Rdat") # Scenario 1 result for model 1
load("Q1MCMC20_k1.Rdat") #  Scenario 2 result for model 1
load("Q1MCMC20_555E.Rdat") # Scenario 3 result for model 1
load("Q1MCMC20_555L.Rdat") # Scenario 4 result for model 1
load("Q1MCMC20_555E100.Rdat") # Scenario 5 result for model 1
load("Q1MCMC20_555L100.Rdat") # Scenario 6 result for model 1

Q1MCMC20$Post2  ## -logLik for Model 1
Q1MCMC20_SV$Post2  # -logLik for Model 2



### Renaming for plotting purposes
#Infected 
InfectI <- Q1MCMC20_555$IPred1  #94% eff1cacy
Infectk1 <- Q1MCMC20_k1$IPred1  #100% efficacy
InfectE <- Q1MCMC20_555E$IPred1  # Early vaccine (day 200)
InfectL <- Q1MCMC20_555L$IPred1  # late vaccine (day 450)
InfectE100 <- Q1MCMC20_555E100$IPred1  # Early vaccine (day 200)
InfectL100 <- Q1MCMC20_555L100$IPred1  # late vaccine (day 450)

#Reinfected
ReinfectI <- Q1MCMC20_555$IIPred1
Reinfectk1 <- Q1MCMC20_k1$IIPred1
ReinfectE <- Q1MCMC20_555E$IIPred1
ReinfectL <- Q1MCMC20_555L$IIPred1
ReinfectE100 <- Q1MCMC20_555E100$IIPred1
ReinfectL100 <- Q1MCMC20_555L100$IIPred1

#Deaths
DeathsI <- Q1MCMC20_555$DPred1
Deathsk1 <- Q1MCMC20_k1$DPred1
DeathsE <- Q1MCMC20_555E$DPred1
DeathsL <- Q1MCMC20_555L$DPred1
DeathsE100 <- Q1MCMC20_555E100$DPred1
DeathsL100 <- Q1MCMC20_555L100$DPred1


## Cumulative Infected at day 540
InfectI.540 <- cumsum(Q1MCMC20_555$IPred1[,540])
Infectk1.540 <- cumsum(Q1MCMC20_k1$IPred1[,540])
InfectE.540 <- cumsum(Q1MCMC20_555E$IPred1[,540])
InfectL.540 <- cumsum(Q1MCMC20_555L$IPred1[,540])
InfectE100.540 <- cumsum(Q1MCMC20_555E100$IPred1[,540])
InfectL100.540 <- cumsum(Q1MCMC20_555L100$IPred1[,540])



## Cumulative Reinfected at day 540
ReinfectI.540 <- cumsum(Q1MCMC20_555$IIPred1[,540])
Reinfectk1.540 <- cumsum(Q1MCMC20_k1$IIPred1[,540])
ReinfectE.540 <- cumsum(Q1MCMC20_555E$IIPred1[,540])
ReinfectL.540 <- cumsum(Q1MCMC20_555L$IIPred1[,540])
ReinfectE100.540 <- cumsum(Q1MCMC20_555E100$IIPred1[,540])
ReinfectL100.540 <- cumsum(Q1MCMC20_555L100$IIPred1[,540])



##Cumulative Infected at day 540
DeathsI.540 <- Q1MCMC20_555$DPred1[,540]
Deathsk1.540 <- Q1MCMC20_k1$DPred1[,540]
DeathsE.540 <- Q1MCMC20_555E$DPred1[,540]
DeathsL.540 <- Q1MCMC20_555L$DPred1[,540]
DeathsE100.540 <- Q1MCMC20_555E100$DPred1[,540]
DeathsL100.540 <- Q1MCMC20_555L100$DPred1[,540]


################################################################################
######%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### Plotting the densities using "Guassian kernel"

####### For individual:
#Infected
den1 <- density(InfectI.540, kernel = "gaussian", from = min(InfectI.540), to = max(InfectI.540), window=kernel, n=50000)

den2 <- density(Infectk1.540,  kernel = "gaussian", from = min(InfectI.540), to = max(InfectI.540), window=kernel, n=50000)

den4 <- density(InfectE.540, kernel = "gaussian", from = min(InfectE.540), to = max(InfectE.540), window=kernel, n=50000)

den5 <- density(InfectL.540, kernel = "gaussian", from = min(InfectL.540), to = max(InfectL.540), window=kernel, n=50000)

den6 <- density(InfectE100.540, kernel = "gaussian", from = min(InfectE100.540), to = max(InfectE100.540), window=kernel, n=50000)

den7 <- density(InfectL100.540, kernel = "gaussian", from = min(InfectL100.540), to = max(InfectL100.540), window=kernel, n=50000)

##### Reinfected 
den11 <- density(ReinfectI.540)

den22 <- density(Reinfectk1.540)

den44 <- density(ReinfectE.540)

den55 <- density(ReinfectL.540)

den66 <- density(ReinfectE100.540)

den77 <- density(ReinfectL100.540)


#### Deaths

den1.1 <- density(DeathsI.540)

den2.2 <- density(Deathsk1.540)

den4.4 <- density(DeathsE.540)

den5.5 <- density(DeathsL.540)

den6.6 <- density(DeathsE100.540)

den7.7 <- density(DeathsL100.540)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Scenario plots for the Cumulative Infected at day 540. See Figure 5(a)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pdf("Infected.all.pdf")
plot(den1, xlim= c(0, 27000000), ylim=c(0e+00, 4e-07), xlab = substitute(paste(bold("Cumulative Infected"))), ylab = substitute(paste(bold("Density"))), col = "blue", lty = 2, lwd = 2, 
     main = "Scenario Comparison of Infected - Model 1")
lines(den2, col = "black" , lty = 2, lwd = 2)
lines(den4, col = "red3" , lty = 2, lwd = 2)
lines(den5, col = "violet", lty = 16, lwd = 2)
lines(den6, col = "darkturquoise" , lty = 4, lwd = 2)
lines(den7, col = "tan2", lty = 16, lwd = 2)
#lines(den8, col = "yellowgreen" , lty = 2, lwd = 2)

linetypes = c(2,1,4, 16, 4, 16, 32)
names(linetypes) = c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5", "Scenario 6", lwd = 2)


legend("topright", c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5", "Scenario 6"),
       cex=0.85, col =c("blue","black","red3","violet", "darkturquoise","tan2"), lty = linetypes, lwd = 2)

dev.off()


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#      HELLINGER DISTANCE CALCULATION OF THE DENSITIES OF DEATHS FOR ALL 6 SCENARIOS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Note: Step-by-step description of how the Hellinger distance is calculated is shown below:

###################################
# For Scenario 1 and Scenario 2:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Step 1: Calculate Density Estimates for Both Posteriors
den1 # density p
den2 # density q

# Step 2: Create Common Breaks:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Important: -To calculate the Hellinger distance, you should use the same breaks (bins) for both density estimates.
#            -You can do this by determining the common range of values and creating a common set of breaks for both densities.
#            -The range function helps you find the range of values, and the seq function can be used to create breaks:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Find the common range of the x values of the densities
common_range <- range(den1$x, den2$x)

# Define the number of breaks (you can adjust this as needed)
num_breaks <- 100

# Choosing equally spaced points along the ranges
z.points <- seq(common_range[1], common_range[2], length.out = num_breaks)

#Step 3: Calculate Densities with Common Breaks
#Now, you need to re-calculate the densities using the common breaks:

# Calculate density estimates with common breaks
den1_common <- density(InfectI.540, from = common_range[1], to = common_range[2], n = num_breaks) #94%
den2_common <- density(Infectk1.540, from = common_range[1], to = common_range[2], n = num_breaks) #100%
hist(Infectk1.540, breaks= 100)

hellinger_distance <- (sum((sqrt(den1_common$y) - sqrt(den2_common$y))^2)*(z.points[2]-z.points[1]))^1/2  #use this one

# Print the Hellinger distance
cat("Hellinger Distance:", hellinger_distance, "\n")
# Ans: Hellinger Distance: 2.034972e-05  (same as we have in the paper)


###################################################
#Scenarios 1 and 3
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Find the common range of values
common_range.14 <- range(den1$x, den4$x)

# Define the number of breaks (you can adjust this)
num_breaks <- 100

# Create common breaks for both densities
z.points <- seq(common_range.14[1], common_range.14[2], length.out = num_breaks)

#Step 3: Calculate Densities with Common Breaks
#Now, you need to re-calculate the densities using the common breaks:

# Calculate density estimates with common breaks
den1_common <- density(InfectI.540, from = common_range.14[1], to = common_range.14[2], n = num_breaks) #94%
den4_common <- density(InfectE.540, from = common_range.14[1], to = common_range.14[2], n = num_breaks) #100%

#Step 4: Calculate Hellinger Distance

hellinger_distance <- (sum((sqrt(den1_common$y) - sqrt(den4_common$y))^2)*(z.points[2]-z.points[1]))^1/2  #use this one

# Print the Hellinger distance
cat("Hellinger Distance:", hellinger_distance, "\n")
#Hellinger Distance: 0.3406251



###################################################
#Early (94 and 100 efficacy): Scenarios 3 and 5
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Find the common range of values
common_range.46 <- range(den4$x, den6$x)

# Define the numbr of breaks (you can adjust this)
num_breaks <- 100

# Create common breaks for both densities
z.points <- seq(common_range.46[1], common_range.46[2], length.out = num_breaks)

# Calculate density estimates with common breaks
den4_common <- density(InfectE.540, from = common_range.46[1], to = common_range.46[2], n = num_breaks) #94%
den6_common <- density(InfectE100.540, from = common_range.46[1], to = common_range.46[2], n = num_breaks) #100%

#Step 4: Calculate Hellinger Distance
hellinger_distance <- (sum((sqrt(den4_common$y) - sqrt(den6_common$y))^2)*(z.points[2]-z.points[1]))^1/2  #use this one

cat("Hellinger Distance:", hellinger_distance, "\n")

#Hellinger Distance: 2.902628e-05


#############################################
#Late (94 and 100 efficacy): Scenarios 4 and 6
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Find the common range of values
common_range.57 <- range(den5$x, den7$x)

# Define the number of breaks (you can adjust this)
num_breaks <- 100

# Create common breaks for both densities
z.points <- seq(common_range.57[1], common_range.57[2], length.out = num_breaks)

# Calculate density estimates with common breaks
den5_common <- density(InfectL.540, from = common_range.57[1], to = common_range.57[2], n = num_breaks) #94%
den7_common <- density(InfectL100.540, from = common_range.57[1], to = common_range.57[2], n = num_breaks) #100%


#Step 4: Calculate Hellinger Distance
hellinger_distance <- (sum((sqrt(den5_common$y) - sqrt(den7_common$y))^2)*(z.points[2]-z.points[1]))^1/2  #use this one

cat("Hellinger Distance:", hellinger_distance, "\n")
####Hellinger Distance: 0.03411431


###################################################
#Early and late at 100% efficacy : Scenarios 5 and 6
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Find the common range of values
common_range.67 <- range(den6$x, den7$x)

# Define the number of breaks (you can adjust this)
num_breaks <- 100

# Create common breaks for both densities
z.points <- seq(common_range.67[1], common_range.67[2], length.out = num_breaks)

# Calculate density estimates with common breaks
den6_common <- density(InfectE100.540, from = common_range.67[1], to = common_range.67[2], n = num_breaks) #94%
den7_common <- density(InfectL100.540, from = common_range.67[1], to = common_range.67[2], n = num_breaks) #100%


#Step 4: Calculate Hellinger Distance
hellinger_distance <- (sum((sqrt(den6_common$y) - sqrt(den7_common$y))^2)*(z.points[2]-z.points[1]))^1/2  #use this one

cat("Hellinger Distance:", hellinger_distance, "\n")
#Hellinger Distance: 0.5098861 


#####################################################################################
########## Scenario plots for the Cumulative Reinfected at day 540. See Figure 5(b)


pdf("Reinfections.all.pdf")
plot(den11,  xlim= c(0, 50), ylim=c(0e+00, 0.17) , col = "blue", lty = 2, lwd = 2,
     xlab = substitute(paste(bold("Cumulative Reinfected"))), ylab = substitute(paste(bold("Density"))),
     main = " Scenario Comparison of Reinfected - Model 1 ")
lines(den22, col = "black" , lty = 2, lwd = 2)
lines(den44, col = "red3" , lty = 2, lwd = 2)
lines(den55, col = "violet", lty = 16, lwd = 2)
lines(den66, col = "darkturquoise" , lty = 4, lwd = 2)
lines(den77, col = "tan2", lty = 16, lwd = 2)
#lines(den3, col = "yellowgreen" , lty = 2, lwd = 2)

linetypes = c(2,1,4, 16, 4, 16, 32)
#names(linetypes) = c("94% Efficacy", "100% Efficacy", "Early vaccine-94%", "Late vaccine-94%", "Early vaccine-100%", "Late vaccine-100%", lwd = 3)
names(linetypes) = c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5", "Scenario 6", lwd = 2)

# legend("topright", c("94% Efficacy", "100% Efficacy", "Early vaccine-94%", "Late vaccine-94%", "Early vaccine-100%", "Late vaccine-100%"),
#        cex=0.75, col =c("blue","black","red3","violet", "darkturquoise","tan2"), lty = linetypes, lwd = 3)
legend("topright", c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5", "Scenario 6"),
       cex=0.85, col =c("blue","black","red3","violet", "darkturquoise","tan2"), lty = linetypes, lwd = 2)

dev.off()

######### Hellinger Distance calculation

###################################################
## 94 and 100 efficacy : Scenarios 1 and 2
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Find the common range of values
common_range.12 <- range(den11$x, den22$x)

# Define the number of breaks (you can adjust this)
num_breaks <- 10

# Create common breaks for both densities
z.points <- seq(common_range.12[1], common_range.12[2], length.out = num_breaks)

# Calculate density estimates with common breaks
den11_common <- density(ReinfectI.540, from = common_range.12[1], to = common_range.12[2], n = num_breaks) #94%
den22_common <- density(Reinfectk1.540, from = common_range.12[1], to = common_range.12[2], n = num_breaks) #100%

#Step 4: Calculate Hellinger Distance
hellinger_distance <- (sum((sqrt(den11_common$y) - sqrt(den22_common$y))^2)*(z.points[2]-z.points[1]))^1/2  #use this one
cat("Hellinger Distance:", hellinger_distance, "\n")
#Hellinger Distance: 0.0001795918 


#################################
#Scenario 1 and scenario 3: 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Find the common range of values
common_range.14 <- range(den11$x, den44$x)

# Define the number of breaks (you can adjust this)
num_breaks <- 10
hist(den11_common$y, breaks=10)

# Create common breaks for both densities
z.points <- seq(common_range.14[1], common_range.14[2], length.out = num_breaks)

#Step 3: Calculate Densities with Common Breaks
#Now, you need to re-calculate the densities using the common breaks:

# Calculate density estimates with common breaks
den11_common <- density(ReinfectI.540, from = common_range.14[1], to = common_range.14[2], n = num_breaks) #94%
den44_common <- density(ReinfectE.540, from = common_range.14[1], to = common_range.14[2], n = num_breaks) #100%

#Step 4: Calculate Hellinger Distance
hellinger_distance <- (sum((sqrt(den11_common$y) - sqrt(den44_common$y))^2)*(z.points[2]-z.points[1]))^1/2  #use this one

cat("Hellinger Distance:", hellinger_distance, "\n")

#Hellinger Distance: Hellinger Distance: 0.3755169 


###########################################
#Scenario 3 and 5: Early (94 and 100)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Find the common range of values
common_range.46 <- range(den44$x, den66$x)

# Define the number of breaks (you can adjust this)
num_breaks <- 10

# Create common breaks for both densities
z.points <- seq(common_range.46[1], common_range.46[2], length.out = num_breaks)

# Calculate density estimates with common breaks
den44_common <- density(ReinfectE.540, from = common_range.46[1], to = common_range.46[2], n = num_breaks) #94%
den66_common <- density(ReinfectE100.540, from = common_range.46[1], to = common_range.46[2], n = num_breaks) #100%

#Step 4: Calculate Hellinger Distance
hellinger_distance <- (sum((sqrt(den44_common$y) - sqrt(den66_common$y))^2)*(z.points[2]-z.points[1]))^1/2  #use this one

cat("Hellinger Distance:", hellinger_distance, "\n")
##Hellinger Distance: 1.235194e-05 



###################################################
#Late (94 and 100 efficacy): Scenario 4 and Scenario 6
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Find the common range of values
common_range.57 <- range(den55$x, den77$x)

# Define the number of breaks (you can adjust this)
num_breaks <- 10

# Create common breaks for both densities
z.points <- seq(common_range.57[1], common_range.57[2], length.out = num_breaks)

# Calculate density estimates with common breaks
den55_common <- density(ReinfectL.540, from = common_range.57[1], to = common_range.57[2], n = num_breaks) #94%
den77_common <- density(ReinfectL100.540, from = common_range.57[1], to = common_range.57[2], n = num_breaks) #100%


#Step 4: Calculate Hellinger Distance
hellinger_distance <- (sum((sqrt(den55_common$y) - sqrt(den77_common$y))^2)*(z.points[2]-z.points[1]))^1/2  #use this one

cat("Hellinger Distance:", hellinger_distance, "\n")
#### Hellinger Distance: 0.03300149   



###################################################
#Late and early (100% efficacy): Scenarios 5 and 6
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Find the common range of values
common_range.67 <- range(den66$x, den77$x)

# Define the number of breaks (you can adjust this)
num_breaks <- 10

# Create common breaks for both densities
z.points <- seq(common_range.67[1], common_range.67[2], length.out = num_breaks)

# Calculate density estimates with common breaks
den66_common <- density(ReinfectE100.540, from = common_range.67[1], to = common_range.67[2], n = num_breaks) #94%
den77_common <- density(ReinfectL100.540, from = common_range.67[1], to = common_range.67[2], n = num_breaks) #100%

#Step 4: Calculate Hellinger Distance
hellinger_distance <- (sum((sqrt(den66_common$y) - sqrt(den77_common$y))^2)*(z.points[2]-z.points[1]))^1/2  #use this one

cat("Hellinger Distance:", hellinger_distance, "\n")

#Hellinger Distance: 0.5310828 


#####################################################################################
############%%%%%%%%%%  Scenario plots for the Cumulative Deaths at day 540. See Figure 5(c)


pdf("Deaths.all.pdf")
plot(den1.1,  xlim= c(480, 780), ylim=c(0, 0.0165), col = "blue", lty = 2, lwd = 2, xlab = substitute(paste(bold("Cumulative Deaths"))),
     ylab = substitute(paste(bold("Density"))), main = "Scenario Comparison of Deaths - Model 1")
lines(den2.2, col = "black" , lty = 2, lwd = 2)
lines(den4.4, ylab = c(0,0.02), col = "red3" , lty = 2, lwd = 2)
lines(den5.5, col = "violet", lty = 16, lwd = 2)
lines(den6.6, col = "darkturquoise" , lty = 4, lwd = 2)
lines(den7.7, col = "tan2", lty = 16, lwd = 2)
linetypes = c(2,1,4, 16, 4, 16, 32, 2)
names(linetypes) = c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5", "Scenario 6", lwd = 3)
legend("topright", c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5", "Scenario 6"),
       cex=0.85,  col =c("blue","black","red3","violet", "darkturquoise","tan2"), lty = linetypes, lwd = 3)
dev.off()

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############# Hellinger Distance calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Select Equally Spaced Points: 
#Evaluate PDFs:For each equally spaced point, calculate the values of the probability density functions (PDFs) of the two distributions being compared.
#Calculate Hellinger Distance: 
?density


common_range <- range(den1.1$x, den2.2$x) 

# Define the number of breaks (you can adjust this)
num_breaks <- 100

# # Choosing equally spaced points along the range
z.points <- seq(common_range[1], common_range[2], length.out = num_breaks)

#Step 3: Calculate Densities with Common Breaks
#Now, you need to re-calculate the densities using the common breaks:

# Calculate density estimates with common breaks
den1.1_common <- density(DeathsI.540, from = common_range[1], to = common_range[2], n = num_breaks) #94%
den2.2_common <- density(Deathsk1.540, from = common_range[1], to = common_range[2], n = num_breaks) #100%

hist(DeathsI.540, breaks=100)

####This part is important for error correction, see below for the reason
h1 <- den1.1_common$bw  #bandwidth
h2 <- den2.2_common$bw  #bandwidth

st1 <- h1*sum(den1.1_common$y)
st2 <- h2*sum(den2.2_common$y)
m1 <- den1.1_common$y/st1  # .......... (*)
m2 <- den2.2_common$y/st2   # .......... (**)

#Step 4: Calculate Hellinger Distance: With the density estimates using common breaks and equally spaced points, you can calculate the Hellinger distance:

# Calculate the Hellinger distance
hellinger_distance1 <- (sum((sqrt(den1.1_common$y) - sqrt(den2.2_common$y))^2)*(z.points[2]-z.points[1]))^1/2 # will use this

#Important: For whatever reason, the density function in R creates densities that sum to 1+1/(2N), give or take some floating point error
#in the calculation. Although this is a small error, we can compensate by dividing the values returned by densites by the values in (*) and (**)

hellinger_distance2 <- (sum((sqrt(m1) - sqrt(m2))^2)*(z.points[2]-z.points[1]))^1/2 # will use this

# Print the Hellinger distance
cat("Hellinger Distance:", hellinger_distance1, "\n")

cat("Hellinger Distance:", hellinger_distance2, "\n")

# Hellinger Distance: 9.606272e-05 
# > cat("Hellinger Distance:", hellinger_distance2, "\n")
# Hellinger Distance: 9.007374e-05 
#notice that the distance is slightly smaller but negligible, so we will use the Hellinger distance1.

#This idea is used to calculate the remaining hellinger distance between other distributions.

############################################################################
#Early vaccination (94 ad 100 %)
# Find the common range of values
common_range.46 <- range(den4.4$x, den6.6$x)

# Define the number of breaks (you can adjust this)
num_breaks <- 100

# Choosing equally spaced points along the range
z.points <- seq(common_range.46[1], common_range.46[2], length.out = num_breaks)

# Calculate density estimates with common breaks
den4.4_common <- density(DeathsE.540, from = common_range.46[1], to = common_range.46[2], n = num_breaks) #94%
den6.6_common <- density(DeathsE100.540, from = common_range.46[1], to = common_range.46[2], n = num_breaks) #100%

# 
# h1 <- den4.4_common$bw  #bandwidth
# h2 <- den6.6_common$bw  #bandwidth
# 
# st1 <- h1*sum(den4.4_common$y)
# st2 <- h2*sum(den6.6_common$y)
# m1 <- den4.4_common$y/st1  
# m2 <- den6.6_common$y/st2  
# 
#hellinger_distance <- (sum((sqrt(m1) - sqrt(m2))^2)*(z.points[2]-z.points[1]))^1/2 # will use this
hellinger_distance <- (sum((sqrt(den4.4_common$y) - sqrt(den6.6_common$y))^2)*(z.points[2]-z.points[1]))^1/2  #use this one

cat("Hellinger Distance:", hellinger_distance, "\n")
####Hellinger Distance: 0.0002776956  


################################################################################
#Late vaccination (94 and 100 %)
# Find the common range of values
common_range.57 <- range(den5.5$x, den7.7$x)

# Define the number of breaks (you can adjust this)
num_breaks <- 100

## Choosing equally spaced points along the range
z.points <- seq(common_range.57[1], common_range.57[2], length.out = num_breaks)

# Calculate density estimates with common breaks
den5.5_common <- density(DeathsL.540, from = common_range.57[1], to = common_range.57[2], n = num_breaks) #94%
den7.7_common <- density(DeathsL100.540, from = common_range.57[1], to = common_range.57[2], n = num_breaks) #100%


hellinger_distance <- (sum((sqrt(den5.5_common$y) - sqrt(den7.7_common$y))^2)*(z.points[2]-z.points[1]))^1/2  #use this one

cat("Hellinger Distance:", hellinger_distance, "\n")
#### Hellinger Distance: 0.03115741   

###################################################
#Early and late 100% efficacy
###################################################
# Find the common range of values
common_range.67 <- range(den6.6$x, den7.7$x)

# Define the number of breaks (you can adjust this)
num_breaks <- 100

# Choosing equally spaced points along the range
z.points <- seq(common_range.67[1], common_range.67[2], length.out = num_breaks)

# Calculate density estimates with common breaks
den6.6_common <- density(DeathsE100.540, from = common_range.67[1], to = common_range.67[2], n = num_breaks) #94%
den7.7_common <- density(DeathsL100.540, from = common_range.67[1], to = common_range.67[2], n = num_breaks) #100%


hellinger_distance <- (sum((sqrt(den6.6_common$y) - sqrt(den7.7_common$y))^2)*(z.points[2]-z.points[1]))^1/2  #use this one
cat("Hellinger Distance:", hellinger_distance, "\n")
#### Hellinger Distance: 0.6443222   

######################## THE END #####################################

# Note: If you are interested in reproducing Figure 6 for research purposes, please contact me for the code. The code is basically the same as I have above, the only difference is the "R.dat files" you might need if you don't perform the scenario analysis by yourself using the "SVEIRDI2R2.Model2.R" code.
