**Model 1:** To obtain the parameter estimates and Psuedo-R^2 in Table 1 (in the paper) and the model fits Figure 3(a-f), you need three major files: "S1EI1R1DVS2I2R2.Model1.R", "MH55.R" and "MH55.cpp". In addition, you need two other files: "StartStepValuesMC8_new.RData" and "Q1MCMC20.Rdat". 

**Instruction:**

   1. Ensure that you can run a CPP file by downloading necessary files
   2. Source "MH55.R" and "MH55.cpp" using the source function in R (see the code )
   3. Read in the data and run the initial/start values
   4. Load the "StartStepValuesMC8_new.RData" 
   5. Load "Q1MCMC20.Rdat" to obtain the 50,000 samples I generated. 
   6. Plot the states and obtain the estimates



**Model 2:** To obtain the parameter estimates and Psuedo-R^2 in Table 1 (in the paper) and the model fits Figure 4(a-f), you need three major files: "SVEIRDI2R2.Model2.R", "SVEIRD.R" and "SVEIRD.cpp". In addition, you need two other files: "StartStepValuesMC8_new.RData" and "Q1MCMC20_SV.Rdat". 

**Instruction:** 

1.  Ensure that you can run a CPP file by downloading necessary files
2.  Source "SVEIRD.R" and "SVEIRD.cpp" using the source function in R (see the code )
3.  Read in the data and run the initial/start values
4.  Load the "StartStepValuesMC8_new.RData" 
5. Load "Q1MCMC20_SV.Rdat" to obtain the 50,000 samples I generated. 
6.  Plot the states and obtain the estimates


**Scenario Analysis and Hellinger Distance (for models 1 and 2):** To reproduce the plots in Figure 5, run the "Scenario_HellingerDist_Analysis.R" file. This code will reproduce the density plots for all the six scenarios and their Hellinger distance calculation.

**Instruction:** Ensure to load all the .Rdat files in the "Scenario_HellingerDist_Analysis.R" code. Those are the results obtained from the scenario analysis I performed. You don't need all the scenario analysis code since it is similar to the main code. All you need is the results I have generated. If you want to perform your own scenario analysis, please use my main code and consider all the six scenarios I included in the "Scenario_HellingerDist_Analysis.R" file.

**In summary:**
- There are a total of **"9 R file" and "16 .Rdat files"**. Carefully run this file in the directory you have downloaded them all.
- I would suggest you create two different folders for the two models so you can run the models seperately to avoid any form of mistake when running the code.
- Also, when performing the scenario analysis, ensure that all the .Rdat files needed are in the same folder. 


