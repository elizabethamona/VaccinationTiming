# Studying Disease Reinfection Rates, Vaccine Efficacy and the Timing of Vaccine Rollout in the context of Infectious Diseases
### Authors: Elizabeth Amona, Indranil Sahoo, Ed Boone, Ryad Ghanam

*What is the goal of this work?* Maybe add the abstract here... 

**Model 1:** To obtain the parameter estimates and Psuedo-R^2 like I have in Table 1 (in the paper) and the model fits Figure 3(a-f), you need three major files: "S1EI1R1DVS2I2R2.Model.R", "MH55.R" and "MH55.cpp". In addition, you need another file called the "StartStepValuesMC8_new.RData". You can generate your own "Q1MCMC20.Rdat" file after you have run the code once.

**Instruction:**

   1. Ensure that you can run a CPP file by downloading necessary files
   2. Source "MH55.R" and "MH55.cpp" using the source function in R (see the code )
   3. Read in the data and run the initial/start values
   4. Load the "StartStepValuesMC8_new.RData" 
   5. Plot the states and obtain your estimates



**Model 2:** To obtain the parameter estimates and Psuedo-R^2 in Table 1 (in the paper) and the model fits Figure 4(a-f), you need three major files: "SVEIRDI2R2.Model2.R", "SVEIRD.R" and "SVEIRD.cpp". In addition, you need two other files: "StartStepValuesMC8_new.RData" and "Q1MCMC20_SV.Rdat". 

**Instruction:** 

1.  Ensure that you can run a CPP file by downloading necessary files
2.  Source "SVEIRD.R" and "SVEIRD.cpp" using the source function in R (see the code )
3.  Read in the data and run the initial/start values
4.  Load the "StartStepValuesMC8_new.RData"
5.  Plot the states and obtain your estimates


**Scenario Analysis and Hellinger Distance:** To reproduce a similar plot in Figure 5, run the "Scenario_HellingerDist_Analysis.R" file. This code will plot the density plots for all the six scenarios and their Hellinger distance calculation.

**Instruction:** Ensure to load all the .Rdat files you obtain when you run your scenario analysis. I already generated my results, so I just "load" the results for all six scenarios in the "Scenario_HellingerDist_Analysis.R" code provided. Note that those are the results obtained from the scenario analysis I performed. If you want to perform your own scenario analysis, please use my main code and consider all the six scenarios I included in the "Scenario_HellingerDist_Analysis.R" file.

**In summary:**
- There are a total of **"7 R files" and "1 .Rdat file"**. Carefully run this file in the directory you have downloaded them all.
- I would suggest you create two different folders for the two models so you can run the models seperately to avoid any form of mistake when running the code.
- Also, when performing the scenario analysis, ensure that the .Rdat file is in the same folder. 


