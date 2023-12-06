# Studying Disease Reinfection Rates, Vaccine Efficacy and the Timing of Vaccine Rollout in the context of Infectious Diseases
### Authors: Elizabeth Amona, Indranil Sahoo, Ed Boone, Ryad Ghanam

**Abstract:**  Qatar has undergone distinct waves of infections, compounded by the emergence of variants, posing additional complexities. This research uniquely delves into the varied efficacy of existing vaccines and the pivotal role of vaccination timing in the context of COVID-19. Departing from conventional modeling, we introduce two models that account for the impact of vaccines on infections, reinfections, and deaths. Recognizing the intricacy of these models, we use the Bayesian framework and specifically utilize the Metropolis-Hastings Sampler for estimation of model parameters. Scenario analyses on the two models are conducted, quantifying the similarities in probability distributions of infected, reinfected, and deaths using the Hellinger distance metric. Comparative analysis, employing the Bayes factor, underscores the plausibility of a model assuming a different susceptibility rate to reinfection, as opposed to assuming the same susceptibility rate for both infections and reinfections. Results highlight the adverse outcomes associated with delayed vaccination, emphasizing the efficacy of early vaccination in reducing infections, reinfections, and deaths. We advocate prioritizing early vaccination, irrespective of vaccine efficacy, as a key strategy in effectively combating future pandemics. This study contributes vital insights for evidence-based public health interventions, providing clarity on vaccination strategies and reinforcing preparedness for challenges posed by infectious diseases.

### The Models Description 
**Model 1:** To obtain the parameter estimates and Psuedo-R^2 like I have in Table 1 (in the paper) and the model fits Figure 3(a-f), you need three major files: "S1EI1R1DVS2I_IR_R.Model.R", "MH55.R" and "MH55.cpp". In addition, you need another file called the "StartStepValuesMC8_new.RData". You can generate your own "Q1MCMC20.Rdat" file after you have run the code once.

**Instruction for Model 1 Analysis:**

   1. Ensure that you can run a CPP file by downloading necessary files
   2. Source "MH55.R" and "MH55.cpp" using the source function in R (see the code )
   3. Read in the data and run the initial/start values
   4. Load the "StartStepValuesMC8_new.RData" 
   5. Plot the states and obtain your estimates



**Model 2:** To obtain the parameter estimates and Psuedo-R^2 in Table 1 (in the paper) and the model fits Figure 4(a-f), you need three major files: "SVEIRDI_IR_R.Model2.R", "SVEIRD.R" and "SVEIRD.cpp". In addition, you need two other files: "StartStepValuesMC8_new.RData" and "Q1MCMC20_SV.Rdat". 

**Instruction for Model 2 Analysis:**

1.  Ensure that you can run a CPP file by downloading necessary files
2.  Source "SVEIRD.R" and "SVEIRD.cpp" using the source function in R (see the code )
3.  Read in the data and run the initial/start values
4.  Load the "StartStepValuesMC8_new.RData"
5.  Plot the states and obtain your estimates


**Scenario Analysis and Hellinger Distance:** To reproduce a similar plot in Figure 5, run the "ScenarioHellingerDistAnalysis.R" file. This code will plot the density plots for all the six scenarios and their Hellinger distance calculation.

**Instruction:** Ensure to load all the .Rdat files you obtain when you run your scenario analysis. I already generated my results, so I just "load" the results for all six scenarios in the "ScenarioHellingerDistAnalysis.R" code provided. Note that those are the results obtained from the scenario analysis I performed. If you want to perform your own scenario analysis, please use my main code and consider all the six scenarios I included in the "ScenarioHellingerDistAnalysis.R" file.

**In summary:**
- There are a total of **"7 R files" and "1 .Rdat file"**. Carefully run this file in the directory you have downloaded them all.
- I would suggest you create two different folders for the two models so you can run the models seperately to avoid any form of mistake when running the code.
- Also, when performing the scenario analysis, ensure that the .Rdat file is in the same folder. 


