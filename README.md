# drTGI
Impact of dosing history on the assessment of TGI
## Files and folders from *Maia/Rscript*
  Simulation files for TGI Claret model:
  - **SimODE.R** : ODE model syntax mrgsolve 
  - **SimPred.R**: Analytical form of the ODE syntax mrgsolve
  - **functions.R**: Rscript with useful functions to built the design matrices for simulation
  - **Simulation_SLD_complete.R**: simulation of SLD time profiles for different scenarios of design matrices
  
  Fit model in Stan for TGI Claret model:
  - **TGI_hierarchicalODE_random_effects.stan**: TGI hierarchical ODE model code in stan 
  - **TGI_hierarchicalODE_random_effects.R** : R script with the model fit with cmdstan 
  
 Extra Rscript:
 - **CheckNonCenterParametrizationStanODE.R**: intuitive check for non-centered parametrization of the parameters in the TGI ODE model 
 - **plot_parameters_model_TGI.R** : plot density parameters, random effects from Simulated dataset
 - **BRMSAnalyticFormFit.R**: fit TGI ODE model on a simple flat dose-cohort with **brms** R package
 - **TGI_SLD_ODE_nlmixr_fit.R**: fit TGI ODE model on a simulated dataset with flat dose cohort with **nlmixr** R package
 - **S2_Dose_reduction_function.R** : function for Scenario 2 Dose reduction design matrix construction for simulations
 - **S3_Dose_omission_function.R** : function for Scenario 3 Dose omission design matrix construction for simulations
 - **Scenario_1_function.R** :Scebario 1 design matrix flat dose-cohort unbalanced matrices for simulations
 - **Simulation_SLD_Scenario_0.R** : Simulated SLD times profile on a balanced design matrice 
 - **Simulation_SLD_Scenario_1.R** : Simulated SLD times profile on an unbalanced design matrice 

Extra R folder:
-**Other version ODE TGI stan codes**: a second version of solving the TGI ODE in stan 
-**Example_hierarchicalODE_theoph**: example hierarchical ODE on theophylline data set : Stan codes
-**First test ode sld integartor Rstan**: check how the output of the stan code for ODE and integartor taht solve the ODE works in stan  
-**TGI Simple version stan ODE** : stan code for simple version TGI ODE mode: for 1 subject, 3 subjects, constant standard error, centered parametrisation 

## Description of the files from *Maia/Rscript*

1. **SimPred.R**: contains the construction of models under the analytical form.
 One can found the models used for the simulations for a flat/constant dose and time-varying dose and also the corresponding design matrices to test thid models. To check the implementation of the model, the simulations where tested for the design matrices with DOSE=10 for the flat dose level, and  DOSE=(10,5) for time varying dose.
 
   **Note**: In mrgsolve, using $PRED or $TABLE to exppress the form of the model, by construction these 2 blocks does not capture the time-varying covariates. Thus, the model structure for the time-varying was implemented recursively.

2. **SimODE.R** : contains the construction of models under the ODE form. In this file, one can found the model structure expressed in the mrgsolve syntax, the design matrices for the flat dose and time-varying dose and the simulations for both types of dose-events. To check the implementation of the model and the design matrices, the simulations were tested for the DOSE=10 flat dose level and for DOSE=(10,5) for time varying dose.

3. **BRMSAnalyticFormFit.R**: fit the brms on the analytical form **SimPred.R**.
 This file represent a first test to understand the synatx, implementation and fit of the **brms** as well as the checking of the outcome of the model. It was use a flat dose (DOSE=10) and a time varying dose (DOSE=5,10), the same as in the SimPred.R.
 
 
 
Once the syntax of the model and the results where satisfactory, the next step is to build different scenarios step-by-step for a better approach to the real trial events.

Built the function to create the scenarios for the dose- manipulation.

 - **Scenarion_1_function.R**:unbalanced design matrix: create a function that generate the number of participants per cohort and the number of observations per  individual folowing the truncated Poisson distribution (with parapeter lambda_N=4) and zero-truncated Poisson distribution ( with parameter lambda_0=4) respectively.


 - **S2_Dose_reduction_function** : Scenario 2 of dose reduction :  using the dataset generated from the Scenario_1_function, built a function for the dose reduction in the specific cohorts (i.e: cohorts 3, 4 and 5).


 - **S3_Dose_omission_function** : Scenario 3 of dose omission : using the same dataset generated from  the scenario_1_function, built a function for the dose omission independet from the dose reduction in the specific cohorts (i.e: cohorts 3, 4 and 5)

 - **functions**: this files contains all the above functions: useful to generate the entire simulation.

 



4. **Simulation_SLD_Scenario_0.R**: balanced design matrix: 
    - dose-cohorts with same number of observations/individual:
       - 5 dose-cohorts : dose levels equal to 1/6, 1/2, 1, 2 and 3 mg
       - generate the number of  individual/cohort ~ **TruncatedPoisson(lambda_N=4)**
       - total number of individul per trial is equal to the sum of the individual/cohort
       - each cohort have a different flat/constant dose => the corresponding 5 dose levels 
       - **all individuals:  equal number of observations** 
       
    -  simulate individual SLD time profiles with mrgsim()
    -  repeate each trail multiple time and store the resulting simulation dataset as a data frame
    -  the simulated individual SLD time profiles are pruned by discarding data points past the first instance of 20% increase from nadir
    -  plot the SLD profile for one of the simulation
    
    **Remark**: source the SimODE.R
    
 
 5. **Simulation_SLD_Scenario_1.R**: unbalanced design matrix:
     - dose-cohorts with different number of participants/cohort and different number of observations/individual:
       - 5 dose-cohorts : dose levels equal to 1/6, 1/2, 1, 2 and 3 mg
       - generate the number of  individual/cohort ~ **TruncatedPoisson(lambda_N=4)**
       - total number of individul per trial is equal to the sum of the individual/cohort
       - generate the number of  observations/individual ~ **ZeroTruncatedPoisson(lambda_0=4)**
       - each cohort have a different flat/constant dose => the corresponding 5 dose levels 
       - **NOTE**: lambda_0 reflect the short follow-up or early treatment discontinuation
       
    -  simulate individual SLD time profiles with mrgsim()
    -  repeate each trail multiple time and store the resulting simulation dataset as a data frame
    -  the simulated individual SLD time profiles are pruned by discarding data points past the first instance of 20% increase from nadir
    -  plot the SLD profile for one of the simulation

   **Remark**:  source the **SimODE.R** and **Scenario_1_function** before compiling this file

 6. **Simulation_SLD_complete.R**: unbalanced design matrix with dose modifications.
 
     **NOTE** This file contains the replicated design matrix and the simulations, implemented with all scenarios and constratints to reflect the actual clinical trial conditions.
 
    A. Built the design matrix
       - Step1: generate the initial unbalanced design matrix: use  **Scenario_1_function**
       - Step2: **Scenario 2**:  Dose reduction : use **dose_reduction_function** on the unbalanced design matrix generated on Step1 to built the dose reduction event on the cohorts of interest. Store the dataset with the modifications.   
       - Step3: **Scenario 3**: Dose omission: use **dose_omission_function** on the unbalanced design matrix generated on Step1 to built the dose omission event on the cohorts  of interest. Store the dataset with the modifications. **Note**: dose_omission is computed independent from the dose_reduction
       - Step4:**Scenario 4**: Combine the dose_reduction dataset and dose_omission dataset from step2 and step3
       - Step5: Merge the dataset from Step4 with the initial unbalanced design matrix from Step1.( from Step1 matrix  selecte only the cohorts which are  apart from the dose modification: i.e cohort 1 and 2 respectively
       
   **Remark 1**: Each scenario implemetation contains multiple flag columns, in order to merge and bind correctly the final design matrix. 
   
   **Remark 2**: We generate(repeat) multiple times the design matrix then on each repeatition simulate the SLD profiles.
              
   
    B. Simulations 
      - simulate individual SLD time profiles with mrgsim() using the design matrix from A for all scenarios independet
      - repeate each trail multiple time and store the resulting simulation dataset as a data frame
      - the simulated individual SLD time profiles are pruned by discarding data points past the first instance of $20%$ increase from nadir
      - plot the SLD profile for one of the simulation
   
   **Remark**: source the **SimODE.R** and **functions.R** before compiling this file.
           
           
  
