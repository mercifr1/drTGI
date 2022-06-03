# drTGI
Impact of dosing history on the assessment of TGI

## Description of the files from *Maia/Rscript*

1. **SimPred.R**: contains the construction of models under the analytical form.
 One can found the models used for the simulations for a flat/constant dose and time-varying dose and also the corresponding design matrices. To check the implementation of the model, the simulations where tested for the design matrices with DOSE=10 for the flat dose level, and  DOSE=(10,5) for time varying dose.
 
 **Note**: In mrgsolve, using $PRED or $TABLE to exppress the form of the model, by construction these 2 blocks does not capture the time-varying covariates. Thus, the model structure for the time-varying was implemented recursively.

2. **SimODE.R** : contains the construction models under the ODE form. In this file, one can found the model structure expressed in the mrgsolve syntax, the design matrices for the flat dose and time-varying dose and the simulations for both types of doses. To check the implementation of the model and the design matrices, the simulations were tested for the DOSE=10 flat dose level and for DOSE=(10,5) for time varying dose.



Once the syntax of the model and the results where satisfactory, the next step is to build different scenarios step-by-step for a better approach to the real trial events.

3. **Scenario_1.R**: balanced design matrix: 
    - dose-cohorts with same number of individual/cohort:
       - 5 cohorts
       - 5 individual/cohort
       - 25 individual in total

    -  each cohort have a different flat/constant dose=> 5 dose levels 
    -  all individuals:  same number of observations 
 
 4. **Scenario_2.R**: unbalanced design matrix: **NOTE**: to run this file it is necessary to run and source the **SimODE.R** before
     - dose-cohorts with different number of individual/cohort:
       - 5 cohorts
       - number of individual per cohort follows a truncated Poisson distr. (lambda_N=4)
       - obtain between 3 and 10 individuals per cohort.

    -  each cohort have a different flat/constant dose=> 5 dose levels 
    -  number of observations per individual is different and it follows a zero-truncated Poisson distr. (lambda_0=4)
    -  lambda_0 reflect the short follow-up or early treatment discontinuation


 5. **Scenario_2A.R**: unbalanced design matrix: Same scenario as **Scenario_2.R**
 
 **NOTE**:  this file contain a proper R code for the future design matrix compared to the **Scenario_2.R**, where I forced some value by hand to take a certain value.
           to run this file it is necessary to run and source the **SimODE.R** before.
           
           
   **NOTE**: Next updated need to be implemented: 
   
             - dose reductions for  3rd ,4th and 5 th cohorts ( toxicity)
             - remove the measurements that incresed more then 20% from nadir.

