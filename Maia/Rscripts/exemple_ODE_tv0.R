library(mrgsolve)
library(dplyr)
library(ggplot2)
library(purrr)



# test with ODE

#create the design matrix
adms <- tibble(ID = 1, time = c(0,24,48), evid = 1, cmt = 1, amt = 100, mdv = 1, COV = c(100, 80, 60), DV = NA_real_)
obss <- tibble(ID = 1, time = 72, evid = 0, cmt = 1, mdv = 0, COV = 60)

data_lag <- bind_rows(adms, obss) %>% 
  bind_rows(mutate(bind_rows(adms, obss), ID = 2, COV = 100))
#A tibble: 8 x 8
#ID  time  evid   cmt   amt   mdv   COV    DV
#<dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#  1     0     1     1   100     1    100   NA
#  1    24     1     1   100     1    80    NA
#  1    48     1     1   100     1    60    NA
#  1    72     0     1    NA     0    60    NA
#  2     0     1     1   100     1   100    NA
#  2    24     1     1   100     1   100    NA
#  2    48     1     1   100     1   100    NA
#  2    72     0     1    NA     0   100    NA

code <- '
$PARAM CL = .5, VC = 50, LAG = 1, COV = 100

$CMT CENTRAL

$TABLE double DV = CENTRAL / VC ;

$MAIN 
double KE = CL/VC ;
double F1 = COV/100 ;

ALAG_CENTRAL = LAG ;

F_CENTRAL = F1 ;

$ODE
dxdt_CENTRAL = - KE * CENTRAL ;
$CAPTURE DV, F1, KE
'
model_lag <- mcode("model_lag", code)
mrg_lag <- model_lag %>% 
  data_set(data_lag) %>% 
  mrgsim_df()
#    ID time   CENTRAL       DV  F1   KE
# 1  1    0   0.00000 0.000000 1.0 0.01
# 2  1   24  63.56269 1.271254 0.8 0.01
# 3  1   48  97.67220 1.953444 0.6 0.01
# 4  1   72 124.50369 2.490074 0.6 0.01
# 5  2    0   0.00000 0.000000 1.0 0.01
# 6  2   24  79.45336 1.589067 1.0 0.01
# 7  2   48 141.95359 2.839072 1.0 0.01
# 8  2   72 191.11801 3.822360 1.0 0.01


# Here I tried to put the close form solution
# The code is not correct, need to be changed
code2 <- '
$PARAM CL = .5, VC = 50, LAG = 1, COV = 100



$TABLE 
double CENTRAL_0=0;
double DV = CENTRAL / VC ;
double CENTRAL=CENTRAL_0*exp(-KE*TIME);


$MAIN 
double KE = CL/VC ;
double F1 = COV/100 ;

$CAPTURE DV CENTRAL F1 KE
'
model_lag_CF <- mcode("model_lag_CF", code2)

data_lag_CF<-data_lag%>%mutate(cmt=0)

mrg_lagCF <- model_lag_CF %>% 
  data_set(data_lag_CF) %>% 
  mrgsim_df()
#    ID time DV CENTRAL  F1   KE
# 1  1    0  0       0 1.0 0.01
# 2  1   24  0       0 0.8 0.01
# 3  1   48  0       0 0.6 0.01
# 4  1   72  0       0 0.6 0.01
# 5  2    0  0       0 1.0 0.01
# 6  2   24  0       0 1.0 0.01
# 7  2   48  0       0 1.0 0.01
# 8  2   72  0       0 1.0 0.01