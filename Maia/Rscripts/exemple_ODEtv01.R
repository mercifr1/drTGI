
# Simple design matrix
data<-tibble(
  ID=1,
  TIME=c(0,24,31,46,48,72,82),
  EVID=c(1,1,0,0,1,1,0),
  CMT=c(1,1,2,2,1,1,2),
  AMT=c(80,60,0,0,40,60,0),
  MDV=1,
  COV=c(80,60,60,60,40,60,60)
  
)

# # A tibble: 7 x 7
#      ID  TIME  EVID   CMT   AMT   MDV   COV
#     <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1     1     0     1     1    80     1    80
# 2     1    24     1     1    60     1    60
# 3     1    31     0     2     0     1    60
# 4     1    46     0     2     0     1    60
# 5     1    48     1     1    40     1    40
# 6     1    72     1     1    60     1    60
# 7     1    82     0     2     0     1    60
# 


#design matrix: replace EVID 0 with EVID 2 
data0<-tibble(
  ID=1,
  TIME=c(0,24,31,46,48,72,82),
  EVID=c(1,1,2,2,1,1,2),
  CMT=c(1,1,2,2,1,1,2),
  AMT=c(80,60,0,0,40,60,0),
  MDV=1,
  COV=c(80,60,60,60,40,60,60)
  
)


#  close form solution of the ODE

code00<-'
$PARAM TVKA = 0.05, CL = 1, VC = 10, COV = 100
$CMT DEPOT CENTRAL
$TABLE
double DV  = CENTRAL / VC  ;
$MAIN
double KA = TVKA * (COV / 100) ;
double KE = CL / VC ;
double V=VC;
$PKMODEL depot = TRUE
$CAPTURE DV, KA,COV
'

# ODE
code01<-'
$PARAM TVKA = 0.05, CL = 1, VC = 10, COV = 100
$CMT DEPOT CENTRAL
$TABLE
double DV  = CENTRAL / VC  ;
$MAIN
double KA = TVKA * (COV / 100) ;
double KE = CL / VC ;
$ODE
dxdt_DEPOT   = - KA * DEPOT ;
dxdt_CENTRAL =   KA * DEPOT  - KE * CENTRAL ;
$CAPTURE DV, KA,COV
'

# Simulation for close form sol
mode00<-mcode("foo",code00)
mrgsim(mode00,data)%>%as.data.frame
#   ID TIME     DEPOT  CENTRAL       DV   KA COV
# 1  1    0  80.00000  0.00000 0.000000 0.04  80
# 2  1   24  98.94018 13.57832 1.357832 0.03  60
# 3  1   31  80.19935 20.05727 2.005727 0.03  60
# 4  1   46  51.13736 18.72215 1.872215 0.03  60
# 5  1   48  89.13224 17.14453 1.714453 0.02  40
# 6  1   72 103.38532 16.68364 1.668364 0.03  60
# 7  1   82  76.58973 22.66174 2.266174 0.03  60


mrgsim(mode00,data0,end=-1)%>%as.data.frame
#    ID TIME     DEPOT  CENTRAL       DV   KA COV
# 1  1    0  80.00000  0.00000 0.000000 0.04  80
# 2  1   24  98.94018 13.57832 1.357832 0.03  60
# 3  1   31  80.19935 20.05727 2.005727 0.03  60
# 4  1   46  51.13736 18.72215 1.872215 0.03  60
# 5  1   48  89.13224 17.14453 1.714453 0.02  40
# 6  1   72 103.38532 16.68364 1.668364 0.03  60
# 7  1   82  76.58973 22.66174 2.266174 0.03  60



#simulation for ODE
mode01<-mcode("foo1",code01)
mrgsim(mode01,data)%>%as.data.frame


#   ID TIME     DEPOT  CENTRAL       DV   KA COV
# 1  1    0  80.00000  0.00000 0.000000 0.04  80
# 2  1   24  98.94018 13.57832 1.357832 0.03  60
# 3  1   31  80.19935 20.05727 2.005727 0.03  60
# 4  1   46  51.13736 18.72216 1.872216 0.03  60
# 5  1   48  89.13224 17.14453 1.714453 0.02  40
# 6  1   72 103.38532 16.68364 1.668364 0.03  60
# 7  1   82  76.58973 22.66174 2.266174 0.03  60



# Simulationm ODE: data matrix with EVID 2
mrgsim(mode01,data0,end=-1)%>%as.data.frame
#    ID TIME     DEPOT  CENTRAL       DV   KA COV
# 1  1    0  80.00000  0.00000 0.000000 0.04  80
# 2  1   24  98.94018 13.57832 1.357832 0.03  60
# 3  1   31  80.19935 20.05727 2.005727 0.03  60
# 4  1   46  51.13736 18.72215 1.872215 0.03  60
# 5  1   48  89.13224 17.14453 1.714453 0.02  40
# 6  1   72 103.38532 16.68364 1.668364 0.03  60
# 7  1   82  76.58973 22.66174 2.266174 0.03  60