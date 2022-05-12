
#DATA#

datatest <- tibble(ID = 1, time = c(0, 24, 48), evid = 1, cmt = 1, amt = 1000, mdv = 1,  COV = 1) %>%
  add_row(ID = 1, time = c(49, 50), evid = 0, cmt = 2, amt = 0, mdv = 1, COV = 2) %>%
  bind_rows(mutate(., ID = 2, COV = 1), mutate(., ID = 3, COV = 2)) %>%
  mutate(DV = NA_real_)
print(datatest)
# A tibble: 15 x 8
#      ID  time  evid   cmt   amt   mdv   COV    DV
#     <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1      1     0     1     1  1000     1     1    NA
# 2     1    24     1     1  1000     1     1    NA
# 3     1    48     1     1  1000     1     1    NA
# 4     1    49     0     2     0     1     2    NA
# 5     1    50     0     2     0     1     2    NA
# 6     2     0     1     1  1000     1     1    NA
# 7     2    24     1     1  1000     1     1    NA
# 8     2    48     1     1  1000     1     1    NA
# 9     2    49     0     2     0     1     1    NA
# 10     2    50     0     2     0     1     1    NA
# 11     3     0     1     1  1000     1     2    NA
# 12     3    24     1     1  1000     1     2    NA
# 13     3    48     1     1  1000     1     2    NA
# 14     3    49     0     2     0     1     2    NA
# 15     3    50     0     2     0     1     2    NA


test0<-'
 $PARAM COV = 1

 $CMT EV CENT

 $TABLE
 double DV = (CENT / V2);

 $MAIN
 double CL = 1 ;
 if (COV == 2) CL = 2 ;
 double V2 = 5 ;
 double K12 = .05 ;
 double K20 = CL / V2 ;
 ALAG_EV = .1 ;

 $ODE
 dxdt_EV = - K12 * EV ;
 dxdt_CENT = K12 * EV - K20 * CENT  ;

 $CAPTURE DV
'
modtest<-mcode("test001",test0)

mrg0 <- mrgsim(modtest, datatest)
mrgh001 <- mrgsim(modtest, datatest, hmax = 0.01)
mrgevid2 <- datatest %>% 
  mutate(evid = ifelse(evid == 0, 2, evid)) %>% 
  mrgsim(x = modtest, end = -1)


datatest %>% 
  select(ID, time, evid, amt, COV) %>% 
  mutate(
    
    MRG0 = mrg0$DV, 
    MRG_h001 = mrgh001$DV, 
    MRG_EVID2 = mrgevid2$DV
  ) %>% 
  print.data.frame()
#    ID time evid  amt COV     MRG0 MRG_h001 MRG_EVID2
# 1   1    0    1 1000   1  0.00000  0.00000   0.00000
# 2   1   24    1 1000   1 19.62053 19.62053  19.62053
# 3   1   48    1 1000   1 25.69410 25.69410  25.69410
# 4   1   49    0    0   2 27.76512 27.76512  27.76512
# 5   1   50    0    0   2 29.29141 29.29141  29.29141
# 6   2    0    1 1000   1  0.00000  0.00000   0.00000
# 7   2   24    1 1000   1 19.62053 19.62053  19.62053
# 8   2   48    1 1000   1 25.69410 25.69410  25.69410
# 9   2   49    0    0   1 32.56425 32.56425  32.56425
# 10  2   50    0    0   1 38.41543 38.41543  38.41543
# 11  3    0    1 1000   2  0.00000  0.00000   0.00000
# 12  3   24    1 1000   2  8.64667  8.64667   8.64667
# 13  3   48    1 1000   2 11.25160 11.25160  11.25160
# 14  3   49    0    0   2 18.08403 18.08403  18.08403
# 15  3   50    0    0   2 22.80198 22.80198  22.80198