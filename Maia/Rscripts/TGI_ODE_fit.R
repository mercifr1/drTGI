library(rstan)
data_for_fit<-design_flat_nadir%>%
  filter(REP==1)%>%
  mutate(TIME=ifelse(TIME<0.1,0,TIME))



mod<-stan_model("name_of_stan_file")


data <- list( N = nrow(data_for_fit),
              J = length(unique(data_for_fit$ID)),
              y = data_for_fit$DV,
              t0 = 0.0, 
              num_obs_Ind = data_for_fit %>% group_by(ID) %>% summarise("maxx" = max(N_TIME)) %>% select(maxx) %>% pull(),
)

