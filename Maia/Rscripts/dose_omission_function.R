library(tidyr)
library(extraDistr)

dose_omission <- function(initial_data, cohort, proportion, omission, n_unif_omis, min_unif_omis, max_unif_omis, seed_unif_omis) {
  set.seed(seed_unif_omis)

  data_1 <- initial_data %>% #' initial_data is from Scenario_1
    dplyr::filter(coh == cohort) %>%
    dplyr::group_by(ID) %>%
    mutate(TIME_unif_omis = round(x = runif(
      n = n_unif_omis,
      min = min_unif_omis,
      max = max_unif_omis
    ), digits = 3)) %>%
    ungroup()


  #' sample the proportion of subjects per cohort to receive the dose_omission
  id_proportion <- proportion * length(unique(data_1$ID))
  nsel <- ifelse(id_proportion - floor(id_proportion) < 0.5, floor(id_proportion), ceiling(id_proportion))

  mysel_omission <- sample(data_1$ID, size = nsel, replace = F)

  data_omis <- data_1 %>%
    group_by(ID) %>%
    mutate(
      IDflag_omis = ifelse(ID %in% mysel_omission, 1, 0), # flag for the suject in question
      TIMEflag_omis = ifelse(IDflag_omis == 1 & TIME >= TIME_unif_omis, 1, 0),
      DOSEflag_omis = {
        TIMEflag_omis == 1
      } %>%
        {
          . * !duplicated(.)
        },
      DOSE_omis = ifelse(DOSEflag_omis == 1, omission, DOSE)
    ) %>%
    ungroup()

  for (i in unique(data_omis$ID)) {
    try1 <- filter(data_omis, ID == i & coh == 5) %>%
      pull(DOSE_omis)
    if (any(try1 == 0.001)) {
      place <- which(try1 == 0.001)
      if (length(try1) > place) {
        try1[place + 1] <- 0.001
      }
    }
    data_omis[data_omis$ID == i & data_omis$coh == 5, "DOSE_omis"] <- try1
  }

  data_omis
}



data_bind_omis <- rbind(
  dose_omission(
    initial_data = full_data,
    cohort = 3,
    proportion = 0.10,
    omission = 0.001,
    n_unif_omis = 1, min_unif_omis = 7 * 12 / 356, max_unif_omis = 7 * 24 / 356, seed_unif_omis = 123
  ),
  dose_omission(
    initial_data = full_data,
    cohort = 4,
    proportion = 0.25,
    omission = 0.001,
    n_unif_omis = 1, min_unif_omis = 7 * 12 / 356, max_unif_omis = 7 * 24 / 356, seed_unif_omis = 124
  ),
  dose_omission(
    initial_data = full_data,
    cohort = 5,
    proportion = 0.50,
    omission = 0.001,
    n_unif_omis = 1, min_unif_omis = 7 * 12 / 356, max_unif_omis = 7 * 24 / 356, seed_unif_omis = 125
  )
)
