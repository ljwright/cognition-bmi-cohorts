library(tidyverse)
library(mice)
library(tictoc)

rm(list = ls())

# 1. Load Data ----
load("Data/df_cog.Rdata")

df_mice <- map(df_cog,
    ~ .x %>%
      rename(childhood_bmi = matches("bmi_(10|11)")) %>%
      select(-matches("_raw"), -cohort,
             -matches("(bmi|self_report)_(0|1).$")) %>%
      rename_with(~ str_replace(.x, "_resid", "")))

rm(df_cog)

# 2. Multiple Imputation ----
tic()
imp <- list()
for(cohort in names(df_mice)){
  df <- df_mice[[cohort]]
  
  pred <- make.predictorMatrix(df)
  pred[, c("id", "male")] <- 0
  
  temp <- list()
  for (male in c(0, 1)){
    df_m <- filter(df, male == !!male)
    
    temp[[male + 1]] <- parlmice(df_m, n.core = 4, n.imp.core = 16,
                                 predictorMatrix = pred,
                                 defaultMethod = c("rf", "logreg", "polyreg", "polr"))
  }
  
  imp[[cohort]] <- rbind(temp[[1]], temp[[2]])
}
toc()

rm(df, df_mice, temp, pred, male, cohort)

save(imp, file = "Data/mice.Rdata")

# 1. Coefficient at each age (repeat OLS)
# 2. Cog splines at each age (repeat OLS)
# 3. GCM
# 4. LGCA
# 5. GAMLSS at each age with cog splines (repeat)
