library(tidyverse)
library(haven)
library(glue)
library(gallimaufr)
library(ridittools)
library(broom)
library(summarytools)
library(magrittr)
library(officer)
library(flextable)
library(furrr)
library(tictoc)
library(childsds)
library(mice)
library(estimatr)
library(quantreg)

rm(list = ls())

options(future.globals.maxSize = 891289600)

# 1. Load Data ----
load("Data/df_raw.Rdata")
load("Data/mice.Rdata")

imp_long <- map(imp,
                ~ complete(.x, "long", TRUE) %>%
                  as_tibble() %>%
                  rename(imp = .imp) %>%
                  select(-.id, -matches("bmi_"), -matches("_report"),
                         -father_edu_level))
rm(imp)


# 2. Add Variables ----
to_ridit <- function(x){
  x_dense <- dense_rank(x)
  
  x_ridit <- toridit(table(x_dense))[x_dense]
  
  as.numeric(x_ridit)
}

df_inv <- map_dfr(imp_long, select, -matches("_..$"))

df_height <- map_dfr(imp_long, 
                     ~ .x %>%
                       select(id, cohort, imp, male, survey_weight, matches("(height|age)_")) %>%
                       pivot_longer(matches("(height|age)_")) %>%
                       separate(name, c("variable", "fup"), sep = "_") %>%
                       pivot_wider(names_from = variable, values_from = value)) %>%
  mutate(height_chart = sds(height, age, male, "height",
                            uk1990.ref, male = 1, female = 0)) %>%
  group_by(cohort, imp, male, fup) %>%
  mutate(height_ridit = to_ridit(height),
         height_z = wtd_scale(height, survey_weight),
         height_rank = percent_rank(height)) %>%
  ungroup() %>%
  rename(height_raw = height) %>%
  select(-male, -survey_weight, -age)

df_cog <- map_dfr(imp_long, 
                  ~ .x %>%
                    select(id, cohort, imp, survey_weight, matches("_resid_")) %>%
                    pivot_longer(matches("_resid_"), 
                                 names_to = "test", 
                                 values_to = "cog")) %>%
  group_by(cohort, imp, test) %>%
  mutate(cog_ridit = to_ridit(cog),
         cog_z = wtd_scale(cog, survey_weight),
         cog_rank = percent_rank(cog)) %>%
  ungroup() %>%
  select(-cog, -survey_weight) %>%
  filter(!str_detect(test, "comp"))

save(df_inv, df_cog, df_height,
     file = "Data/df_analysis.Rdata")

rm(imp_long)
gc()

# 3. Join Datasets ----
load("Data/df_analysis.Rdata")

keys <- df_cog %>%
  group_by(cohort) %>%
  group_keys() %>%
  pull(1)

split_df <- function(df, names, values){
  df %>%
    group_split(cohort, .keep = FALSE) %>%
    map(~ .x %>%
          pivot_wider(names_from = all_of(!!names), values_from = matches(!!values)) %>%
          arrange(id, imp)) %>%
    set_names(keys)
}

df_all <- df_cog %>%
  mutate(test = str_replace(test, "_resid", "")) %>%
  split_df("test", "^cog_")

df_all <- split_df(df_height, "fup", "^height_") %>%
  map(~ .x %>% select(-id, -imp)) %>%
  map2(df_all, ., bind_cols)

df_all <- df_inv %>%
  group_split(cohort, .keep = FALSE) %>%
  map(~ .x %>% arrange(id, imp) %>% select(-id, -imp)) %>%
  set_names(keys) %>%
  map2(df_all, ., bind_cols)

rm(df_cog, df_height, df_inv, keys, split_df)

save(df_all, file = "Data/df_all.Rdata")


# 4. Main Regression Functions ----
load("Data/df_all.Rdata")
load("Data/df_raw.Rdata")

mod_covars <- list(basic = "1",
                   sep = c("father_class", "mother_edu_level"),
                   height = c("father_class", "mother_edu_level", 
                              "mother_height", "father_height")) %>%
  map(glue_collapse, " + ")

mod_sex <- list(all = c(0, 1), male = 1, female = 0)

get_combn <- function(cohort){
  nms <- names(df_all[[cohort]])
  
  expand_grid(cog_score = str_subset(nms, "^cog_"),
              height_score = str_subset(nms, "^height_"))
}

mod_specs <- tibble(cohort = names(df_all)) %>%
  mutate(combn = map(cohort, get_combn)) %>%
  unnest(combn) %>%
  expand_grid(sex = names(mod_sex),
              mod = names(mod_covars),
              type = c("cc", "mi")) %>%
  filter(str_sub(cog_score, -2) == str_sub(height_score, -2),
         !str_detect(cog_score, "_comp_")) %>%
  mutate(spec_id = row_number(), .before = 1)

get_df <- function(spec){
  df_all[[spec$cohort]] %>%
    select(id, imp, survey_weight,
           cog_score = all_of(!!spec$cog_score),
           height_score = all_of(!!spec$height_score),
           male, mother_edu_level, father_class,
           mother_height, father_height) %>%
    filter(male %in% mod_sex[[!!spec$sex]])
}


get_result <- function(spec_id){
  spec <- slice(mod_specs, !!spec_id)
  
  df_reg <- get_df(spec)
  
  mod_form <- paste("height_score ~ cog_score", 
                    mod_covars[[spec$mod]],
                    sep = " + ")
  if (spec$sex == "all") mod_form <- glue("{mod_form} + male")
  
  run_mod <- function(df) lm_robust(as.formula(mod_form), df, survey_weight)
  
  if (spec$type == "cc"){
    
    mod <- df_reg %>% 
      filter(imp == 0) %>%
      run_mod()
    
  } else{
    
    mod <- df_reg %>% 
      filter(imp > 0) %>%
      group_split(imp) %>%
      map(run_mod) %>%
      pool()
    
  }
  
  tidy(mod, conf.int = TRUE) %>%
    filter(term == "cog_score") %>%
    select(beta = estimate, se = std.error,
           lci = conf.low, uci = conf.high)
}


# 4. Main Regression Results ----
mod_dict <- c(basic = "Sex-Adjusted",
              sep = "+ SEP",
              height = "+ Parental Height")

set.seed(1)
tic()
plan(multisession, workers = 4)
res_main <- mod_specs %>%
  sample_frac() %>%
  mutate(future_map_dfr(spec_id, get_result, .progress = TRUE)) %>%
  arrange(spec_id)
future:::ClusterRegistry("stop")
toc()

res_main <- res_main %>%
  separate(cog_score, into = c("cog", "cog_score", "test", "fup"), sep = "_") %>%
  select(-cog) %>%
  mutate(cohort_clean = case_when(cohort == "1970c_measure" ~ "1970c (Measured)",
                                  cohort == "2001c_white" ~ "2001c (White Only)",
                                  TRUE ~ cohort),
         mod_clean = factor(mod_dict[mod], mod_dict),
         fup_clean = ifelse(as.integer(fup) < 14, 11, 16),
         test_clean = glue("{str_to_title(test)} @ age {fup_clean}"),
         cog_clean = str_to_title(cog_score),
         height_score = str_sub(height_score, 8, -4),
         height_clean = str_to_title(height_score),
         string = glue("{round(beta, 2)} ({round(lci, 2)}, {round(uci, 2)})"),
         across(c(cog_clean, height_clean), ~ ifelse(.x == "Z", "Z-Score", .x)))

save(res_main, file = "Data/main_results.Rdata")


# 5. MNAR Regression Functions ----
mnar_specs <- mod_specs %>%
  select(-spec_id, -type) %>%
  distinct() %>%
  filter(mod == "basic", sex == "all") %>%
  filter(str_detect(cog_score, "^cog_ridit"), 
         str_detect(height_score, "^height_chart")) %>%
  expand_grid(diff = seq(from = -1, to = 1, length.out = 21)) %>%
  mutate(spec_id = row_number(), .before = 1)

get_mnar <- function(spec_id){
  spec <- slice(mnar_specs, !!spec_id)
  
  df_reg <- get_df(spec) %>%
    group_by(id) %>%
    mutate(miss = max(is.na(height_score))) %>%
    ungroup() %>%
    mutate(height_mnar = ifelse(miss == 1, height_score + spec$diff, height_score)) %>% 
    filter(imp > 0) %>%
    group_split(imp)
  
  mod_form <- paste("height_mnar ~ cog_score", 
                    mod_covars[[spec$mod]],
                    sep = " + ")
  if (spec$sex == "all") mod_form <- glue("{mod_form} + male")
  
  run_mod <- function(df) lm_robust(as.formula(mod_form), df, survey_weight)
  
  map(df_reg, run_mod) %>%
    pool() %>%
    tidy(conf.int = TRUE) %>%
    filter(term == "cog_score") %>%
    select(beta = estimate, se = std.error,
           lci = conf.low, uci = conf.high)
}


# 5. MNAR Regression Results ----
set.seed(1)
tic()
plan(multisession, workers = 4)
res_mnar <- mnar_specs %>%
  sample_frac() %>%
  mutate(future_map_dfr(spec_id, get_mnar, .progress = TRUE)) %>%
  arrange(spec_id)
future:::ClusterRegistry("stop")
toc()

get_diff <- function(beta_x, se_x, beta_y, se_y){
  beta <- beta_y - beta_x
  z <- beta / sqrt((se_x^2)+(se_y^2))
  p <- pnorm(z)
  
  tibble(beta = beta, p = p)
}

res_heat <- res_mnar %>%
  separate(cog_score, into = c("cog", "cog_score", "test", "fup"), sep = "_") %>%
  mutate(fup = ifelse(as.numeric(fup) < 14, 11, 16)) %>%
  unite("test", test, fup) %>%
  select(cohort, test, diff, beta, se)

res_heat <- expand_grid(rename_with(res_heat, ~ glue("{.x}_x")),
                        rename_with(res_heat, ~ glue("{.x}_y"))) %>%
  filter(test_x == test_y) %>%
  # filter(cohort_x < cohort_y) %>%
  mutate(get_diff(beta_x, se_x, beta_y, se_y),
         across(c(diff_x, diff_y), factor)) %>%
  select(-beta_x, -beta_y, -se_x, -se_y, -test_y) %>%
  rename(test = test_x)

save(res_mnar, res_heat, file = "Data/mnar_results.Rdata")


# 6. Difference in Height at Follow-Up ----
get_att_df <- function(cohort){
  df_all[[cohort]] %>%
    filter(imp == 0) %>%
    select(id, matches("^height_chart")) %>%
    rename_with(~ str_replace(.x, "^height_chart_", "")) %>%
    pivot_longer(-id, names_to = "fup", values_to = "height_chart") %>%
    mutate(miss = ifelse(is.na(height_chart), 1, 0),
           fup = as.integer(fup)) %>%
    left_join(df_all[[cohort]] %>%
                filter(imp == 0) %>%
                select(id, male, survey_weight),
              by = "id")
}

df_att <- names(df_all) %>%
  set_names(., .) %>%
  map_dfr(get_att_df, .id = "cohort")

att_specs <- expand_grid(
  df_att %>%
    distinct(cohort, fup) %>%
    rename_with(~ glue("{.x}_x")),
  df_att %>%
    distinct(cohort, fup) %>%
    rename_with(~ glue("{.x}_y")),
  male = c(0, 1),
) %>%
  filter(cohort_x == cohort_y,
         str_length(cohort_x) == 5,
         fup_x < fup_y) %>%
  select(cohort = cohort_x, fup_height = fup_x, fup_miss = fup_y, male) %>%
  mutate(spec_id = row_number(), .before = 1)

get_att <- function(spec_id){
  spec <- slice(att_specs, !!spec_id)
  
  df_h <- df_att %>%
    filter(cohort == !!spec$cohort,
           fup == !!spec$fup_height,
           male == !!spec$male) %>%
    select(id, height_chart, survey_weight)
  
  df_m <- df_att %>%
    filter(cohort == !!spec$cohort,
           fup == !!spec$fup_miss,
           male == !!spec$male) %>%
    select(id, miss)
  
  df_reg <- left_join(df_h, df_m, by = "id")
  
  lm_robust(height_chart ~ miss, df_reg, survey_weight) %>%
    tidy(conf.int = TRUE) %>%
    filter(term == "miss") %>%
    select(beta = estimate, lci = conf.low, uci = conf.high)
}

res_att <- att_specs %>%
  mutate(map_dfr(spec_id, get_att)) %>%
  mutate(sex = factor(male, labels = c("Female", "Male")),
         fup_height = ifelse(fup_height < 10, 
                             glue("0{fup_height}"),
                             as.character(fup_height)),
         cohort_clean = glue("{cohort}: Height @ Age {fup_height}"),
         fup_miss = factor(fup_miss))

save(res_att, file = "Data/attrition_results.Rdata")


# 7. Quantile Regressions ----
qreg_specs <- mod_specs %>%
  select(-type, -spec_id) %>%
  filter(sex == "all",
         str_detect(cog_score, "cog_ridit"),
         str_detect(height_score, "height_chart"),
         str_length(cohort) == 5) %>%
  distinct() %>%
  mutate(spec_id = row_number(), .before = 1)

get_qreg <- function(spec_id){
  spec <- slice(qreg_specs, !!spec_id)
  
  df_reg <- get_df(spec) %>%
    select(imp, height_score, cog_score, male, father_class,
           mother_edu_level, mother_height, father_height, survey_weight) %>%
    group_split(imp, .keep = FALSE)
  
  mod_form <- paste("height_score ~ cog_score + male +", mod_covars[[spec$mod]])
  
  get_boots <- function(df_imp){
    form <- paste(mod_form, "+ survey_weight")
    
    df_frame <- df_imp %>%
      model.frame(as.formula(form), .)
    
    map_dfr(1:200, get_boot, df_frame)
  }
  
  get_boot <- function(boot, df_frame){
    set.seed(boot)
    df_boot <- sample_frac(df_frame, replace = TRUE)
    
    mod <- rq(as.formula(mod_form), 1:9/10,
              df_boot, weights = survey_weight)
    
    tibble(tau = 1:9/10, estimate = set_names(coef(mod)[2, ], NULL))
  }
  
  get_ci <- function(estimate){
    quantile(estimate, c(.5, .025, .975)) %>%
      set_names(c("beta", "lci", "uci")) %>%
      as_tibble_row()
  }
  
  map_dfr(df_reg, get_boots, .id = "type") %>%
    mutate(type = ifelse(type == "1", "cc", "mi")) %>%
    group_by(type, tau) %>%
    summarise(get_ci(estimate),
              .groups = "drop")
}

set.seed(1)
tic()
plan(multisession, workers = 2)
res_qreg <- qreg_specs %>%
  sample_frac() %>%
  # slice(1:2) %>%
  mutate(res = future_map(spec_id, get_qreg, .progress = TRUE)) %>%
  unnest(res) %>%
  arrange(spec_id)
future:::ClusterRegistry("stop")
toc()

res_qreg <- res_qreg %>%
  mutate(test = str_sub(cog_score, 11, -4),
         fup = str_sub(cog_score, -2) %>% as.numeric(),
         fup_clean = ifelse(fup < 14, 11, 16),
         test_clean = glue("{str_to_title(test)} @ Age {fup_clean}"))

save(res_qreg, file = "Data/quantile_results.Rdata")
