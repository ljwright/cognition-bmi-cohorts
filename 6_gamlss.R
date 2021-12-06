library(tidyverse)
library(gamlss)
library(gallimaufr)
library(ridittools)
library(glue)
library(tictoc)
library(furrr)

rm(list = ls())


# 1. Load Data ----
load("Data/df_cog.Rdata")

to_ridit <- function(x){
  x_dense <- dense_rank(x)
  
  x_ridit <- toridit(table(x_dense))[x_dense]
  
  as.numeric(x_ridit)
}

df_bmi <- map_dfr(df_cog, 
                  ~ .x %>%
                    dplyr::select(id, matches("^bmi")) %>%
                    pivot_longer(-id) %>%
                    mutate(age = str_sub(name, -2) %>% as.integer(),
                           name = str_sub(name, 1, -4)),
                  .id = "cohort") %>%
  pivot_wider(names_from = name, values_from = value) %>%
  drop_na() %>%
  filter(age >= 20)

df_cov <- map_dfr(df_cog,
                  ~ .x %>%
                    dplyr::select(id, male, birth_weight,
                                  father_class, # any_of("maternal_age"),
                                  matches("bmi_(10|11)")) %>%
                    rename_with(~ "childhood_bmi", .cols = matches("bmi_(10|11)")) %>%
                    drop_na() %>%
                    group_by(male) %>%
                    mutate(childhood_bmi = wtd_scale(childhood_bmi)) %>%
                    ungroup(),
                  .id = "cohort")

df_cog <- map_dfr(df_cog,
                  ~ .x %>%
                    dplyr::select(id, matches("_resid$")) %>%
                    rename_with(~ str_replace(.x, "_resid", "")) %>%
                    pivot_longer(-id, values_to = "cog_score", names_to = "cog_var"),
                  .id = "cohort") %>%
  drop_na() %>%
  group_by(cohort, cog_var) %>%
  mutate(cog_std = wtd_scale(cog_score),
         cog_ridit = to_ridit(cog_score),
         cog_rank = percent_rank(cog_score)) %>%
  ungroup() %>%
  dplyr::select(-cog_score)

df_panel <- df_bmi %>%
  group_by(cohort) %>%
  complete(id, age) %>%
  arrange(cohort, id, age) %>%
  group_by(cohort, id) %>%
  mutate(last_age = ifelse(is.na(lead(bmi)) & !is.na(bmi), age, NA)) %>%
  filter(last_age == min(last_age, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(cohort, id, last_age)

gc()

# 2. Model Objects ----
sexes <- list(all = c(0, 1), female = 0, male = 1)

gamlss_forms <- list(basic = c(),
                  sep = "father_class",
                  all = c("birth_weight", "father_class", "childhood_bmi"))

gamlss_specs <- df_bmi %>%
  distinct(cohort, age) %>%
  left_join(., ., by = "cohort") %>%
  rename(dep_age = age.x, last_age = age.y) %>%
  group_by(cohort) %>%
  filter(last_age == min(last_age) | 
           last_age == max(last_age)) %>%
  ungroup() %>%
  left_join(df_cog %>% 
              pivot_longer(cog_std:cog_rank, names_to = "cog_score") %>%
              distinct(cohort, cog_var, cog_score),
            by = "cohort") %>%
  expand_grid(sex = c("all", "male", "female"),
              mod = names(gamlss_forms)) %>%
  group_by(cohort) %>%
  filter(sex == "all" |
           last_age == min(last_age)) %>%
  ungroup() %>%
  mutate(spec_id = row_number(), .before = 1)
 
get_spec <- function(spec_id){
  gamlss_specs %>%
    slice(!!spec_id)
}

make_df <- function(spec_id){
  spec <- get_spec(spec_id)
  
  df_panel %>%
    filter(cohort == !!spec$cohort,
           last_age >= !!spec$last_age) %>%
    dplyr::select(id) %>%
    inner_join(df_cov %>%
                 filter(male %in% sexes[[spec$sex]]) %>%
                 dplyr::select(-cohort),
               by = "id") %>%
    inner_join(df_cog %>%
                 filter(cohort == !!spec$cohort,
                        cog_var == !!spec$cog_var) %>%
                 dplyr::select(id, cog = all_of(spec$cog_score)),
               by = "id") %>%
    inner_join(df_bmi %>%
                 filter(cohort == !!spec$cohort,
                        age == !!spec$dep_age) %>%
                 dplyr::select(id, bmi),
               by = "id")
}

save(df_bmi, df_cog, df_cov, df_panel,
     gamlss_specs, sexes, gamlss_forms,
     get_spec, make_df,
     file = "Data/gamlss_objects.Rdata")


# 3. Run gamlss ----
load("Data/gamlss_objects.Rdata")

run_gamlss <- function(spec_id){
  spec <- get_spec(spec_id)
  
  df <- make_df(spec_id)
  
  mod_form <- c("~ cog", gamlss_forms[[spec$mod]]) %>%
    glue_collapse(" + ")
  
  if (spec$sex == "all") mod_form <- glue("{mod_form} + male")
  
  main_form <- glue("bmi {mod_form}") %>% as.formula()
  mod_form <- as.formula(mod_form)
  
  run_mod <- function(data){
    gamlss(main_form,
           sigma.formula = mod_form,
           nu.formula = mod_form,
           family = BCCG, 
           data = data,
           trace = FALSE) %>%
      coefAll() %>%
      map_dfr(~ enframe(.x, name = "term", value = "beta"), .id = "param") %>%
      pivot_wider(names_from = param, values_from = beta)
  }
  
  main <- run_mod(df)
  
  boots <- map_dfr(1:200, 
              ~ df %>%
                sample_frac(replace = TRUE) %>%
                run_mod()) %>% 
    chop(mu:nu) 
  
  tibble(res = list(main), boots = list(boots), n = nrow(df))
}

set.seed(1)
tic()
plan(multisession)
gamlss_res <- gamlss_specs %>%
  filter(str_detect(cog_var, "g_closer"), 
         cog_score != "cog_ridit")  %>%
  dplyr::select(spec_id) %>%
  sample_frac() %>%
  mutate(res = future_map(spec_id, run_gamlss, .progress = TRUE))
future:::ClusterRegistry("stop")
toc()

save(gamlss_res, file = "Data/gamlss_results.Rdata")

# 4. Plot Results ----






df_gamlss <- df_x %>%
  mutate(mod = map(data, get_gamlss)) %>%
  dplyr::select(-data, -test_resid) %>%
  unnest(mod)

df_gamlss %>%
  filter(test == "maths @ 11") %>%
  ggplot() +
  aes(x = centile, y = bmi, color = cog_score) +
  facet_wrap(cohort ~ age, scales = "free_y") +
  geom_line() +
  geom_hline(yintercept = c(18.5, 25, 30))

df_gamlss %>%
  arrange(cog_score) %>%
  group_by(test, cohort, age, centile) %>%
  summarise(bmi = diff(bmi),
            .groups = "drop") %>%
  filter(!str_detect(test, "^comp"),
         age >= 10) %>%
  ggplot() +
  aes(x = centile, y = bmi, color = age, group = age) +
  facet_grid(test ~ cohort) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed")

df_gamlss %>%
  arrange(cog_score) %>%
  group_by(test, cohort, age, centile) %>%
  summarise(bmi = diff(bmi),
            .groups = "drop") %>%
  filter(!str_detect(test, "^comp"),
         age >= 10) %>%
  filter(centile %in% c(10, 25, 50, 75, 90)) %>%
  ggplot() +
  aes(x = age, y = bmi, color = cohort, group = cohort) +
  facet_grid(test ~ centile) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed")