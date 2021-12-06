library(tidyverse)
library(mice)
library(gallimaufr)
library(estimatr)
library(ridittools)
library(glue)
library(splines)
library(ggeffects)
library(furrr)
library(tictoc)
library(broom)

rm(list = ls())

options(future.globals.maxSize = 891289600)

# 1. Load Data ----
load("Data/mice.Rdata")
load("Data/df_cog.Rdata")

imp_long <- map(imp,
                ~ .x %>%
                  complete("long", include = TRUE) %>%
                  as_tibble() %>%
                  select(-.id) %>%
                  rename(imp = .imp))

df_bmi <- map_dfr(imp_long,
                  ~ .x %>%
                    filter(imp == 0) %>%
                    select(id, matches("^bmi_")) %>%
                    pivot_longer(-id, names_to = "bmi_var", 
                                 values_to = "bmi_value"),
                  .id = "cohort")

nest_list <- function(df){
  split(df, df$cohort) %>%
    map(~ split(.x[[3]], .x[[2]]))
}

sample_list <- df_bmi %>%
  count(cohort, id, wt = is.na(bmi_value), name = "n_miss") %>%
  uncount(2, .id = "sample") %>%
  filter(sample == 1 | n_miss == 0) %>%
  mutate(sample = ifelse(sample == 1, "all", "balanced")) %>%
  select(cohort, sample, id) %>%
  nest_list()

miss_list <- df_bmi %>%
  filter(is.na(bmi_value)) %>%
  select(cohort, bmi_var, id) %>%
  nest_list()

bmi_sd <- df_bmi %>%
  group_by(cohort, bmi_var) %>%
  summarise(sd = sd(bmi_value, na.rm = TRUE),
            .groups = "drop") %>%
  nest_list()

rm(imp, df_bmi)


# 2. Simple Linear Regression ----
# Model Parameters
sexes <- list(female = 0, male = 1, all = 0:1)

types <- list(cc = 0, mi = 1:64)

mods <- lst(basic = "1",
            sep = c("father_class", "mother_edu_level"),
            early = c("birth_weight", "childhood_bmi"), # "maternal_age", 
            all = c(early, sep))

covars <- c("male", mods[["all"]])

mod_specs <- expand_grid(cohort = names(imp_long),
                         type = c("cc", "mi"),
                         sample = c("all", "balanced"),
                         sex = names(sexes),
                         mod = names(mods),
                         bmi_score = c("bmi_raw", "bmi_std", "bmi_rank"),
                         cog_score = c("cog_ridit", "cog_rank", "cog_std")) %>%
  mutate(bmi_var = map(cohort, 
                       ~ imp_long[[.x]] %>%
                         names() %>%
                         str_subset("^bmi_"))) %>%
  unnest(bmi_var) %>%
  mutate(cog_var = map(cohort, 
                       ~ imp_long[[.x]] %>%
                         names() %>%
                         str_subset("^(maths|verbal|vocab|g)_"))) %>%
  unnest(cog_var)

mod_specs <- mod_specs %>%
  filter(str_detect(cog_var, "g_all"),
         cog_score == "cog_std",
         bmi_score %in% c("bmi_raw", "bmi_std")) %>%
  mutate(spec_id = row_number(),
         .before = 1)


# Model Functions
to_ridit <- function(x){
  x_dense <- dense_rank(x)
  
  x_ridit <- toridit(table(x_dense))[x_dense]
  
  as.numeric(x_ridit)
}

get_scores <- function(x, stub){
  tibble(raw = x,
         std = wtd_scale(x),
         rank = percent_rank(x),
         ridit = to_ridit(x)) %>%
    rename_with(~ glue("{stub}_{.x}"))
}

prepare_df <- function(spec_id){
  spec <- slice(mod_specs, !!spec_id)
  
  imp_long[[spec$cohort]] %>%
    filter(male %in% sexes[[spec$sex]],
           imp %in% types[[spec$type]]) %>%
    rename(bmi = all_of(spec$bmi_var),
           cog = all_of(spec$cog_var)) %>%
    group_by(imp) %>%
    mutate(get_scores(cog, "cog"),
           get_scores(bmi, "bmi")) %>%
    ungroup() %>%    
    select(imp, id, survey_weight,
           bmi = all_of(spec$bmi_score),
           cog = all_of(spec$cog_score),
           all_of(covars)) %>%
    filter(id %in% sample_list[[spec$cohort]][[spec$sample]])
}

get_form <- function(spec_id, cog = "cog", spec_df = mod_specs){
  spec <- slice(spec_df, !!spec_id)
  
  ind_vars <- ifelse(spec$sex == "all", "male", "1") %>%
    c(mods[[spec$mod]]) %>%
    glue_collapse(" + ")
  
  glue("bmi ~ {cog} + {ind_vars}")
}

get_pool <- function(spec_id){
  spec <- slice(mod_specs, !!spec_id)
  
  df <- prepare_df(spec_id)
  if (spec$type == "cc"){
    df <- uncount(df, 2, .id = "imp")
  }
  
  mod_form <- get_form(spec_id) %>% as.formula()
  
  group_split(df, imp) %>%
    map(~ lm_robust(mod_form, .x, weights = survey_weight)) %>%
    pool() %>%
    tidy(conf.int = TRUE) %>%
    filter(term == "cog") %>%
    select(beta = estimate, se = std.error, 
           pval = p.value, 
           lci = conf.low, uci = conf.high)
}


# Run Models
set.seed(1)
tic()
plan(multisession)
ols_res <- mod_specs %>%
  sample_frac() %>%
  mutate(future_map_dfr(spec_id, get_pool, .progress = TRUE)) %>%
  arrange(spec_id)
future:::ClusterRegistry("stop")
toc()

save(ols_res, file = "Data/ols_res.Rdata")


# 3. Splines ----
# Model Function
get_splines <- function(spec_id){
  spec <- slice(mod_specs, !!spec_id)
  
  df <- prepare_df(spec_id)
  if (spec$type == "cc"){
    df <- uncount(df, 2, .id = "imp")
  }
  
  get_form <- function(cog){
    ind_vars <- ifelse(spec$sex == "all", "male", "1") %>%
      c(mods[[spec$mod]]) %>%
      glue_collapse(" + ")
    
    glue("bmi ~ {cog} + {ind_vars}")
  }
  
  get_ggeffect <- function(data, mod_form){
    mod <- mod_form %>%
      as.formula() %>%
      lm_robust(data, weights = survey_weight)
    
    ggeffect(mod, "cog [-3:3 by=0.1]", data = data)
  }
  
  get_effects <- function(mod_form){
    group_split(df, imp) %>%
      map(get_ggeffect, mod_form) %>%
      pool_predictions() %>%
      as_tibble() %>%
      select(-group)
  }
  
  bind_rows(linear = get_form("cog") %>% get_effects(),
            splines = get_form("ns(cog, 2)") %>% get_effects(),
            .id = "mod_form")
}

# Run Models
set.seed(1)
tic()
plan(multisession)
spline_res <- mod_specs %>%
  sample_frac() %>%
  mutate(res = future_map(spec_id, get_splines, .progress = TRUE)) %>%
  unnest(res) %>%
  arrange(spec_id)
future:::ClusterRegistry("stop")
toc()

save(spline_res, file = "Data/spline_res.Rdata")


# 4. MAR Pattern-Mixture ----
mnar_specs <- mod_specs %>%
  filter(type == "mi", mod == "basic", sex == "all", cog_score == "cog_std",
         bmi_score == "bmi_raw", sample == "all") %>%
  expand_grid(mnar = seq(from = -1, to = 1, by = 0.1)) %>%
  mutate(spec_id = row_number())

prepare_mnar <- function(spec_id){
  spec <- slice(mnar_specs, !!spec_id)
  
  imp_long[[spec$cohort]] %>%
    filter(male %in% sexes[[spec$sex]],
           imp %in% types[[spec$type]]) %>%
    rename(bmi = all_of(spec$bmi_var),
           cog = all_of(spec$cog_var)) %>%
    mutate(bmi = ifelse(id %in% miss_list[[spec$cohort]][[spec$bmi_var]],
                        bmi + spec$mnar * bmi_sd[[spec$cohort]][[spec$bmi_var]],
                        bmi)) %>%
    group_by(imp) %>%
    mutate(get_scores(cog, "cog"),
           get_scores(bmi, "bmi")) %>%
    ungroup() %>%    
    select(imp, id, survey_weight,
           bmi = all_of(spec$bmi_score),
           cog = all_of(spec$cog_score),
           all_of(covars)) %>%
    filter(id %in% sample_list[[spec$cohort]][[spec$sample]])
}

get_mnar <- function(spec_id){
  spec <- slice(mnar_specs, !!spec_id)
  
  df <- prepare_mnar(spec_id)
  if (spec$type == "cc"){
    df <- uncount(df, 2, .id = "imp")
  }
  
  mod_form <- get_form(spec_id, spec_df = mnar_specs) %>% as.formula()
  
  group_split(df, imp) %>%
    map(~ lm_robust(mod_form, .x, weights = survey_weight)) %>%
    pool() %>%
    tidy(conf.int = TRUE) %>%
    filter(term == "cog") %>%
    select(beta = estimate, se = std.error, 
           pval = p.value, 
           lci = conf.low, uci = conf.high)
}

# Run Models
set.seed(1)
tic()
plan(multisession)
mnar_res <- mnar_specs %>%
  slice(1:5) %>%
  mutate(map_dfr(spec_id, get_mnar)) %>%
  # sample_frac() %>%
  # mutate(future_map_dfr(spec_id, get_mnar, .progress = TRUE)) %>%
  arrange(spec_id)
future:::ClusterRegistry("stop")
toc()

save(mnar_res, file = "Data/mnar_res.Rdata")

