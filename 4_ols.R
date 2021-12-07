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
                  dplyr::select(-.id) %>%
                  rename(imp = .imp))

df_bmi <- map_dfr(imp_long,
                  ~ .x %>%
                    filter(imp == 0) %>%
                    dplyr::select(id, matches("^bmi_")) %>%
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
  dplyr::select(cohort, sample, id) %>%
  nest_list()

miss_list <- df_bmi %>%
  filter(is.na(bmi_value)) %>%
  dplyr::select(cohort, bmi_var, id) %>%
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

mods <- lst(basic = c(),
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
  filter(sample == "all" | str_detect(cog_var, "^g_all"),
         bmi_score == "bmi_raw" | str_detect(cog_var, "^g_all"),
         sex == "all" | str_detect(cog_var, "^g_all")) %>%
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
    dplyr::select(imp, id, survey_weight,
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
    dplyr::select(beta = estimate, se = std.error, 
           pval = p.value, 
           lci = conf.low, uci = conf.high)
}


# Run Models
set.seed(1)
tic()
plan(multisession, workers = 6)
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
  
  get_ggeffect <- function(data, mod_form){
    df_mod <- data
    
    mod <- mod_form %>%
      as.formula() %>%
      lm_robust(df_mod, weights = survey_weight)
    
    # ggeffect(mod, "cog [-3:3 by=0.1]", data = data)
    ggpredict(mod, "cog [-3:3 by=0.05]", data = df_mod)
  }
  
  get_effects <- function(mod_form){
    group_split(df, imp) %>%
      map(get_ggeffect, mod_form) %>%
      pool_predictions() %>%
      as_tibble() %>%
      dplyr::select(-group)
  }
  
  bind_rows(linear = get_form(spec_id, "cog") %>% get_effects(),
            splines = get_form(spec_id, "splines::ns(cog, 2)") %>% get_effects(),
            .id = "mod_form")
}

# Run Models
set.seed(1)
tic()
plan(multisession, workers = 6)
spline_res <- mod_specs %>%
  filter(cog_score != "cog_ridit",
         str_detect(cog_var, "^g_(all|closer)"),
         sex == "all", sample == "all", 
         bmi_score == "bmi_raw") %>%
  sample_frac() %>%
  mutate(res = future_map(spec_id, get_splines, .progress = TRUE)) %>%
  unnest(res) %>%
  arrange(spec_id)
future:::ClusterRegistry("stop")
toc()

save(spline_res, file = "Data/spline_res.Rdata")


# 4. MAR Pattern-Mixture ----
mnar_specs <- mod_specs %>%
  filter(type == "mi", mod == "basic", sex == "all",
         str_detect(cog_var, "g_(all|closer)"),
         cog_score == "cog_std", bmi_score == "bmi_raw",
         sample == "all") %>%
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
    dplyr::select(imp, id, survey_weight,
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
    dplyr::select(beta = estimate, se = std.error, 
           pval = p.value, 
           lci = conf.low, uci = conf.high)
}

# Run Models
set.seed(1)
tic()
plan(multisession, workers = 6)
mnar_res <- mnar_specs %>%
  # slice(1:5) %>%
  # mutate(map_dfr(spec_id, get_mnar)) %>%
  sample_frac() %>%
  mutate(future_map_dfr(spec_id, get_mnar, .progress = TRUE)) %>%
  arrange(spec_id)
future:::ClusterRegistry("stop")
toc()

save(mnar_res, file = "Data/mnar_res.Rdata")


# 5. GAMLSS ----
item_obs <- imp_long %>%
  map(~ .x %>%
        filter(imp == 0) %>%
        dplyr::select(id, male, all_of(mods[["all"]])) %>%
        drop_na() %>%
        pull(id))

gamlss_specs <- mod_specs %>%
  filter(type == "cc",
         str_detect(cog_var, "g_(all|closer)"),
         bmi_score == "bmi_raw",
         cog_score == "cog_std") %>%
  expand_grid(item_type = c("all", "obs"))

prepare_df_gamlss <- function(spec_id, item_type = "all"){
  spec <- slice(mod_specs, !!spec_id)
  
  df <- prepare_df(spec_id)
  
  if (item_type == "obs"){
    df <- df %>%
      filter(id %in% item_obs[[spec$cohort]])
  }
  
  df %>%
    dplyr::select(id, survey_weight, bmi, cog, male,
           all_of(mods[[spec$mod]])) %>%
    drop_na()
}

get_gamlss <- function(spec_id, item_type = "all"){
  spec <- slice(mod_specs, !!spec_id)
  
  df <- prepare_df_gamlss(spec_id, item_type)
  
  mod_form <- get_form(spec_id) %>% str_replace("bmi ", "") %>% as.formula()
  main_form <- get_form(spec_id) %>% as.formula()
  
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
  
  set.seed(1)
  boots <- map_dfr(1:200, 
                   ~ df %>%
                     sample_frac(replace = TRUE) %>%
                     run_mod()) %>% 
    chop(mu:nu) 
  
  tibble(res = list(main), boots = list(boots), n = nrow(df))
}

# Run Models
set.seed(1)
tic()
plan(multisession, workers = 6)
gamlss_res <- gamlss_specs %>%
  # mutate(map_dfr(spec_id, get_mnar)) %>%
  sample_frac() %>%
  mutate(future_map2_dfr(spec_id, item_type, get_gamlss, .progress = TRUE)) %>%
  arrange(spec_id)
future:::ClusterRegistry("stop")
toc()

save(gamlss_res, file = "Data/gamlss_res.Rdata")
