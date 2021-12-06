library(tidyverse)
library(gallimaufr)
library(ridittools)
library(glue)
library(tictoc)
library(furrr)
library(splines)
library(lme4)
library(broom)
library(splines)

rm(list = ls())

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

transform_age <- function(age, func, stub){
  func(age, 2) %>%
    as_tibble() %>%
    mutate(across(everything(), as.numeric)) %>%
    rename_with(~ glue("{stub}{.x}"))
}

df_s <- df_bmi %>%
  mutate(age = str_sub(bmi_var, -2) %>% as.integer()) %>%
  distinct(cohort, age) %>%
  group_by(cohort) %>%
  mutate(# transform_age(age, poly, "age_p"),
    transform_age(age, ns, "age_ns"),
    age_ns0 = 1) %>%
  arrange(cohort, age) %>%
  ungroup()

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

rm(df_bmi, imp)

# 2. Model Objects ----
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
                         cog_score = c("cog_ridit", "cog_rank", "cog_std")) %>%
  mutate(cog_var = map(cohort, 
                       ~ imp_long[[.x]] %>%
                         names() %>%
                         str_subset("^(maths|verbal|vocab|g)_"))) %>%
  unnest(cog_var) %>%
  filter(str_detect(cog_var, "^g_all")) %>%
  mutate(spec_id = row_number(),
         .before = 1)

# Functions
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
    rename(cog = all_of(spec$cog_var)) %>%
    group_by(imp) %>%
    mutate(get_scores(cog, "cog")) %>%
    ungroup() %>%    
    select(imp, id, survey_weight,
           matches("^bmi_"),
           cog = all_of(spec$cog_score),
           all_of(covars)) %>%
    filter(id %in% sample_list[[spec$cohort]][[spec$sample]]) %>%
    pivot_longer(matches("^bmi_"), names_to = "age", values_to = "bmi") %>%
    mutate(age = str_sub(age, -2) %>% as.integer()) %>%
    left_join(df_s %>%
                filter(cohort == spec$cohort),
              by = "age")
}

get_form <- function(spec_id){ # NEEDS TO CHANGE TO HAVE CORRECT COVARIATES - CAN'T ADD + 1
  spec <- slice(mod_specs, !!spec_id)
  
  ind_vars <- ifelse(spec$sex == "all", "male", "1") %>%
    c(mods[[spec$mod]]) %>%
    glue_collapse(" + ")
  
  glue("bmi ~ -1 + age_ns0 + age_ns1 + age_ns2 +
       cog + cog*age_ns1 + cog*age_ns2 +
       (-1 + age_ns0 + age_ns1 + age_ns2 | id)")
}

run_lmer <- function(spec_id){
  spec <- slice(mod_specs, !!spec_id)
  
  df <- prepare_df(spec_id)
  
  mod_form <- get_form(spec_id) %>% as.formula()
  
  get_glht <- function(pred_type, mod){
    reg_ex <- ifelse(pred_type == "margin", "cog", "(cog|age_ns)")
    
    run_glht <- function(spec){
      res <- multcomp::glht(mod, spec)
      
      bind_cols(tidy(confint(res)), se = tidy(res)$std.error) %>%
        select(beta = estimate, se = se, lci = conf.low, uci = conf.high)
    }
    
    df_s %>%
      filter(cohort == spec$cohort) %>%
      expand_grid(cog = -1:1,
                  term = names(fixef(mod))) %>%
      mutate(beta = ifelse(str_detect(term, !!reg_ex), 1, 0),
             beta = ifelse(str_detect(term, "age_ns1"), beta * age_ns1, beta),
             beta = ifelse(str_detect(term, "age_ns2"), beta * age_ns2, beta),
             beta = ifelse(str_detect(term, "cog"), beta * cog, beta)) %>%
      select(age, cog, term, beta) %>%
      group_by(age, cog) %>%
      filter(max(abs(beta)) != 0) %>%
      ungroup() %>%
      nest(spec = c(term, beta)) %>%
      mutate(spec = map(spec, ~ deframe(.x) %>% t()),
             res = map(spec, run_glht)) %>%
      unnest(res) %>%
      select(-spec)
  }
  
  get_preds <- function(data){
    mod <- lmer(mod_form, data)
    
    bind_rows(margin = get_glht("margin", mod),
              prediction = get_glht("prediction", mod),
              .id = "pred_type")
  }
  
  group_split(df, imp) %>%
    map_dfr(get_preds, .id = "imp")
}

spec_id <- 73

data <- df %>% filter(imp == min(imp))
mod <- lmer(mod_form, data)


get_glht <- function(pred_type){
  reg_ex <- ifelse(pred_type == "margin", "cog", "(cog|age_ns)")
  
  run_glht <- function(spec){
    res <- multcomp::glht(mod, spec)
    
    bind_cols(tidy(confint(res)), se = tidy(res)$std.error) %>%
      select(beta = estimate, se = se, lci = conf.low, uci = conf.high)
  }
  
  df_s %>%
    filter(cohort == spec$cohort) %>%
    expand_grid(cog = -1:1,
                term = names(fixef(mod))) %>%
    mutate(beta = ifelse(str_detect(term, !!reg_ex), 1, 0),
           beta = ifelse(str_detect(term, "age_ns1"), beta * age_ns1, beta),
           beta = ifelse(str_detect(term, "age_ns2"), beta * age_ns2, beta),
           beta = ifelse(str_detect(term, "cog"), beta * cog, beta)) %>%
    select(age, cog, term, beta) %>%
    group_by(age, cog) %>%
    filter(max(abs(beta)) != 0) %>%
    ungroup() %>%
    nest(spec = c(term, beta)) %>%
    mutate(spec = map(spec, ~ deframe(.x) %>% t()),
           res = map(spec, get_glht)) %>%
    unnest(res) %>%
    select(-spec)
}

bind_rows(margin = get_glht("margin"),
          prediction = get_glht("prediction"),
          .id = "pred_type")






mod <- lmer(bmi ~ splines::ns(age, 2) + cog*splines::ns(age, 2) + (1 | id), data)

margins(mod, variables = "cog", at = list(age = unique(data$age)))

get_ggeffect <- function(data, mod_form){
  mod <- lmer(mod_form, data)
  
  ggeffect(mod, c("cog [-1:1 by=1]", "age_s"), data = data)
}

mod <- df %>%
  filter(imp <= 10) %>%
  group_split(imp) %>%
  map(get_ggeffect, mod_form)

pool_scalar <- function(estimate, se){
  res <- pool.scalar(estimate, se^2)
  
  tibble(beta = res$qbar, se = sqrt(res$t)) %>%
    mutate(lci = qnorm(.025, beta, se),
           uci = qnorm(.975, beta, se))
}

map_dfr(mod, as_tibble, .id = "imp") %>%
  mutate(group = as.character(group) %>% as.numeric(),
         imp = as.numeric(imp)) %>%
  arrange(group, x, imp) %>%
  group_by(x, group) %>%
  summarise(pool_scalar(predicted, std.error),
            .groups = "drop") %>%
  filter(group == 0)

pool_predictions(mod) %>%
  as_tibble() %>%
  mutate(group = as.character(group) %>% as.numeric()) %>%
  filter(group == 0)



library(ggeffects)

ggeffect(mod, c("cog [-1:1 by=1]", "age")) %>%
  x
pool_predictions(mod) %>%
  as_tibble() %>%
  mutate(group = as.character(group) %>% as.numeric(),
         x = as.factor(x)) %>%
  ggplot() +
  aes(x = group, y = predicted, ymin = conf.low, ymax = conf.high, 
      color = x, fill = x) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_line() 

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


























# OLD -------
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

df_s <- df_bmi %>% 
  arrange(cohort, age) %>%
  group_by(cohort) %>%
  distinct(age) %>%
  mutate(ns(age, 2) %>%
           as_tibble() %>%
           mutate(across(everything(), as.double)) %>%
           rename_with(~ glue("age_ns{.x}"))) %>%
  ungroup()

df_bmi <- df_bmi %>%
  left_join(df_s, by = c("cohort", "age"))


# 2. Run Models ----
cohort <- "1958c"
cog_score <- "cog_rank"

df <- df_cog %>%
  filter(str_detect(cog_var, "g_closer")) %>%
  select(cohort, id, cog = all_of(cog_score)) %>%
  left_join(df_bmi, by = c("cohort", "id")) %>%
  left_join(df_cov, by = c("cohort", "id")) %>%
  filter(cohort == !!cohort) %>%
  add_count(id) %>%
  filter(n >= 3) %>%
  select(-n)

mod <- lmer(bmi ~ age_ns1 + age_ns2 + childhood_bmi +
              cog + I(cog*age_ns1) + I(cog*age_ns2) + 
              (1 + age_ns1 + age_ns2 | id), df)

cog_df <- tibble(cog_var = c("cog_std", "cog_rank"),
                 cog = list(c(-1:1), 0:2/2)) %>%
  unnest(cog)

df_s %>%
  expand_grid(filter(cog_df, cog_var == "cog_std")) %>%
  filter(cohort == !!cohort) %>%
  select(-cog_var) %>%
  expand_grid(term = names(fixef(mod))) %>%
  mutate(beta = ifelse(str_detect(term, "(Intercept|cog|age_ns)"), 1, 0),
         beta = ifelse(str_detect(term, "age_ns1"), beta * age_ns1, beta),
         beta = ifelse(str_detect(term, "age_ns2"), beta * age_ns2, beta),
         beta = ifelse(str_detect(term, "cog"), beta * cog, beta)) %>%
  select(cohort, age, cog, term, beta) %>%
  nest(spec = -c(cohort:cog)) %>%
  mutate(spec = map(spec,
                    ~ .x %>%
                      pivot_wider(names_from = term, values_from = beta) %>%
                      as.matrix()),
         res = map(spec, ~ multcomp::glht(mod, .x) %>%
                     confint() %>% tidy()))  %>%
  unnest(res) %>%
  ggplot() +
  aes(x = age, y = estimate, ymin = conf.low, ymax = conf.high, 
      color =  cog, fill = cog, group = cog) +
  geom_ribbon(color = NA, alpha = 0.2) +
  geom_line()




summary(mod)

df_s %>%
  filter(cohort == !!cohort) %>%
  expand_grid()

enframe(fixef(mod)) %>%
  mutate(value = ifelse(row_number() == 3, 1, 0)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  as.matrix() %>%
  multcomp::glht(mod, .)




poly(1:100, degree = 2) %>% 
  as_tibble() %>%
  rename_with(~ paste0("x", .x)) %>%
  mutate(x = row_number()) %>%
  ggplot() +
  aes(x = x, y = x2) +
  geom_line()


df_s %>%
  ggplot() +
  aes(x = age, y = age_ns2) +
  geom_line()

# 2. Run