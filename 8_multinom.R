library(tidyverse)
library(broom)
library(tictoc)
library(paletteer)
library(glue)
library(furrr)
library(here)
library(splines)
library(nnet)
library(Hmisc)
library(mc2d)
library(mice)
library(lcmm)
library(mice)
library(miceadds)
library(gallimaufr)
library(flextable)
library(officer)
library(beepr)

rm(list = ls())

# 1. Load Data ----
load("Data/df_cog.Rdata")
load("Data/df_lcmm.Rdata")
load("Data/lcmm_fit.Rdata")
load("Data/lgca_mods.Rdata")

mod <- lgca[[4]]

save(mod, file = "Data/lcmm_choice.Rdata")

rm(get_summary, plot_fit, df_s, df_gmm, lgca)

# 2. Probability Table ----
probs_tbl <- get_probs(mod) %>%
  select(class_ord, mean_prob,
         prop1, prop2, prop3, prop4) %>%
  mutate(mean_prob = glue("{round(mean_prob, 1)}%")) %>%
  flextable() %>%
  set_header_labels(class_ord = "", mean_prob = "",
                    prop1 = "Class 1", prop2 = "Class 2",
                    prop3 = "Class 3", prop4 = "Class 4") %>%
  add_header(class_ord = "Class", mean_prob = "Class Proportions",
             prop1 = "Average Class Probability",
             prop2 = "Average Class Probability", 
             prop3 = "Average Class Probability", 
             prop4 = "Average Class Probability") %>%
  border_remove() %>%
  merge_v(1) %>%
  merge_h(part = "header") %>%
  border_inner_h(border = fp_border(color="grey50", width = 1, style = "dashed"),
                 part = "body") %>%
  hline_top(border = fp_border(color="black", width = 2), part = "all") %>%
  hline_bottom(border = fp_border(color="black", width = 2), part = "all") %>%
  fix_border_issues(part = "all") %>%
  align(j = 1:6, align="center", part = "all") %>%
  valign(j = 1, valign = "center") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  autofit()
probs_tbl


# 3. Get Pseudo Classes ----
df_prob <- mod$pprob %>%
  as_tibble() %>%
  select(-class)

get_vmultinom <- function(x){
  df_prob %>%
    select(-gmm_id) %>%
    rmultinomial(nrow(.), 1, .) %>%
    as_tibble() %>%
    mutate(gmm_id = df_prob$gmm_id) %>%
    pivot_longer(-gmm_id, names_to = "pseudo") %>%
    filter(value == 1) %>%
    mutate(pseudo = str_sub(pseudo, 2) %>% as.integer()) %>%
    select(gmm_id, pseudo)
}

set.seed(1)
tic()
pseudo_draws <- map_dfr(1:60, get_vmultinom, .id = "imp") %>%
  left_join(get_ord(mod), by = c("pseudo" = "class")) %>%
  mutate(imp = as.numeric(imp),
         class_ord = factor(class_ord)) %>%
  left_join(id_lookup, by = "gmm_id") %>%
  select(imp, cohort, id, class_ord)
toc()

save(pseudo_draws, file = "Data/pseudo_draws.Rdata")

rm(df_prob, get_vmultinom)


# 3. Run Imputations ----
load("Data/pseudo_draws.Rdata")
load("Data/df_cog.Rdata")

covars <- c("male", "maternal_age", 
            "bmi_10", "bmi_11",
            "birth_weight", "puberty_girl",
            "puberty_boy", "mother_edu_level",
            "father_class", "survey_weight")

get_mice <- function(imp, cohort, male){
  df_mice <- df_cog[[cohort]] %>%
    filter(male == !!male) %>%
    select(id, g = matches("g_all_1._resid"),
           any_of(covars)) %>%
    rename(childhood_bmi = matches("bmi_(10|11)")) %>%
    left_join( # inner_join(
      pseudo_draws %>%
        filter(imp == !!imp,
               cohort == !!cohort) %>%
        select(-imp, -cohort), 
      by = "id"
    )
  
  preds <- make.predictorMatrix(df_mice)
  preds[, "id"] <- 0
  
  mice(df_mice, m = 1, maxit = 10, 
       seed = imp, print = FALSE) %>%
    complete() %>%
    as_tibble()
}

plan(multisession, workers = 4)
tic()
imp_long <- pseudo_draws %>%
  distinct(imp, cohort) %>%
  expand_grid(male = 0:1) %>%
  mutate(imputed = future_pmap(list(imp, cohort, male), get_mice, .progress = TRUE))
toc()
future:::ClusterRegistry("stop")

save(imp_long, file = "Data/imp_long.Rdata")

# 4. Run Models ----
mod_covars <- c("g", "class_ord", "survey_weight", # "childhood_bmi",
                "birth_weight", "mother_edu_level", "father_class")

imp_mod <- imp_long %>%
  mutate(imputed = map(imputed, 
                       ~ .x %>%
                         select(all_of(mod_covars)))) %>%
  unnest(imputed) %>%
  mutate(cohort = factor(cohort)) %>%
  group_by(imp, cohort, male) %>%
  mutate(g = wtd_scale(g, survey_weight)) %>%
  ungroup() %>%
  nest(data = -imp)


mod_forms <- list(basic = "male",
                  sep = c("male", "father_class", "mother_edu_level"),
                  all = c("male", "birth_weight", #"childhood_bmi",
                          "mother_edu_level", "father_class"))

get_mult <- function(imp, mod_name){
  data <- imp_mod$data[[imp]]
  
  mod_form <- glue_collapse(mod_forms[[mod_name]], " + ") %>%
    paste("class_ord ~ cohort*g +", .) %>%
    as.formula()
  
  multinom(mod_form, data, weights = survey_weight) %>%
    tidy() %>%
    select(y.level, term, beta = estimate, se = std.error)
}

get_pool <- function(terms, beta, se){
  pool_mi(beta, se = se) %>%
    summary() %>%
    as_tibble() %>%
    select(beta = 1, se = 2, lci = 5, uci = 6) %>%
    mutate(term = terms[[1]], .before = 1)
}

pool_mult <- function(mod_name){
  map_dfr(imp_mod$imp, get_mult, mod_name, .id = "imp") %>%
    arrange(y.level, imp, term) %>%
    chop(c(term, beta, se)) %>%
    group_by(y.level) %>%
    summarise(get_pool(term, beta, se),
              .groups = "drop")
}

tic()
df_mult <- tibble(mod_name = names(mod_forms)) %>%
  mutate(res = map(mod_name, pool_mult))
toc()
future:::ClusterRegistry("stop")

save(df_mult, file = here("Data", "df_mult.Rdata"))

x <- pool_mult("basic")

x %>%
  filter(str_sub(term, -1) == "g") %>%
  mutate(across(c(beta, lci, uci), exp)) %>%
  ggplot() +
  aes(x = term, y = beta, ymin = lci, ymax = uci) +
  facet_wrap(~ y.level) +
  geom_hline(yintercept = 1) +
  geom_pointrange() +
  scale_y_log10()


plot_lgca(mod) + 
  annotate(geom = "rect", ymin = 15, ymax = 20, xmin = -Inf, xmax = Inf, alpha = 0.1) +
  annotate(geom = "rect", ymin = 25, ymax = 30, xmin = -Inf, xmax = Inf, alpha = 0.1)
