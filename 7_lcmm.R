library(tidyverse)
library(lcmm)
library(splines)
library(glue)
library(broom)
library(patchwork)

rm(list = ls())

# 1. Load Data ----
load("Data/df_cog.Rdata")

df_gmm <- map_dfr(df_cog,
                  ~ .x %>%
                    select(id, matches("^bmi")) %>%
                    pivot_longer(-id, names_to = "age", values_to = "bmi"),
                  .id = "cohort") %>%
  mutate(age = str_sub(age, -2) %>% as.numeric()) %>%
  filter(age >= 15) %>%
  drop_na()

df_s <- df_gmm %>% 
  arrange(age) %>%
  distinct(age) %>%
  mutate(ns(age, 2) %>%
           as_tibble() %>%
           mutate(across(everything(), as.double)) %>%
           rename_with(~ glue("age_ns{.x}")),
         intercept = 1)

id_lookup <- df_gmm %>%
  group_by(cohort, id) %>%
  summarise(gmm_id = cur_group_id(),
            .groups = "drop")

df_gmm %>%
  left_join(id_lookup, by = c("cohort", "id")) %>%
  left_join(df_s, by = "age") %>%
  add_count(gmm_id, name = "fups")  %>% 
  count(cohort, fups) %>%
  mutate(n = n/fups) %>%
  add_count(cohort, wt = n, name = "total") %>%
  mutate(prop = 100*n/total) %>%
  group_by(cohort) %>%
  mutate(cum_prop = ifelse(row_number() == 1, 100, 100 - lag(cumsum(prop)))) %>% 
  ungroup() %>%
  select(cohort, fups, cum_prop) %>%
  pivot_wider(names_from = fups, values_from = cum_prop)

df_gmm <- df_gmm %>%
  left_join(id_lookup, by = c("cohort", "id")) %>%
  left_join(df_s, by = "age") %>%
  add_count(gmm_id) %>%
  filter(n > 3) %>%
  select(cohort, gmm_id, bmi, age_ns1, age_ns2) %>%
  structure(class = "data.frame")

rm(df_cog)

save(df_s, df_gmm, id_lookup,
     file = "Data/df_lcmm.Rdata")

# 2. Model Fit ----
get_ord <- function(mod){
  as_tibble(mod$pprob) %>%
    pivot_longer(matches("prob")) %>%
    group_by(name) %>%
    summarise(n = sum(value), .groups = "drop") %>%
    mutate(class = str_sub(name, -1) %>% as.integer()) %>%
    mutate(class_prop = round(n*100/sum(n), 2)) %>%
    arrange(desc(class_prop)) %>%
    mutate(class_ord = glue("Class {row_number()}"),
           class_long = glue("{class_ord} ({round(class_prop, 1)}%)")) %>%
    select(class, class_ord, class_prop, class_long)
}

get_probs <- function(mod){
  av_prob <- mod$pprob %>%
    pivot_longer(-c(gmm_id, class)) %>%
    group_by(name) %>%
    summarise(mean_prob = 100*mean(value),
              .groups = "drop") %>%
    mutate(name = str_replace(name, "prob", "") %>% as.numeric()) %>%
    rename(class = name)
  
  prob_ord <- get_ord(mod) %>%
    mutate(name = glue("prob{class}"),
           name_ord = glue("prop{str_sub(class_ord, -1)}")) %>%
    select(name, name_ord)
  
  prob_matrix <- mod$pprob %>%
    pivot_longer(-c(gmm_id, class)) %>%
    group_by(class, name) %>%
    summarise(mean = mean(value) %>% round(3), .groups = "drop") %>%
    left_join(prob_ord, by = "name") %>%
    arrange(class, name_ord) %>%
    select(-name) %>%
    pivot_wider(names_from = "name_ord", values_from = "mean")
  
  get_ord(mod) %>%
    left_join(av_prob, by = "class") %>%
    left_join(prob_matrix, by = "class") %>%
    select(-class, -class_long, -class_prop)
}

get_summary <- function(mod){
  entropy <- mod$pprob %>%
    pivot_longer(-c(gmm_id, class)) %>%
    mutate(value = value*log(value)) %>%
    summarise(value = sum(value), 
              N = length(unique(gmm_id)),
              G = mod$ng) %>%
    mutate(e = 1 + value/(N*log(G))) %>%
    pull(e)
  
  tibble(bic = mod$BIC, aic = mod$AIC, 
         entropy = entropy, log_lik = mod$loglik)
}

plot_fit <- function(mods){
  clean_names <- c(bic = "BIC", aic = "AIC", 
                   entropy = "Entropy",
                   log_lik = "Log-Likelihood")
  
  map_dfr(mods, get_summary, .id = "groups") %>%
    mutate(groups = as.numeric(groups)) %>%
    pivot_longer(-groups) %>%
    mutate(clean_name = clean_names[name]) %>%
    ggplot() +
    aes(x = groups, y = value) +
    facet_wrap(~ clean_name, scales = "free_y") +
    geom_line() +
    geom_point() +
    theme_minimal() +
    labs(x = "Classes", y = NULL)
}

save(plot_fit, get_summary, get_ord, get_probs,
     file = "Data/lcmm_fit.Rdata")

# 3. Run Models ----
# LGCA
lgca <- vector(mode = "list", length = 5)

lgca[[1]] <- hlme(fixed = bmi ~ age_ns1 + age_ns2,
                  subject = 'gmm_id',
                  data = df_gmm)

for (groups in 2:5){
  lgca[[groups]] <- hlme(fixed = bmi ~ age_ns1 + age_ns2,
                         mixture = ~ age_ns1 + age_ns2,
                         subject = "gmm_id", B = lgca[[1]],
                         data = df_gmm,
                         ng = groups)
}

save(lgca, file = "Data/lgca_mods.Rdata")

# GMM
gmm <- vector(mode = "list", length = 5)

gmm[[1]] <- hlme(fixed = bmi ~ age_ns1 + age_ns2,
                 random = ~ age_ns1 + age_ns2,
                 subject = 'gmm_id',
                 data = df_gmm)

for (groups in 2:5){
  gmm[[groups]] <- hlme(fixed = bmi ~ age_ns1 + age_ns2,
                         mixture = ~ age_ns1 + age_ns2,
                         random = ~ age_ns1 + age_ns2,
                         subject = "gmm_id", B = gmm[[1]], 
                         data = df_gmm,
                         ng = groups)
}

# Grid Search (GMM)
m1 <- hlme(fixed = bmi ~ age_ns1 + age_ns2,
           random = ~ age_ns1 + age_ns2,
           subject = 'gmm_id',
           data = df_gmm)

cl <-  makeCluster(8)
mods <- vector(mode = "list", length = 4)

for (groups in 2:4){
  groups <- i + 1
  clusterExport(cl=cl, "groups")
  
  mods[[i]] <- gridsearch(
    hlme(fixed = bmi ~ age_ns1 + age_ns2,
         mixture = ~ age_ns1 + age_ns2,
         random = ~ age_ns1 + age_ns2,
         subject = 'gmm_id', 
         ng = groups, nwg = TRUE,
         data = df_gmm),
    rep = 100, maxiter = 30, minit = m1, cl = cl
  )
}
Sys.time()
stopCluster(cl)


# 4. LGCA Predictions ----
get_lgca_pred <- function(mod){
  mod_param <- list()
  
  mod_param$fe <- enframe(coef(mod), name = "term", value = "fe") %>%
    filter(!str_detect(term, "^(cholesky|stderr)")) %>%
    separate(term, c("term", "class"), sep = " ") %>%
    mutate(class = str_sub(class, -1) %>% as.integer())
  
 mod_param$class <- mod$pprob %>%
    as_tibble() %>%
    select(gmm_id, class)
  
  mod_param$pred <- mod_param$fe %>%
    left_join(get_ord(mod), by = "class") %>%
    left_join(df_s %>%
                pivot_longer(-age, names_to = "term"),
              by = "term") %>%
    select(age, class_long, fe, value) %>%
    group_by(age, class_long) %>%
    summarise(pred = sum(fe*value),
              .groups = "drop")
  
  return(mod_param)
}

plot_lgca <- function(mod){
  get_lgca_pred(mod) %>%
    pluck("pred") %>%
    mutate(class_long = factor(class_long)) %>%
    ggplot() +
    aes(x = age, y = pred, color = class_long) +
    geom_line(alpha = 1) +
    coord_cartesian(ylim = c(15, 55)) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    labs(x = "Age", y = "Predicted BMI", color = NULL) +
    theme(legend.position = "right")
}

map(lgca[2:5], plot_lgca) %>% reduce(`+`)

# 5. GMM Predictions ----
get_gmm_pred <- function(mod){
  mod_param <- list()
  
  mod_param$fe <- enframe(coef(mod), name = "term", value = "fe") %>%
    filter(!str_detect(term, "^(cholesky|stderr)")) %>%
    separate(term, c("term", "class"), sep = " ") %>%
    mutate(class = str_sub(class, -1) %>% as.integer())
  
  mod_param$re <- mod$predRE %>%
    as_tibble() %>%
    pivot_longer(-gmm_id,
                 names_to = "term",
                 values_to = "re")
  
  mod_param$class <- mod$pprob %>%
    as_tibble() %>%
    select(gmm_id, class)
  
  mod_param$pred <- mod_param$fe %>%
    left_join(mod_param$re, by = "term") %>%
    right_join(mod_param$class, by = c("gmm_id", "class")) %>%
    left_join(get_ord(mod), by = "class") %>%
    mutate(beta = fe + re) %>%
    mutate(beta = fe) %>%
    left_join(df_s %>%
                pivot_longer(-age, names_to = "term"),
              by = "term") %>%
    select(gmm_id, age, class_ord, beta, value) %>%
    group_by(gmm_id, age, class_ord) %>%
    select(age, class_ord, beta, value) %>%
    group_by(age, class_ord) %>%
    summarise(pred = sum(beta*value),
              .groups = "drop")
  
  return(mod_param)
}

plot_gmm <- function(mod){
  get_prediction(mod) %>%
    pluck("pred") %>%
    group_by(class_ord) %>%
    filter(gmm_id %in% sample(unique(gmm_id), 1000)) %>%
    ungroup() %>%
    mutate(class_ord = factor(class_ord)) %>%
    ggplot() +
    aes(x = age, y = pred, color = class_ord, group = gmm_id) +
    geom_line(alpha = 1) + #0.1) +
    coord_cartesian(ylim = c(15, 55)) +
    theme_minimal()
}
