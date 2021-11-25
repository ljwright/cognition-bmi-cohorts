library(tidyverse)
library(lcmm)
library(splines)
library(glue)
library(broom)

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

df_gmm <- df_gmm %>%
  left_join(id_lookup, by = c("cohort", "id")) %>%
  left_join(df_s, by = "age") %>%
  add_count(gmm_id) %>%
  filter(n > 1) %>%
  select(gmm_id, bmi, age_ns1, age_ns2) %>%
  structure(class = "data.frame")

rm(df_cog)


# 2. Run Models ----
m1 <- hlme(fixed = bmi ~ age_ns1 + age_ns2,
           random = ~ age_ns1 + age_ns2,
           subject = 'gmm_id',
           data = df_gmm)

cl <-  makeCluster(8)
mods <- vector(mode = "list", length = 6)
for (i in 1:6){
  groups <- i + 1
  clusterExport(cl=cl, "groups")
  
  mods[[i]] <- gridsearch(
    hlme(fixed = bmi ~ age_ns1 + age_ns2,
         mixture = ~ age_ns1 + age_ns2,
         random = ~ age_ns1 + age_ns2,
         subject = 'gmm_id', 
         ng = groups,
         nwg = TRUE,
         data = df_gmm),
    rep = 100, maxiter = 30, minit = m1, cl = cl
  )
}
Sys.time()
stopCluster(cl)

mod <- hlme(fixed = bmi ~ age_ns1 + age_ns2,
            mixture = ~ age_ns1 + age_ns2,
            random = ~ age_ns1 + age_ns2,
            subject = "gmm_id", 
            data = df_gmm,
            ng = 2, nwg = TRUE)


# 2. Model Fit ----
get_ord <- function(mod){
  as_tibble(mod$pprob) %>%
    pivot_longer(matches("prob")) %>%
    group_by(name) %>%
    summarise(n = sum(value), .groups = "drop") %>%
    mutate(class = str_sub(name, -1) %>% as.integer()) %>%
    mutate(class_prop = round(n*100/sum(n), 2)) %>%
    arrange(desc(class_prop)) %>%
    mutate(class_ord = glue("Class {row_number()}")) %>%
    select(class, class_ord, class_prop)
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
  
  mods %>%
    mutate(fit = map(mod, get_summary)) %>%
    select(-mod) %>%
    unnest(fit) %>%
    pivot_longer(-c(groups, model)) %>%
    mutate(clean_name = clean_names[name]) %>%
    ggplot() +
    aes(x = groups, y = value) +
    facet_wrap(~ clean_name, scales = "free_y") +
    geom_line() +
    geom_point() +
    theme_minimal() +
    labs(x = "Classes", y = NULL)
}


# 3. Get Predictions ----
get_prediction <- function(mod){
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
    left_join(df_s %>%
                pivot_longer(-age, names_to = "term"),
              by = "term") %>%
    select(gmm_id, age, class_ord, beta, value) %>%
    group_by(gmm_id, age, class_ord) %>%
    summarise(pred = sum(beta*value),
              .groups = "drop")
  
  return(mod_param)
}

plot_mod <- function(mod){
  get_prediction(mod) %>%
    pluck("pred") %>%
    group_by(class) %>%
    filter(gmm_id %in% sample(unique(gmm_id), 1000)) %>%
    ungroup() %>%
    mutate(class = factor(class)) %>%
    ggplot() +
    aes(x = age, y = pred, color = class, group = gmm_id) +
    geom_line(alpha = 0.1) +
    coord_cartesian(ylim = c(15, 55)) +
    theme_minimal()
}
