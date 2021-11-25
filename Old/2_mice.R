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
library(scales)
library(estimatr)
library(splines)

rm(list = ls())

# 1. Load Data ----
load_dta <- function(cohort){
  glue("Data/{cohort}_cleaned.dta") %>%
    read_dta() %>%
    as_factor() %>%
    zap_label() %>%
    zap_formats()
}

df_raw <- c("1946c", "1958c", "1970c", "2001c") %>%
  set_names(., .) %>%
  map(load_dta)

map_dbl(df_raw, ~ sum(is.na(.x$male))) # Missing Sex

df_raw <- map(df_raw,
              ~ .x %>%
                drop_na(male))

save(df_raw, file = "Data/df_raw.Rdata")
rm(load_dta)

# 2. Linear Regressions
get_bmi <- function(df){
  df %>%
    select(id, matches("^bmi")) %>%
    pivot_longer(-id, names_to = "age", values_to = "bmi") %>%
    mutate(age = str_sub(age, -2) %>% as.numeric()) %>%
    drop_na()
}

df_bmi <- map_dfr(df_raw, get_bmi, .id = "cohort")

get_cog <- function(df){
  df %>%
    select(id, matches("resid")) %>%
    pivot_longer(-id, names_to = "test", values_to = "cog_score") %>%
    drop_na() %>%
    group_by(test) %>%
    mutate(cog_rank = percent_rank(cog_score)) %>%
    ungroup()
}

df_cog <- map_dfr(df_raw, get_cog, .id = "cohort")

get_inv <- function(df){
  df %>%
    select(id, male, survey_weight, father_class)
}

df_inv <- map_dfr(df_raw, get_inv, .id = "cohort")


get_lm <- function(df){
  lm(bmi ~ cog_rank, df) %>%
    tidy(conf.int = TRUE) %>%
    filter(term == "cog_rank")
}

df_x <- df_bmi %>%
  left_join(df_cog, by = c("cohort", "id")) %>%
  drop_na() %>%
  # left_join(df_inv, by = c("cohort", "id")) %>%
  nest(data = -c(cohort, age, test)) %>%
  separate(test, into = c("test_type", "test_resid", "test_age"), sep = "_", convert = TRUE) %>%
  mutate(test_age = ifelse(test_age < 14, 11, 16),
         test = glue("{test_type} @ {test_age}"))

df_x %>%
  mutate(mod = map(data, get_lm)) %>%
  dplyr::select(-data, -test_resid) %>%
  unnest(mod) %>%
  ggplot() +
  aes(x = age, y = estimate, ymin = conf.low, ymax = conf.high, 
      color = cohort, fill = cohort) +
  facet_wrap(~ test) +
  geom_hline(yintercept = 0) +
  geom_ribbon(color = NA, alpha = 0.2) +
  geom_line()

# 3. GAMLSS ----
library(gamlss)
get_gamlss <- function(df){
  mod <- gamlss(bmi ~ cog_rank,
                sigma.formula = ~ cog_rank,
                nu.formula = ~ cog_rank,
                family = BCCGo, 
                data = df,
                trace = FALSE)
  
  coefs <- coefAll(mod)
  coefs_low <- map(coefs, 1)
  coefs_high <- map(coefs, sum)
  
  bind_rows(low = as_tibble(coefs_low),
            high = as_tibble(coefs_high),
            .id = "cog_score") %>%
    uncount(99, .id = "centile") %>%
    mutate(tau = centile/100,
           bmi = exp(mu) * (1 + nu * exp(sigma) * qnorm(tau))^(1/nu))
}

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


# 4. Splines ----
add_splines <- function(df){
  get_splines <- function(x){
    ns(x, 3) %>%
      as_tibble() %>%
      mutate(across(everything(), as.double)) %>%
      rename_with(~ glue("cog_ns{.x}"))
  }
  
  df %>%
    mutate(get_splines(cog_rank))
}

get_lm_splines <- function(df){
  df_spline <- add_splines(df)
  
  mod <- lm(bmi ~ cog_ns1 + cog_ns2 + cog_ns3, df_spline)
  
  get_glht <- function(string, mod){
    multcomp:::glht(mod, string) %>%
      confint() %>% 
      tidy() %>%
      dplyr::select(beta = 2, lci = 3, uci = 4)
  }
  
  df_spline %>%
    arrange(cog_rank) %>%
    distinct(across(cog_rank:cog_ns3)) %>%
    pivot_longer(-cog_rank) %>%
    filter(value != 0) %>%
    mutate(string = glue("{value}*{name}")) %>%
    group_by(cog_rank) %>%
    summarise(string = glue_collapse(string, " + ")) %>%
    mutate(string = glue("{string} = 0"),
           res = map(string, get_glht, mod)) %>%
    dplyr::select(-string) %>%
    unnest(res)
}

df_spline <- df_x %>%
  mutate(mod = map(data, get_lm_splines)) %>%
  dplyr::select(-data, -test_resid) %>%
  unnest(mod)

df_spline %>%
  arrange(cog_rank) %>%
  filter(!str_detect(test, "^comp"),
         age >= 10) %>%
  ggplot() +
  aes(x = cog_rank, y = beta, color = age, group = age) +
  facet_grid(test ~ cohort) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed")
