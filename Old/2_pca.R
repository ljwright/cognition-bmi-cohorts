library(tidyverse)
library(haven)
library(glue)
library(psych)
library(mice)
library(broom)
library(corrr)
library(splines)

rm(list = ls())


# 1. Load Data ----
cohorts <- c("1946c", "1958c", "1970c") %>%
  set_names()

cog_fld <- "D:/Projects/Cognition Measures/Data"

df_raw <- map(cohorts, 
              ~ glue("{cog_fld}/{.x}_cognition.dta") %>%
                read_dta(col_select = c(id, matches("_1."))))

df_bmi <- map(cohorts,
              ~ glue("Data/{.x}_cleaned.dta") %>%
                read_dta() %>%
                as_factor())



age <- "age_11"
cog <- "verbal_11"
cohort <- "1946c"

get_n <- function(age, cog, cohort){
  df <- df_raw[[cohort]] %>%
    select(cog = all_of(!!cog), age = all_of(!!age))
  
  n <- df %>%
    filter(is.na(age)) %>%
    drop_na(cog) %>%
    nrow()
  
  n*100/nrow(df)
}

get_r2 <- function(age, cog, cohort){
  
  mod_form <- "cog ~ ns(age, 2)"
  if (cohort == "1970c" & age == "age_16"){
    mod_form <- "cog ~ ns(age, 2) + home_test_16"
  }
  mod_form <- as.formula(mod_form)
  
  df <- df_raw[[cohort]] %>%
    rename(cog = all_of(!!cog), age = all_of(!!age))
  
  lm(cog ~ ns(age, 2), df) %>%
      glance() %>%
      pull(1)
}

tibble(cohort = cohorts) %>%
  mutate(age = map(cohort, ~ str_subset(names(df_raw[[.x]]), "^age")),
         cog = map(cohort, ~ str_subset(names(df_raw[[.x]]), "_..$"))) %>%
  unnest(age) %>%
  unnest(cog) %>%
  filter(!str_detect(cog, "^(age|home)"), 
         str_sub(age, -2) == str_sub(cog, -2)) %>%
  mutate(r2 = pmap_dbl(list(age, cog, cohort), get_r2),
         n = pmap_dbl(list(age, cog, cohort), get_n)) %>%
  arrange(desc(n))


# 2. Prepare Data for PCA ----
df_pca <- list()

get_resid <- function(cog, age){
  df <- tibble(cog = cog, age = age)
  
  lm(cog ~ age, df) %>%
    augment(newdata = df) %>%
    pull(.resid)
}

clean_cog <- function(df, cohort){
  df %>%
    semi_join(df_bmi[[cohort]], by = "id") %>%
    select(-matches("_1.$")) %>%
    pivot_longer(-id) %>%
    mutate(type = case_when(str_detect(name, "_raw$") ~ "raw",
                            str_detect(name, "_resid$") ~ "resid",
                            TRUE ~ "std"),
           name = str_replace(name, glue("_{type}"), "")) %>%
    pivot_wider(names_from = name, values_from = value)
}

# 1946c
df_pca$`1946c` <- df_raw$`1946c` %>%
  select(id,
         matches("_11"), -general_11, -matches("^iq"),
         age_15, reading_15, arithmetic_15) %>%
  mutate(across(matches("_11") & !matches("age_11"),
                list(raw = ~ .x,
                     resid = ~ get_resid(.x, age_11))),
         across(matches("_15") & !matches("age_15"),
                list(raw = ~ .x,
                     resid = ~ get_resid(.x, age_15)))) %>%
  clean_cog("1946c")

# 1958c
df_pca$`1958c` <- df_raw$`1958c` %>%
  select(id,
         matches("_11"), -general_11, -matches("^iq"),
         age_16, comprehension_16, maths_16) %>%
  mutate(across(matches("_11") & !matches("age_11"),
                list(raw = ~ .x,
                     resid = ~ get_resid(.x, age_11))),
         across(matches("_16") & !matches("age_16"),
                list(raw = ~ .x,
                     resid = ~ get_resid(.x, age_16)))) %>%
  clean_cog("1958c")

# 1970c
get_resid_home <- function(cog, age, home_test){
  df <- tibble(cog = cog, age = age, home_test = home_test)
  
  lm(cog ~ age + home_test, df) %>%
    augment(newdata = df) %>%
    pull(.resid)
}

df_pca$`1970c` <- df_raw$`1970c` %>%
  select(id, matches("_10"),
         age_16, home_test_16,
         apu_arithmetic_16, apu_vocab_a_16) %>%
  mutate(maths_10_raw = maths_10,
         maths_10_resid = get_resid(maths_10, age_fmt_10),
         across(matches("_10$") & !matches("^(age|maths)_"),
                list(raw = ~ .x,
                     resid = ~ get_resid(.x, age_bas_10))),
         across(matches("_16") & !matches("(age|home_test)_16"),
                list(raw = ~ .x,
                     resid = ~ get_resid_home(.x, age_16, home_test_16)))) %>%
  clean_cog("1970c")


# 3. Run PCA ----
pca_forms <- list(
  `1946c` = list(g_all_11 = str_subset(names(df_pca$`1946c`), "_11"),
                 g_closer_11 = c("verbal_11", "nonverbal_11",
                                 "arithmetic_11", "reading_11")),
  `1958c` = list(g_all_11 = str_subset(names(df_pca$`1958c`), "_11"),
                 g_closer_11 = c("verbal_11", "nonverbal_11",
                                 "arithmetic_11", "comprehension_11")),
  `1970c` = list(g_all_10 = str_subset(names(df_pca$`1970c`), "_10") %>%
                   str_subset("total", TRUE),
                 g_bas_10 = str_subset(names(df_pca$`1970c`), "^bas.*_10$"),
                 g_closer_10 = c("comprehension_total_10", "maths_10", 
                                 "bas_matrices_10", "bas_sim_10"))
)

get_pca <- function(cohort, type, name){
  new_name <- glue("{name}_{type}")
  
  df <- df_pca[[cohort]] %>%
    filter(type == !!type) %>%
    select(all_of(pca_forms[[cohort]][[name]]))
  
  pca_mod <- pca(df)
  
  lds <- loadings(pca_mod) %>% 
    enframe(name = "var", value = "loading") %>%
    mutate(loading = as.double(loading))
  
  prd <- predict(pca_mod, df) %>%
    as_tibble() %>%
    mutate(PC1 = wtd_scale(PC1)) %>%
    rename(!!new_name := 1)
  
  tibble(loadings = list(lds), pred = list(prd))
}

pca_res <- expand_grid(cohort = names(pca_forms),
                       type = c("raw", "resid")) %>%
  mutate(name = map(cohort, ~ names(pca_forms[[.x]]))) %>%
  unnest(name) %>%
  arrange(cohort, type) %>%
  mutate(res = pmap(list(cohort, type, name), get_pca)) %>%
  unnest(res)

pca_res %>%
  unnest(loadings) %>%
  group_by(cohort, type, name) %>%
  summarise(eigenvalue = sum(loading^2),
            vars = n(),
            .groups = "drop") %>%
  mutate(prop = eigenvalue/vars)


# 4. Final Data ----
clean_long <- function(df, cohort){
  df_g <- pca_res %>% 
    filter(cohort == !!cohort) %>%
    pull(pred) %>%
    reduce(bind_cols)
  
  df %>%
    pivot_wider(names_from = type, 
                values_from = matches("_..$")) %>%
    bind_cols(df_g)
}

df_cog <- list()

df_cog$`1946c` <- df_pca$`1946c` %>%
  select(id, type, 
         maths_11 = arithmetic_11, verbal_11,
         maths_15 = arithmetic_15, vocab_15 = reading_15) %>%
  clean_long("1946c")

df_cog$`1958c` <- df_pca$`1958c` %>%
  select(id, type, 
         maths_11 = arithmetic_11, verbal_11,
         maths_16, vocab_16 = comprehension_16) %>%
  clean_long("1958c")

df_cog$`1970c` <- df_pca$`1970c` %>%
  select(id, type, 
         maths_10, verbal_10 = bas_sim_10,
         maths_16 = apu_arithmetic_16, vocab_16 = apu_vocab_a_16) %>%
  clean_long("1970c")

# Check Missing Dates
df_cog$`1946c` %>%
  filter(is.na(g_all_11_resid), !is.na(g_all_11_raw)) %>%
  nrow()

df_cog$`1958c` %>%
  filter(is.na(g_all_11_resid), !is.na(g_all_11_raw)) %>%
  nrow()

df_cog$`1970c` %>%
  filter(is.na(g_all_10_resid), !is.na(g_all_10_raw)) %>%
  nrow()

df_cog$`1946c` %>%
  filter(is.na(vocab_15_resid), !is.na(vocab_15_raw)) %>%
  nrow()

df_cog$`1958c` %>%
  filter(is.na(vocab_16_resid), !is.na(vocab_16_raw)) %>%
  nrow()

df_cog$`1970c` %>%
  filter(is.na(vocab_16_resid), !is.na(vocab_16_raw)) %>%
  nrow()

df_cog$`1970c` %>%
  filter(is.na(vocab_16_resid), !is.na(vocab_16_raw)) %>%
  nrow()
