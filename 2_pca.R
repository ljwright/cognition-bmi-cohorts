library(tidyverse)
library(haven)
library(glue)
library(mice)
library(broom)
library(splines)
library(psych)
library(gallimaufr)
library(tictoc)

rm(list = ls())


# 1. Load Data ----
cohorts <- c("1946c", "1958c", "1970c") %>%
  set_names()

df_raw <- map(cohorts,
              ~ glue("Data/{.x}_cleaned.dta") %>%
                read_dta() %>%
                as_factor() %>%
                zap_formats() %>%
                zap_label())

# Drop Missing Sex
map_dbl(df_raw, ~ sum(is.na(.x$male)))

df_raw <- map(df_raw, ~ drop_na(.x, male))


# 2. Imputed Missing Age Data ----
impute_age <- function(cohort, age){
  df <- df_raw[[cohort]] %>%
    select(id, matches(age), 
           -matches("(bmi|self_report)"))
  
  no_data <- df %>% 
    select(-matches("^age")) %>%
    pivot_longer(-id) %>%
    group_by(id) %>%
    summarise(n = sum(!is.na(value)),
              .groups = "drop") %>%
    filter(n == 0) %>%
    pull(id)
  
  pred <- make.predictorMatrix(df)
  pred[, "id"] <- 0
  
  mice(df, m = 1, maxit = 10, predictorMatrix = pred,
       printFlag = FALSE, defaultMethod = "rf") %>% #pmm
    mice::complete(action = "long") %>%
    as_tibble() %>%
    select(id, matches("^age")) %>%
    mutate(across(matches("age"),
                  ~ ifelse(id %in% no_data, NA, .x)))
}

imp_age <- expand_grid(cohort = cohorts,
                       age = c("_(10|11)", "_(15|16)")) %>%
  mutate(res = map2(cohort, age, impute_age))

join_age <- function(cohort){
  imp_age %>%
    filter(cohort == !!cohort) %>%
    pull(res) %>%
    reduce(full_join, by = "id") %>%
    full_join(df_raw[[cohort]] %>%
                select(-matches("^age")),
              by = "id")
}

df_age <- map(cohorts, join_age)


# 3. Prepare Data for PCA ----
get_resid <- function(cog, age){
  df <- tibble(cog = cog, age = age)
  
  lm(cog ~ ns(age, 2), df) %>%
    augment(newdata = df) %>%
    pull(.resid)
}

clean_cog <- function(df, cohort){
  df %>%
    select(-matches("_1.$")) %>%
    pivot_longer(-id) %>%
    mutate(type = case_when(str_detect(name, "_raw$") ~ "raw",
                            str_detect(name, "_resid$") ~ "resid",
                            TRUE ~ "std"),
           name = str_replace(name, glue("_{type}"), "")) %>%
    pivot_wider(names_from = name, values_from = value)
}

# 1946c
df_pca <- list()

df_pca$`1946c` <- df_age$`1946c` %>%
  select(id, matches("_1.$"), -matches("^(bmi|self)")) %>%
  mutate(across(matches("_11") & !matches("age_11"),
                list(raw = ~ .x,
                     resid = ~ get_resid(.x, age_11))),
         across(matches("_15") & !matches("age_15"),
                list(raw = ~ .x,
                     resid = ~ get_resid(.x, age_15)))) %>%
  clean_cog("1946c")

# 1958c
df_pca$`1958c` <- df_raw$`1958c` %>%
  select(id, matches("_1.$"), -matches("^(bmi|self)"))  %>%
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
  select(id, matches("_1.$"), -matches("^(bmi|self)"))  %>%
  mutate(maths_10_raw = maths_10,
         maths_10_resid = get_resid(maths_10, age_fmt_10),
         across(matches("_10$") & !matches("^(age|maths)_"),
                list(raw = ~ .x,
                     resid = ~ get_resid(.x, age_bas_10))),
         across(matches("_16") & !matches("(age|home_test)_16"),
                list(raw = ~ .x,
                     resid = ~ get_resid_home(.x, age_16, home_test_16)))) %>%
  clean_cog("1970c")


# 4. Run PCA and Extract g ----
pca_forms <- list(
  `1946c` = list(g_all_11 = str_subset(names(df_pca$`1946c`), "_11"),
                 g_closer_11 = c("verbal_11", "nonverbal_11",
                                 "maths_11", "reading_11")),
  `1958c` = list(g_all_11 = str_subset(names(df_pca$`1958c`), "_11"),
                 g_closer_11 = c("verbal_11", "nonverbal_11",
                                 "maths_11", "comprehension_11")),
  `1970c` = list(g_all_10 = str_subset(names(df_pca$`1970c`), "_10") %>%
                   str_subset("total", negate = TRUE),
                 g_bas_10 = str_subset(names(df_pca$`1970c`), "^bas.*_10$") %>%
                   c("verbal_10"),
                 g_closer_10 = c("comprehension_total_10", "maths_10", 
                                 "bas_matrices_10", "verbal_10"))
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


# 5. Combine Data ---
combine_dfs <- function(cohort){
  df_g <- pca_res %>% 
    filter(cohort == !!cohort) %>%
    pull(pred) %>%
    reduce(bind_cols)
  
  df_cog <- df_pca[[cohort]] %>%
    select(id, type,
           matches("^maths_1.$"),
           matches("^verbal_1(0|1)"),
           matches("^vocab_1(5|6)")) %>%
    pivot_wider(names_from = type, 
                values_from = matches("_..$")) %>%
    bind_cols(df_g)
  
  df_age[[cohort]] %>%
    select(-(matches("_1.$") & !matches("^(bmi|self)"))) %>%
    full_join(df_cog, by = "id")
}

df_cog <- map(cohorts, combine_dfs)

save(df_cog, cohorts, file = "Data/df_cog.Rdata")


# 6. Run MICE ----
load("Data/df_cog.Rdata")

df_mice <- map(df_cog,
               ~ .x %>%
                 rename(childhood_bmi = matches("bmi_(10|11)")) %>%
                 select(-matches("_raw"), -cohort,
                        -matches("(bmi|self_report)_(0|1).$")) %>%
                 rename_with(~ str_replace(.x, "_resid", "")))

rm(df_cog)

tic()
imp <- list()
for(cohort in names(df_mice)){
  df <- df_mice[[cohort]]
  
  pred <- make.predictorMatrix(df)
  pred[, c("id", "male")] <- 0
  
  temp <- list()
  for (male in c(0, 1)){
    df_m <- filter(df, male == !!male)
    
    temp[[male + 1]] <- parlmice(df_m, n.core = 8, n.imp.core = 8,
                                 predictorMatrix = pred,
                                 defaultMethod = c("rf", "logreg", "polyreg", "polr"))
  }
  
  imp[[cohort]] <- rbind(temp[[1]], temp[[2]])
}
toc()

rm(df, df_mice, temp, pred, male, cohort)

save(imp, file = "Data/mice.Rdata")
