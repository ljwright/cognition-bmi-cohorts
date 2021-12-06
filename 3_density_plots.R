library(tidyverse)
library(summarytools)
library(glue)
library(ggpp)

rm(list = ls())

# 1. Load Data ----
load("Data/df_cog.Rdata")

# 2. Distribution of Cognition Measures ----
cog_dict <- c(g_all = "g (All)", g_closer = "g (CLOSER)",
              g_bas = "g (BAS)", maths = "Maths", verbal = "Verbal",
              vocab = "Vocabulary")


df_desc <- map_dfr(df_cog,
                   ~ .x %>%
                     select(id, matches("_(raw|resid)$"),
                            -matches("^(bmi|self)")) %>%
                     pivot_longer(-id, names_to = "cog_var", values_to = "cog_score"),
                   .id = "cohort") %>%
  drop_na() %>%
  mutate(type = ifelse(str_detect(cog_var, "_raw$"), "raw", "resid"),
         cog_var = str_replace(cog_var, "_(raw|resid)", ""),
         age = str_sub(cog_var, -2) %>% as.numeric(),
         age_cat = ifelse(age <= 11, "10/11", "15/16"),
         cog_stub = str_sub(cog_var, 1, -4),
         cog_clean = glue("{cog_dict[cog_stub]} @ Age {age_cat}")) %>%
  filter(!str_detect(cog_var, "_bas_")) %>%
  select(cohort, cog_clean, type, cog_score)

df_stat <- df_desc %>%
  group_by(cohort, cog_clean, type) %>%
  descr(cog_score) %>%
  tb() %>%
  select(cohort, cog_clean, type,
         n = n.valid, mean, sd) %>%
         # median = med, min, max, iqr) %>%
  pivot_longer(n:sd, names_to = "stat", values_to = "value") %>%
  mutate(stat = ifelse(stat %in% c("sd", "iqr"),
                       str_to_upper(stat),
                       str_to_title(stat)),
         value = round(value, 2),
         value = ifelse(stat == "N",
                        format(value, big.mark = ",") %>%
                          trimws() %>%
                          str_replace("\\.00", ""),
                        as.character(round(value, 2)))) %>%
  rename(` ` = stat, `  ` =  value) %>%
  nest(tb = 4:5)


# 3. Plot Figures ----
plot_cog <- function(type){
  df_desc %>%
    filter(type == !!type) %>%
    ggplot() +
    aes(x = cog_score) +
    facet_wrap(cohort ~ cog_clean, ncol = 6, scales = "free") +
    geom_density() +
    geom_table(data = filter(df_stat, type == !!type),
               aes(x = Inf, y = Inf, label = tb),
               table.theme = ttheme_gtminimal, hjust = 1.1, vjust = 1.1) +
    theme_bw() +
    labs(x = NULL, y = NULL)
}

plot_cog("raw")
plot_cog("resid")

