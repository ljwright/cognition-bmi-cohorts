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
library(estimatr)

rm(list = ls())


# 1. Load Data ----
load("Data/main_results.Rdata")
load("Data/mnar_results.Rdata")
load("Data/attrition_results.Rdata")
load("Data/quantile_results.Rdata")


# 2. Main ----
# Text
res_string <- res_main %>%
  filter(cog_score == "ridit",
         height_score == "chart",
         str_length(cohort) == 5,
         sex == "all", type == "mi") %>%
  select(cohort, test_clean, mod_clean, string)

res_string %>%
  filter(test_clean == "Verbal @ age 11", mod_clean == "Sex-Adjusted") %>%
  mutate(string = str_replace(string, "\\(", "SD (95% CI = "),
         string = glue("{string} in the {cohort}")) %>%
  pull(string) %>%
  glue_collapse(", ")


# Plots
plot_main <- function(type){
  p <- res_main %>%
    filter(cog_score == "ridit",
           height_score == "chart",
           str_length(cohort) == 5,
           sex == "all") %>%
    filter(type == !!type) %>%
    mutate(cohort_clean = fct_rev(cohort_clean)) %>%
    ggplot() +
    aes(x = cohort_clean, y = beta, ymin = lci, 
        ymax = uci, label = string) +
    facet_grid(test_clean ~ mod_clean, switch = "y", scales = "free_y") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_pointrange(position = position_dodge(0.5)) +
    geom_text(aes(as.numeric(cohort_clean) + 0.2),
              size = 4, y = Inf, vjust = 0, hjust = 1.1) +
    scale_y_continuous(expand = c(0.05, 0.5)) +
    coord_flip() +
    theme_bw() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          strip.background.y = element_blank(),
          axis.title = element_text(size = rel(0.8)),
          legend.position = "none") +
    labs(x = NULL, y = "Marginal Effect")
  
  glue("Images/main_{type}.png") %>%
    ggsave(plot = p, width = 29.7, height = 21, units = "cm")
  
  return(p)
}

map(c("cc", "mi"), plot_main)

plot_sex <- function(type){
  p <- res_main %>%
    filter(cog_score == "ridit",
           height_score == "chart",
           str_length(cohort) == 5,
           sex != "all") %>%
    filter(type == !!type) %>%
    mutate(cohort_clean = fct_rev(cohort_clean),
           sex_clean = str_to_title(sex)) %>%
    ggplot() +
    aes(x = cohort_clean, y = beta, ymin = lci, 
        ymax = uci, label = string, color = sex_clean,
        shape = sex_clean) +
    facet_grid(test_clean ~ mod_clean, switch = "y", scales = "free_y") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_pointrange(position = position_dodge(0.6)) +
    scale_color_brewer(palette = "Dark2") +
    coord_flip() +
    theme_bw() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          strip.background.y = element_blank(),
          axis.title = element_text(size = rel(0.8)),
          legend.position = "bottom") +
    labs(x = NULL, y = "Marginal Effect", 
         color = NULL, shape = NULL)
  
  glue("Images/sex_{type}.png") %>%
    ggsave(plot = p, width = 29.7, height = 21, units = "cm")
  
  return(p)
}

map(c("cc", "mi"), plot_sex)


# 3. MNAR ----
plot_mnar <- function(test){
  res_g <- res_heat %>%
    filter(cohort_x < cohort_y) %>%
    filter(test == !!test,
           str_length(cohort_x) == 5,
           str_length(cohort_y) == 5)
  
  limits <- max(abs(res_g$beta)) * c(-1, 1)
  
  p <- ggplot(res_g) +
    aes(x = diff_x, y = diff_y, fill = beta) + # alpha = 1-p) +
    facet_grid(cohort_y ~ cohort_x, switch = "both") +
    geom_tile() +
    geom_abline(linetype = "dashed", color = "grey60") +
    theme_bw() +
    scale_fill_distiller(type = "div", limit = limits) +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          legend.position = "bottom",
          strip.background = element_rect(fill = "white", colour = "white"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    guides(fill = guide_colorbar(title.position = 'top', 
                                 title.hjust = .5,                                
                                 barwidth = unit(20, 'lines'), 
                                 barheight = unit(.5, 'lines'))) +
    labs(x = NULL, y = NULL, fill = "Effect Size Difference")
  
  glue("Images/mnar_{test}.png") %>%
    ggsave(plot = p, width = 29.7, height = 21, dpi = 600, units = "cm")
  
  return(p)
}

unique(res_heat$test) %>%
  set_names(., .) %>%
  map(plot_mnar)


# 4. Attrition ----
ggplot(res_att) + 
  aes(x = fup_miss, y = beta, ymin = lci, ymax = uci, color = sex) +
  facet_wrap(~ cohort_clean, scales = "free_x") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(position = position_dodge(0.5)) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  labs(x = "Attrition By Age", color = NULL,
       y = "Difference in\nHeight Chart\nZ-Score")
ggsave("Images/attrition_height.png", width = 21,
       height = 16, dpi = 600, units = "cm")


# 5. Quantile Regression ----
res_qreg <- res_qreg %>%
  mutate(tau_clean = glue("{tau*100}th") %>% factor(),
         cohort_f = cohort)

res_qreg %>%
  filter(type == "mi", cohort == "1958c", test == "verbal",
         mod == "basic", tau == c(1, 5, 9)/10) %>%
  mutate(across(beta:uci, round, 2),
         string = glue("{beta} SD ({lci}, {uci}) at the {tau*100}th percentile")) %>%
  pull(string)

qreg_string <- res_qreg %>%
  mutate(across(beta:uci, round, 2),
         string = glue("{beta} SD (95% CI = {lci}, {uci}) at the {tau*100}th percentile")) %>%
  arrange(tau) %>%
  filter(tau %in% c(0.1, 0.9)) %>%
  group_by(mod, test_clean, type, cohort) %>%
  summarise(string = glue_collapse(string, " and "),
            .groups = "drop_last") %>% 
  mutate(string = glue("{string} in the {cohort}")) %>%
  summarise(string = glue_collapse(string, ", "),
            .groups = "drop")

qreg_flx <- flextable(qreg_string) %>%
  valign(j = 1:3, valign = "center") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all")

save_as_docx(qreg_flx, path = "Tables/qreg_string.docx")


# Plot
plot_qreg <- function(mod, type){
  df <- res_qreg %>%
    filter(mod == !!mod, type == !!type)
  
  p <- ggplot(df) +
    aes(x = tau_clean, y = beta, ymin = lci, ymax = uci, group = cohort) +
    facet_grid(test_clean ~ cohort, switch = "y", scales = "free") +
    geom_hline(yintercept = 0) +
    geom_ribbon(alpha = 0.2, color = NA) +
    geom_line() +
    geom_line(aes(group = cohort_f),
              data = select(df, -cohort),
              linetype = "dashed") +
    guides(color = "none", fill = "none") +
    labs(x = "Percentile", y = NULL) +
    theme_bw() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          strip.background.y = element_blank(),
          axis.title = element_text(size = rel(0.8)),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none")
  
  glue("Images/qreg_{mod}_{type}.png") %>%
    ggsave(plot = p, width = 29.7, height = 21, dpi = 600, units = "cm")
  
  return(p)
}

res_qreg %>%
  distinct(mod, type) %$%
  map2(mod, type, plot_qreg)


# 6. Sensitivity Analysis ----
# White-Only or Measured Height
plot_sens <- function(df, stub){
  p <- df %>%
    filter(cog_score == "ridit",
           height_score == "chart",
           sex == "all") %>%
    filter(type == "mi") %>%
    mutate(cohort_clean = fct_rev(cohort_clean)) %>%
    ggplot() +
    aes(x = cohort_clean, y = beta, ymin = lci, 
        ymax = uci, label = string) +
    facet_grid(test_clean ~ mod_clean, switch = "y", scales = "free_y") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_pointrange(position = position_dodge(0.5)) +
    geom_text(aes(as.numeric(cohort_clean) + 0.2),
              size = 4, y = Inf, vjust = 0, hjust = 1.1) +
    scale_y_continuous(expand = c(0.05, 0.5)) +
    coord_flip() +
    theme_bw() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          strip.background.y = element_blank(),
          axis.title = element_text(size = rel(0.8)),
          legend.position = "none") +
    labs(x = NULL, y = "Marginal Effect")
  
  glue("Images/sens_{stub}.png") %>%
    ggsave(plot = p, width = 29.7, height = 21, units = "cm")
  
  return(p)
}

res_main %>%
  filter(!(cohort %in% c("2001c", "1970c_measure"))) %>%
  plot_sens("2001c_white")

res_main %>%
  filter(!(cohort %in% c("1970c", "2001c_white"))) %>%
  plot_sens("1970c_measure")

# Different Cognitive/Height Measures
plot_measure <- function(type, cog_score){
  p <- res_main %>%
    filter(sex == "all", str_length(cohort) == 5) %>%
    filter(type == !!type, cog_score == !!cog_score) %>%
    mutate(cohort_clean = fct_rev(cohort_clean) %>% ordered(),
           mod_clean = fct_rev(mod_clean)) %>%
    ggplot() +
    aes(x = mod_clean, y = beta, ymin = lci, 
        ymax = uci, label = string, color = cohort_clean) +
    facet_grid(test_clean ~ height_clean, switch = "y", scales = "free_x") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_pointrange(position = position_dodge(0.8)) +
    coord_flip() +
    theme_bw() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          strip.background.y = element_blank(),
          axis.title = element_text(size = rel(0.8)),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "bottom") +
    labs(x = NULL, y = "Marginal Effect", color = NULL)
  
  glue("Images/sens_{cog_score}.png") %>%
    ggsave(plot = p, width = 29.7, height = 21, units = "cm")
  
  return(p)
}

res_main %>%
  distinct(type, cog_score) %$%
  map2(type, cog_score, plot_measure)

