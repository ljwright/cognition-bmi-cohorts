df_s <- df_bmi %>%
  mutate(age = str_sub(bmi_var, -2) %>% as.integer()) %>%
  distinct(cohort, age) %>%
  group_by(cohort) %>%
  mutate(# transform_age(age, poly, "age_p"),
    transform_age(age, ns, "age_ns"),
    age_ns0 = 1) %>%
  arrange(cohort, age) %>%
  ungroup()

gcm_specs <- mod_specs %>%
  select(-spec_id, -bmi_var, -bmi_score) %>%
  distict() %>%
  filter(str_detect(cog_var, "^g_all"),
         cog_score != "cog_ridit") %>%
  mutate(spec_id = row_number(),
         .before = 1)

prepare_df_gcm <- function(spec_id){
  spec <- slice(gcm_specs, !!spec_id)
  
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

get_form_gcm <- function(spec_id){
  spec <- slice(gcm_specs, !!spec_id)
    
    ind_vars <- glue_collapse(mods[[spec$mod]], " + ")
    if (spec$sex == "all" & spec$mod == "basic") ind_vars <- "male"
    if (spec$sex == "all" & spec$mod != "basic") ind_vars <- glue("{ind_vars} + male")
    if (length(ind_vars) != 0) ind_vars <- glue("+ {ind_vars}")
    
    paste("bmi ~ -1 + age_ns0 + age_ns1 + age_ns2 +
       cog + cog*age_ns1 + cog*age_ns2 +
       (-1 + age_ns0 + age_ns1 + age_ns2 | id)",
          ind_vars)
}

run_gcm<- function(spec_id){
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
      expand_grid(cog = c(-3:3, 0.5),
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
  
  res <- group_split(df, imp) %>%
    map_dfr(get_preds, .id = "imp")
  if (spec$type == "cc") res <- uncount(res, 2, .id = "imp")
  
  return(res)
}