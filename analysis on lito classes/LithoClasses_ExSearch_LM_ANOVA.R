# Key findings:
# - Exhaustive variable selection identified `sqrt_slope`, `inv_vs30`, and the relevant lithology group dummy as the best predictors for each classification.
# - Diagnostic plots confirm acceptable model assumptions with no major violations.
# - Nested ANOVA revealed significant group effects for Litho_Group_2 and Litho_Group_3, but not for Litho_Group_1.
#-AIC,BIC and adjusted R^2 is very similar for groupping 2 and 3.
library(dplyr)
library(leaps)
library(fastDummies)
library(readr)
library(tibble)

itac <- read_csv("FINAL_DATASET.csv") %>%
  mutate(
    Lito_Bucci_simple = as.integer(Lito_Bucci_simple),
    Litho_Group_1 = case_when(
      Lito_Bucci_simple %in% c(6,7,8,1)       ~ "STIFF ROCKS",
      Lito_Bucci_simple %in% c(5)             ~ "SOFT ROCKS",
      Lito_Bucci_simple %in% c(2,3,4,9,10,11) ~ "SOILS",
      TRUE                                    ~ NA_character_
    ),
    Litho_Group_2 = case_when(
      Lito_Bucci_simple %in% c(1,6,7,8)       ~ "STIFF ROCKS",
      Lito_Bucci_simple %in% c(9,10,11,5)     ~ "SOFT ROCKS",
      Lito_Bucci_simple %in% c(2,3,4)         ~ "SOILS",
      TRUE                                    ~ NA_character_
    ),
    Litho_Group_3 = case_when(
      Lito_Bucci_simple %in% c(7,6,1)         ~ "STIFF ROCKS",
      Lito_Bucci_simple %in% c(8,4,2,10,3,9)  ~ "SOFT ROCKS",
      Lito_Bucci_simple %in% c(11,5)          ~ "SOILS",
      TRUE                                    ~ NA_character_
    ),
    # slope & vs30 transforms 
    log_slope   = log(slope + 1e-5),
    sqrt_slope  = sqrt(slope),
    inv_slope   = 1 / (slope + 1e-5),
    sqrt_vs30   = sqrt(Vs30_Brunelli_simple),
    inv_vs30    = 1 / (Vs30_Brunelli_simple + 1e-5)
  )

#Exhaustive search 
run_litho_exh <- function(df, group_col) {
  dfp <- df %>%
    filter(!is.na(.data[[group_col]]), !is.na(dS2S_SA_T10_000)) %>%
    dummy_cols(
      select_columns         = group_col,
      remove_first_dummy     = TRUE,    # k???1 dummy
      remove_selected_columns = TRUE
    ) %>%
    rename_with(make.names)
  
  slope_vars <- c("sqrt_slope", "inv_slope")
  vs30_vars  <- c("sqrt_vs30",  "inv_vs30")
  dummy_vars <- grep(paste0("^", make.names(group_col), "_"), names(dfp), value = TRUE)
  
  best_adj  <- -Inf
  best_vars <- NULL
  
  for(s in slope_vars) for(v in vs30_vars) {
    vars <- c(s, v, dummy_vars)
    sub  <- dfp %>% select(all_of(vars), dS2S_SA_T10_000) %>% na.omit()
    if(nrow(sub) > 50) {
      ex  <- regsubsets(x = sub[, vars], y = sub$dS2S_SA_T10_000,
                        nbest = 1, nvmax = length(vars), method = "exhaustive")
      sm  <- summary(ex)
      bi  <- which.max(sm$adjr2)
      sel <- names(which(sm$which[bi,]))[-1]
      if(sm$adjr2[bi] > best_adj) {
        best_adj  <- sm$adjr2[bi]
        best_vars <- sel
      }
    }
  }
  
  # Print selected predictors and model summary
  cat(sprintf("\n???? [%s] Best predictors (Adj R?? = %.4f):\n  %s\n",
              group_col, best_adj, paste(best_vars, collapse = ", ")))
  
  form  <- as.formula(paste("dS2S_SA_T10_000 ~", paste(best_vars, collapse = " + ")))
  model <- lm(form, data = dfp)
  
  cat(sprintf("\n--- Model Summary for %s ---\n", group_col))
  print(summary(model))
  
  list(model = model, data = dfp, vars = best_vars)
}

out1 <- run_litho_exh(itac, "Litho_Group_1")
out2 <- run_litho_exh(itac, "Litho_Group_2")
out3 <- run_litho_exh(itac, "Litho_Group_3")


#Diagnostics
par(mfrow = c(2,2))
plot(out1$model, main = "Diagnostics: Group1")
plot(out2$model, main = "Diagnostics: Group2")
plot(out3$model, main = "Diagnostics: Group3")
par(mfrow = c(1,1))  


#ANOVA - group dummy vs. without 
cont_all <- c("sqrt_slope","inv_slope","sqrt_vs30","inv_vs30")

# Group1
cont1 <- intersect(out1$vars, cont_all)          # ??rn: c("sqrt_slope","inv_vs30")
dum1  <- setdiff(out1$vars, cont1)               # ??rn: "Litho_Group_1_STIFF.ROCKS"

form_base1 <- as.formula(
  paste("dS2S_SA_T10_000 ~", paste(cont1, collapse = " + "))
)
form_full1 <- as.formula(
  paste("dS2S_SA_T10_000 ~", paste(c(cont1, dum1), collapse = " + "))
)

base1 <- lm(form_base1, data = out1$data)
full1 <- lm(form_full1, data = out1$data)

cat("=== Litho_Group_1 Nested ANOVA ===\n")
print(anova(base1, full1))


# Group2
cont2 <- intersect(out2$vars, cont_all)
dum2  <- setdiff(out2$vars, cont2)

form_base2 <- as.formula(
  paste("dS2S_SA_T10_000 ~", paste(cont2, collapse = " + "))
)
form_full2 <- as.formula(
  paste("dS2S_SA_T10_000 ~", paste(c(cont2, dum2), collapse = " + "))
)

base2 <- lm(form_base2, data = out2$data)
full2 <- lm(form_full2, data = out2$data)

cat("\n=== Litho_Group_2 Nested ANOVA ===\n")
print(anova(base2, full2))


# Group3
cont3 <- intersect(out3$vars, cont_all)
dum3  <- setdiff(out3$vars, cont3)

form_base3 <- as.formula(
  paste("dS2S_SA_T10_000 ~", paste(cont3, collapse = " + "))
)
form_full3 <- as.formula(
  paste("dS2S_SA_T10_000 ~", paste(c(cont3, dum3), collapse = " + "))
)

base3 <- lm(form_base3, data = out3$data)
full3 <- lm(form_full3, data = out3$data)

cat("\n=== Litho_Group_3 Nested ANOVA ===\n")
print(anova(base3, full3))
#Periodic AIC / BIC / Adj R?? Table

library(dplyr)
library(tidyr)
library(readr)

y_vars <- c(
  "dS2S_SA_T0_010", "dS2S_SA_T0_025", "dS2S_SA_T0_040", "dS2S_SA_T0_050", "dS2S_SA_T0_070",
  "dS2S_SA_T0_100", "dS2S_SA_T0_150", "dS2S_SA_T0_200", "dS2S_SA_T0_250", "dS2S_SA_T0_300",
  "dS2S_SA_T0_350", "dS2S_SA_T0_400", "dS2S_SA_T0_450", "dS2S_SA_T0_500", "dS2S_SA_T0_600",
  "dS2S_SA_T0_700", "dS2S_SA_T0_750", "dS2S_SA_T0_800", "dS2S_SA_T0_900", "dS2S_SA_T1_000",
  "dS2S_SA_T1_200", "dS2S_SA_T1_400", "dS2S_SA_T1_600", "dS2S_SA_T1_800", "dS2S_SA_T2_000",
  "dS2S_SA_T2_500", "dS2S_SA_T3_000", "dS2S_SA_T3_500", "dS2S_SA_T4_000", "dS2S_SA_T4_500",
  "dS2S_SA_T5_000", "dS2S_SA_T6_000", "dS2S_SA_T7_000", "dS2S_SA_T8_000", "dS2S_SA_T9_000",
  "dS2S_SA_T10_000", "dS2S_SA_pga"
)

# Getting related model metrics
get_periodic_model_metrics <- function(df, group_var) {
  res <- tibble::tibble()
  dummy_vars <- grep(paste0("^", group_var, "_"), names(df), value = TRUE)
  
  for (y in y_vars) {
    preds <- c("sqrt_slope", "inv_vs30", dummy_vars)
    preds_bt <- paste0("`", preds, "`")
    form <- as.formula(paste0("`", y, "` ~ ", paste(preds_bt, collapse = " + ")))
    
    dfm <- df %>% select(all_of(c(y, preds))) %>% na.omit()
    if (nrow(dfm) > 50) {
      m <- lm(form, data = dfm)
      res <- bind_rows(res, tibble::tibble(
        period = y,
        group  = group_var,
        AIC    = AIC(m),
        BIC    = BIC(m),
        adj_r2 = summary(m)$adj.r.squared
      ))
    }
  }
  return(res)
}


metrics1 <- get_periodic_model_metrics(out1$data, "Litho_Group_1")
metrics2 <- get_periodic_model_metrics(out2$data, "Litho_Group_2")
metrics3 <- get_periodic_model_metrics(out3$data, "Litho_Group_3")


all_metrics <- bind_rows(metrics1, metrics2, metrics3)

wide_table <- all_metrics %>%
  pivot_wider(names_from = group, values_from = c(adj_r2, AIC, BIC)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
print(wide_table)
write_csv(wide_table, "model_comparison_table_rounded.csv")



# #Exhaustive search -- for short period behaviour 
# run_litho_exh <- function(df, group_col) {
#   dfp <- df %>%
#     filter(!is.na(.data[[group_col]]), !is.na(dS2S_SA_pga)) %>%
#     dummy_cols(
#       select_columns         = group_col,
#       remove_first_dummy     = TRUE,    # k???1 dummy
#       remove_selected_columns = TRUE
#     ) %>%
#     rename_with(make.names)
#   
#   slope_vars <- c("sqrt_slope", "inv_slope")
#   vs30_vars  <- c("sqrt_vs30",  "inv_vs30")
#   dummy_vars <- grep(paste0("^", make.names(group_col), "_"), names(dfp), value = TRUE)
#   
#   best_adj  <- -Inf
#   best_vars <- NULL
#   
#   for(s in slope_vars) for(v in vs30_vars) {
#     vars <- c(s, v, dummy_vars)
#     sub  <- dfp %>% select(all_of(vars), dS2S_SA_pga) %>% na.omit()
#     if(nrow(sub) > 50) {
#       ex  <- regsubsets(x = sub[, vars], y = sub$dS2S_SA_pga,
#                         nbest = 1, nvmax = length(vars), method = "exhaustive")
#       sm  <- summary(ex)
#       bi  <- which.max(sm$adjr2)
#       sel <- names(which(sm$which[bi,]))[-1]
#       if(sm$adjr2[bi] > best_adj) {
#         best_adj  <- sm$adjr2[bi]
#         best_vars <- sel
#       }
#     }
#   }
#   
#   # Print selected predictors and model summary
#   cat(sprintf("\n???? [%s] Best predictors (Adj R?? = %.4f):\n  %s\n",
#               group_col, best_adj, paste(best_vars, collapse = ", ")))
#   
#   form  <- as.formula(paste("dS2S_SA_pga ~", paste(best_vars, collapse = " + ")))
#   model <- lm(form, data = dfp)
#   
#   cat(sprintf("\n--- Model Summary for %s ---\n", group_col))
#   print(summary(model))
#   
#   list(model = model, data = dfp, vars = best_vars)
# }
# 
# out1 <- run_litho_exh(itac, "Litho_Group_1")
# out2 <- run_litho_exh(itac, "Litho_Group_2")
# out3 <- run_litho_exh(itac, "Litho_Group_3")
