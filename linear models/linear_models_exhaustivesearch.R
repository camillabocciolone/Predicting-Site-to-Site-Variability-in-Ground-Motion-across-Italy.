#FORWARD SELECTION WITH AIC 
if (!require("olsrr")) install.packages("olsrr")
if (!require("dplyr")) install.packages("dplyr")
if (!require("readr")) install.packages("readr")
if (!require("fastDummies")) install.packages("fastDummies")

library(olsrr)
library(dplyr)
library(readr)
library(fastDummies)

set.seed(122)

# Caricamento dataset
df <- read_csv("/Users/mattiadebartolomeis/Desktop/APPLIED_STASTIC_PROJECT/ITACAs2s_2_flatfile/FINAL_DATASET.csv")

# Encoding variabili categoriche litologiche
df <- df %>% select(-c(Lito_ISPRA_bilinear, Lito_Bucci_bilinear))
df$Lito_ISPRA_simple <- factor(df$Lito_ISPRA_simple)
df$Lito_Bucci_simple <- factor(df$Lito_Bucci_simple)
df <- dummy_cols(df, select_columns = c("Lito_ISPRA_simple", "Lito_Bucci_simple"),
                 remove_first_dummy = TRUE, remove_selected_columns = TRUE)

# Rimozione colonne vietate e preparazione target
excluded_prefixes <- c("phi_ss", "fa_", "dS2S_SA_")
numerical_vars <- df %>% select(where(is.numeric))
valid_vars <- numerical_vars %>% select(-matches(paste(excluded_prefixes, collapse = "|")))
target_col <- "dS2S_SA_T10_000"
df_model <- bind_cols(valid_vars, target = df[[target_col]]) %>% na.omit()

# Standardizzazione
df_scaled <- df_model %>%
  mutate(across(-target, ~ as.numeric(scale(.))))

# Rimozione colonne tutte NA o con varianza nulla
df_scaled <- df_scaled %>%
  select(where(~!all(is.na(.)))) %>%
  select(where(~sd(., na.rm = TRUE) > 0))

# Suddivisione in train/test
n <- nrow(df_scaled)
train_indices <- sample(seq_len(n), size = floor(0.8 * n))
train_data <- df_scaled[train_indices, ]
test_data  <- df_scaled[-train_indices, ]

# Funzione per selezione feature per gruppi
filter_by_group <- function(group_pattern, df, threshold = 0.3) {
  group_features <- grep(group_pattern, colnames(df), ignore.case = TRUE, value = TRUE)
  if (length(group_features) == 0) return(character(0))
  corrs <- sapply(group_features, function(f) cor(df[[f]], df$target, use = "complete.obs"))
  corrs <- corrs[!is.na(corrs)]
  best <- names(which.max(abs(corrs)))
  if (!is.na(best) && abs(corrs[best]) >= threshold) return(best)
  return(character(0))
}

# Selezione feature
best_vs30  <- filter_by_group("vs30|vseq_ntc|vs_eq_ntc", train_data)
best_lito  <- filter_by_group("lito", train_data)
best_slope <- filter_by_group("slope", train_data)
best_depth <- filter_by_group("depth_to_bedrock", train_data)
nrec_cols  <- c("nrec_min", "nrec_max")
high_corr_features <- names(which(abs(cor(train_data))["target", ] >= 0.3))
nrec_used <- intersect(nrec_cols, high_corr_features)

all_group_cols <- c(best_vs30, best_lito, best_slope, best_depth, nrec_used)
final_features <- unique(all_group_cols)
if (length(final_features) == 0) stop("Nessuna feature selezionata dopo i filtri.")

# Separazione numeriche e dummy
numeric_features <- setdiff(final_features, grep("Lito", final_features, value = TRUE))
dummy_features <- grep("Lito", final_features, value = TRUE)

# Offset per trasformazioni
epsilon <- 1e-3
min_vals <- abs(min(train_data[, numeric_features], na.rm = TRUE)) + epsilon

# Feature derivate + interazioni tra numeriche
numeric_features_extended <- c(
  numeric_features,
  paste0("I(", numeric_features, "^2)"),
  paste0("I(log(", numeric_features, " + ", min_vals, "))"),
  paste0("I(1 / (", numeric_features, " + ", min_vals, "))"),
  paste0("I(sqrt(", numeric_features, " + ", min_vals, "))"),
  paste0("I(exp(", numeric_features, "))")
)

# Aggiunta interazioni tra variabili numeriche (semplici)
combinations_numeric <- combn(numeric_features, 2, simplify = TRUE)
interaction_terms <- apply(combinations_numeric, 2, function(x) paste0("I(", x[1], " * ", x[2], ")"))
numeric_features_extended <- c(numeric_features_extended, interaction_terms)

# Formula del modello
features_extended <- c(numeric_features_extended, dummy_features)
form <- as.formula(paste("target ~", paste(features_extended, collapse = " + ")))

# Modello completo + selezione
model_full <- lm(form, data = train_data)
model_selected <- ols_step_forward_aic(model_full, progress = TRUE)

# Predizione sul test
preds_test <- predict(model_selected$model, newdata = test_data)

# Valutazione performance
# Valutazione performance (ignorando righe con predizioni NA)
actual_test <- test_data$target

# Indici delle predizioni valide (non NA)
valid_idx <- which(!is.na(preds_test))

# Calcolo RMSE e R² solo sui valori validi
rmse_test <- sqrt(mean((preds_test[valid_idx] - actual_test[valid_idx])^2))
ss_total <- sum((actual_test[valid_idx] - mean(actual_test[valid_idx]))^2)
ss_res <- sum((actual_test[valid_idx] - preds_test[valid_idx])^2)
r2_test <- 1 - (ss_res / ss_total)


# Output finale
cat("Modello finale selezionato:\n")
print(summary(model_selected$model))

cat("Valutazione modello sul test set:\n")
cat("RMSE test:", rmse_test, "\n")
cat("R^2 test:", r2_test, "\n")



#WITH FORWARD SELECTION WITH CROSS VALIDATION 
# # Pacchetti richiesti
# if (!require("dplyr")) install.packages("dplyr")
# if (!require("readr")) install.packages("readr")
# if (!require("fastDummies")) install.packages("fastDummies")
# 
# library(dplyr)
# library(readr)
# library(fastDummies)
# 
# set.seed(122)
# 
# # Caricamento dataset
# df <- read_csv("/Users/mattiadebartolomeis/Desktop/APPLIED_STASTIC_PROJECT/ITACAs2s_2_flatfile/FINAL_DATASET.csv")
# 
# # Encoding variabili categoriche litologiche
# df <- df %>% select(-c(Lito_ISPRA_bilinear, Lito_Bucci_bilinear))
# df$Lito_ISPRA_simple <- factor(df$Lito_ISPRA_simple)
# df$Lito_Bucci_simple <- factor(df$Lito_Bucci_simple)
# df <- dummy_cols(df, select_columns = c("Lito_ISPRA_simple", "Lito_Bucci_simple"),
#                  remove_first_dummy = TRUE, remove_selected_columns = TRUE)
# 
# # Rimozione colonne vietate e preparazione target
# excluded_prefixes <- c("phi_ss", "fa_", "dS2S_SA_")
# numerical_vars <- df %>% select(where(is.numeric))
# valid_vars <- numerical_vars %>% select(-matches(paste(excluded_prefixes, collapse = "|")))
# target_col <- "dS2S_SA_T10_000"
# df_model <- bind_cols(valid_vars, target = df[[target_col]]) %>% na.omit()
# 
# # Standardizzazione
# df_scaled <- df_model %>% mutate(across(-target, ~ as.numeric(scale(.))))
# 
# # Rimozione colonne tutte NA o con varianza nulla
# df_scaled <- df_scaled %>% select(where(~!all(is.na(.)))) %>% select(where(~sd(., na.rm = TRUE) > 0))
# 
# # Suddivisione in train/test
# n <- nrow(df_scaled)
# train_indices <- sample(seq_len(n), size = floor(0.8 * n))
# train_data <- df_scaled[train_indices, ]
# test_data  <- df_scaled[-train_indices, ]
# 
# # Funzione per selezione feature per gruppi
# filter_by_group <- function(group_pattern, df, threshold = 0.3) {
#   group_features <- grep(group_pattern, colnames(df), ignore.case = TRUE, value = TRUE)
#   if (length(group_features) == 0) return(character(0))
#   corrs <- sapply(group_features, function(f) cor(df[[f]], df$target, use = "complete.obs"))
#   corrs <- corrs[!is.na(corrs)]
#   best <- names(which.max(abs(corrs)))
#   if (!is.na(best) && abs(corrs[best]) >= threshold) return(best)
#   return(character(0))
# }
# 
# # Selezione feature per gruppo
# best_vs30  <- filter_by_group("vs30|vseq_ntc|vs_eq_ntc", train_data)
# best_lito  <- filter_by_group("lito", train_data)
# best_slope <- filter_by_group("slope", train_data)
# best_depth <- filter_by_group("depth_to_bedrock", train_data)
# nrec_cols  <- c("nrec_min", "nrec_max")
# high_corr_features <- names(which(abs(cor(train_data))['target', ] >= 0.3))
# nrec_used <- intersect(nrec_cols, high_corr_features)
# 
# all_group_cols <- c(best_vs30, best_lito, best_slope, best_depth, nrec_used)
# final_features <- unique(all_group_cols)
# if (length(final_features) == 0) stop("Nessuna feature selezionata dopo i filtri.")
# 
# # Feature numeriche e dummy
# numeric_features <- setdiff(final_features, grep("Lito", final_features, value = TRUE))
# dummy_features <- grep("Lito", final_features, value = TRUE)
# 
# # Trasformazioni e interazioni
# epsilon <- 1e-3
# min_vals <- abs(min(train_data[, numeric_features], na.rm = TRUE)) + epsilon
# numeric_features_extended <- c(
#   numeric_features,
#   paste0("I(", numeric_features, "^2)"),
#   paste0("I(log(", numeric_features, " + ", min_vals, "))"),
#   paste0("I(1 / (", numeric_features, " + ", min_vals, "))"),
#   paste0("I(sqrt(", numeric_features, " + ", min_vals, "))"),
#   paste0("I(exp(", numeric_features, "))")
# )
# combinations_numeric <- combn(numeric_features, 2, simplify = TRUE)
# interaction_terms <- apply(combinations_numeric, 2, function(x) paste0("I(", x[1], " * ", x[2], ")"))
# numeric_features_extended <- c(numeric_features_extended, interaction_terms)
# features_extended <- c(numeric_features_extended, dummy_features)
# 
# # Funzione di validazione incrociata
# cv_glm_error <- function(data, formula, k = 5) {
#   n <- nrow(data)
#   folds <- sample(rep(1:k, length.out = n))
#   errors <- numeric(k)
# 
#   for (i in 1:k) {
#     train_fold <- data[folds != i, ]
#     test_fold  <- data[folds == i, ]
#     model <- lm(formula, data = train_fold)
#     preds <- predict(model, newdata = test_fold)
#     valid_idx <- which(!is.na(preds))
#     errors[i] <- mean((preds[valid_idx] - test_fold$target[valid_idx])^2)
#   }
#   return(mean(errors))
# }
# 
# # Forward selection con cross-validation
# selected_names <- c()
# remaining <- features_extended
# best_cv_error <- Inf
# 
# repeat {
#   best_feature <- NULL
#   best_feature_error <- Inf
# 
#   for (feat in remaining) {
#     current_features <- c(selected_names, feat)
#     current_formula <- as.formula(paste("target ~", paste(current_features, collapse = " + ")))
#     cv_error <- cv_glm_error(train_data, current_formula, k = 5)
# 
#     if (cv_error < best_feature_error) {
#       best_feature <- feat
#       best_feature_error <- cv_error
#     }
#   }
# 
#   if (length(selected_names) == 0 && !is.null(best_feature)) {
#     selected_names <- c(selected_names, best_feature)
#     remaining <- setdiff(remaining, best_feature)
#     best_cv_error <- best_feature_error
#   } else if (!is.null(best_feature) && best_feature_error < best_cv_error) {
#     selected_names <- c(selected_names, best_feature)
#     remaining <- setdiff(remaining, best_feature)
#     best_cv_error <- best_feature_error
#   } else {
#     break
#   }
# }
# 
# # Costruzione modello finale e valutazione su test
# if (length(selected_names) == 0) stop("Nessuna feature selezionata dalla forward selection.")
# form_final <- as.formula(paste("target ~", paste(selected_names, collapse = " + ")))
# model_final <- lm(form_final, data = train_data)
# preds_test <- predict(model_final, newdata = test_data)
# actual_test <- test_data$target
# valid_idx <- which(!is.na(preds_test))
# rmse_test <- sqrt(mean((preds_test[valid_idx] - actual_test[valid_idx])^2))
# ss_total <- sum((actual_test[valid_idx] - mean(actual_test[valid_idx]))^2)
# ss_res <- sum((actual_test[valid_idx] - preds_test[valid_idx])^2)
# r2_test <- 1 - (ss_res / ss_total)
# 
# cat("Modello finale con cross-validation:\n")
# print(summary(model_final))
# cat("\nValutazione sul test set:\n")
# cat("RMSE:", rmse_test, "\n")
# cat("R^2:", r2_test, "\n")