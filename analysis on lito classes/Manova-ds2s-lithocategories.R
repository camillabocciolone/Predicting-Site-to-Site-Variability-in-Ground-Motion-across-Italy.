
# manova to see if the classes were statistically significant over the prediction 
# of ds2s 

# we also used bonferroni confidence intervals to compare pair of lito classes
library(dplyr)

# INSERT THE RIGHT PATH TO THE DATASET
df_itaca = read.csv("ITACA_and_maps_V1.csv", sep=" ")

column_names = names(df_itaca) # Get all the column numbers
column_names

variabili_dipendenti <- c(
  "dS2S_SA_T0_010", "dS2S_SA_T0_025", "dS2S_SA_T0_040", "dS2S_SA_T0_050", "dS2S_SA_T0_070",
  "dS2S_SA_T0_100", "dS2S_SA_T0_150", "dS2S_SA_T0_200", "dS2S_SA_T0_250", "dS2S_SA_T0_300",
  "dS2S_SA_T0_350", "dS2S_SA_T0_400", "dS2S_SA_T0_450", "dS2S_SA_T0_500", "dS2S_SA_T0_600",
  "dS2S_SA_T0_700", "dS2S_SA_T0_750", "dS2S_SA_T0_800", "dS2S_SA_T0_900", "dS2S_SA_T1_000",
  "dS2S_SA_T1_200", "dS2S_SA_T1_400", "dS2S_SA_T1_600", "dS2S_SA_T1_800", "dS2S_SA_T2_000",
  "dS2S_SA_T2_500", "dS2S_SA_T3_000", "dS2S_SA_T3_500", "dS2S_SA_T4_000", "dS2S_SA_T4_500",
  "dS2S_SA_T5_000", "dS2S_SA_T6_000", "dS2S_SA_T7_000", "dS2S_SA_T8_000", "dS2S_SA_T9_000",
  "dS2S_SA_T10_000", "dS2S_SA_pgv"
)

#MANOVA on ispra categories

fattore_lito <- "Lito_ISPRA_simple" 

ggplot(df_itaca, aes(x = as.factor(Lito_ISPRA_simple), y = dS2S_SA_T1_000, fill = as.factor(Lito_ISPRA_simple))) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.6) +
  theme_minimal() +
  labs(title = "Boxplot della proxy dS2S_SA_T1_000 per le 50 classi litologiche ISPRA",
       x = "Classe litologica ISPRA",
       y = "dS2S_SA_T1_000") +
  theme(axis.text.x = element_text(angle = 90, size = 8),
        legend.position = "none")


# Crea un subset con le variabili necessarie
df <- df_itaca[, c(fattore_lito, variabili_dipendenti)]

# Elimina eventuali NA
df <- na.omit(df)

# Assicurati che il fattore sia "factor"
df[[fattore_lito]] <- as.factor(df[[fattore_lito]])

# Statistiche base
g <- nlevels(df[[fattore_lito]])
p <- length(variabili_dipendenti)
n <- nrow(df)

### VERIFICA DELLE ASSUNZIONI

# 1. Normalità multivariata (se possibile per dimensioni gruppi)
library(MVN)
Ps <- sapply(levels(df[[fattore_lito]]), function(lv) {
  gruppo <- df[df[[fattore_lito]] == lv, variabili_dipendenti]
  if (nrow(gruppo) >= p + 1) {
    mvn.test <- mvn(data = gruppo, mvnTest = "hz")
    return(mvn.test$multivariateNormality$`p value`)
  } else {
    return(NA)  # Troppi pochi dati per test
  }
})
Ps

# Test di omogeneità delle matrici di covarianza (Box's M test)
#summary(boxM(df[, variabili_dipendenti], grouping = df[[fattore_lito]]))


#MANOVA multivariate 
# Formula
form <- as.formula(paste("cbind(", paste(variabili_dipendenti, collapse = ", "), ") ~ ", fattore_lito))

manova_model <- manova(formula = form, data = df)

summary(manova_model, test = "Wilks")

library(dplyr)

# Calcolo dimensioni
n_tot <- nrow(df)
gruppi <- levels(df[[fattore_lito]])
g <- length(gruppi)
p <- length(variabili_dipendenti)

# Estrai medie e numerosità per ciascun gruppo
medie_per_gruppo <- df %>%
  group_by(!!as.symbol(fattore_lito)) %>%
  summarise(across(all_of(variabili_dipendenti), mean))

n_i <- df %>%
  group_by(!!as.symbol(fattore_lito)) %>%
  summarise(n = n()) %>%
  pull(n)

names(n_i) <- gruppi

# Matrice dei residui
W <- summary(manova_model)$SS$Residuals

# Bonferroni correction
alpha <- 0.05
k <- p * choose(g, 2)  # numero totale di confronti
qT <- qt(1 - alpha / (2 * k), df = n_tot - g)

# Funzione per confronto a coppie
confronti_bonferroni <- function(gr1, gr2) {
  m1 <- unlist(medie_per_gruppo[medie_per_gruppo[[fattore_lito]] == gr1, -1])
  m2 <- unlist(medie_per_gruppo[medie_per_gruppo[[fattore_lito]] == gr2, -1])
  n1 <- n_i[gr1]
  n2 <- n_i[gr2]
  
  diff <- m1 - m2
  se <- sqrt(diag(W) / (n_tot - g) * (1/n1 + 1/n2))
  
  lower <- diff - qT * se
  upper <- diff + qT * se
  
  sig <- !(lower < 0 & upper > 0)
  
  tibble(
    gruppo1 = gr1,
    gruppo2 = gr2,
    variabile = variabili_dipendenti,
    diff = diff,
    lower = lower,
    upper = upper,
    significativo = sig
  )
}

# Esegui confronti su tutte le coppie
risultati <- purrr::map_dfr(combn(gruppi, 2, simplify = FALSE),
                            ~confronti_bonferroni(.x[1], .x[2]))

# Filtra solo confronti significativi
risultati_significativi <- risultati %>% filter(significativo == TRUE)

# Visualizza un estratto
head(risultati_significativi)

# Filtra i confronti non significativi
non_significativi <- risultati[risultati$significativo == FALSE, ]

# Mostra i confronti non significativi
head(non_significativi)

# Raggruppa per coppie di gruppi e verifica se c'è almeno una variabile significativa per ciascuna
library(dplyr)

coppie_significative <- risultati %>%
  group_by(gruppo1, gruppo2) %>%
  summarise(almeno_una_significativa = any(significativo), .groups = "drop") %>%
  filter(almeno_una_significativa == TRUE) %>%
  select(gruppo1, gruppo2)

coppie_significative

# ANOVA univariate per ciascuna variabile
summary.aov(manova_model)

#MANOVA on BUcci classes

ggplot(df_itaca, aes(x = as.factor(Lito_Bucci_simple), y = dS2S_SA_T1_000, fill = as.factor(Lito_Bucci_simple))) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.6) +
  theme_minimal() +
  labs(title = "Boxplot della proxy dS2S_SA_T1_000 per le 11 classi litologiche Bucci",
       x = "Classe litologica Bucci",
       y = "dS2S_SA_T1_000") +
  theme(axis.text.x = element_text(angle = 90, size = 8),
        legend.position = "none")

fattore_lito <- "Lito_Bucci_simple" 

# Crea un subset con le variabili necessarie
df <- df_itaca[, c(fattore_lito, variabili_dipendenti)]

# Elimina eventuali NA
df <- na.omit(df)

# Assicurati che il fattore sia "factor"
df[[fattore_lito]] <- as.factor(df[[fattore_lito]])

#MANOVA  
# Formula
form <- as.formula(paste("cbind(", paste(variabili_dipendenti, collapse = ", "), ") ~ ", fattore_lito))

# MANOVA
manova_model <- manova(formula = form, data = df)

# Risultato con test di Wilks
summary(manova_model, test = "Wilks")

# ANOVA univariate per ciascuna variabile
summary.aov(manova_model)


#vs30 from profile vs ISPRA
df_itaca_clean <- df_itaca[!is.na(df_itaca$vs30_from_profile)& df_itaca$vs30_from_profile<1500, ]


ggplot(df_itaca_clean, aes(x = as.factor(Lito_ISPRA_simple), y = vs30_from_profile, fill = as.factor(Lito_ISPRA_simple))) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.6) +
  theme_minimal() +
  labs(title = "Boxplot della proxy vs_30_from_profile per le 50 classi litologiche ISPRA",
       x = "Classe litologica ISPRA",
       y = "vs_30_from_profile") +
  theme(axis.text.x = element_text(angle = 90, size = 8),
        legend.position = "none")



# ANOVA vera e propria
df_itaca_clean$Lito_ISPRA_simple <- as.factor(df_itaca_clean$Lito_ISPRA_simple)
anova_model_ispra <- aov(vs30_from_profile ~ Lito_ISPRA_simple, data = df_itaca_clean)
summary(anova_model_ispra)

----#vs30 from profile vs Bucci

ggplot(df_itaca_clean, aes(x = as.factor(Lito_Bucci_simple), y = vs30_from_profile, fill = as.factor(Lito_Bucci_simple))) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.6) +
  theme_minimal() +
  labs(title = "Boxplot della proxy vs_30_from_profile per le 11 classi litologiche Bucci",
       x = "Classe litologica Bucci",
       y = "vs_30_from_profile") +
  theme(axis.text.x = element_text(angle = 90, size = 8),
        legend.position = "none")

# ANOVA vera e propria
df_itaca_clean$Lito_Bucci_simple <- as.factor(df_itaca_clean$Lito_Bucci_simple)
anova_model <- aov(vs30_from_profile ~ Lito_Bucci_simple, data = df_itaca_clean)
summary(anova_model)


#vs30_m_s_from_topography vs ISPRA
df_itaca_clean <- df_itaca[!is.na(df_itaca$vs30_m_s_from_topography), ]


ggplot(df_itaca_clean, aes(x = as.factor(Lito_ISPRA_simple), y = vs30_m_s_from_topography, fill = as.factor(Lito_ISPRA_simple))) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.6) +
  theme_minimal() +
  labs(title = "Boxplot della proxy vs30_m_s_from_topography per le 50 classi litologiche ISPRA",
       x = "Classe litologica ISPRA",
       y = "vs30_m_s_from_topography") +
  theme(axis.text.x = element_text(angle = 90, size = 8),
        legend.position = "none")

# ANOVA vera e propria
df_itaca_clean$Lito_ISPRA_simple <- as.factor(df_itaca_clean$Lito_ISPRA_simple)
anova_model_ispra <- aov(vs30_m_s_from_topography ~ Lito_ISPRA_simple, data = df_itaca_clean)
summary(anova_model_ispra)

#vs30_m_s_from_topography vs Bucci
ggplot(df_itaca_clean, aes(x = as.factor(Lito_Bucci_simple), y = vs30_m_s_from_topography, fill = as.factor(Lito_Bucci_simple))) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.6) +
  theme_minimal() +
  labs(title = "Boxplot della proxy vs30_m_s_from_topography per le 11 classi litologiche Bucci",
       x = "Classe litologica Bucci",
       y = "vs30_m_s_from_topography") +
  theme(axis.text.x = element_text(angle = 90, size = 8),
        legend.position = "none")

# ANOVA vera e propria
df_itaca_clean$Lito_Bucci_simple <- as.factor(df_itaca_clean$Lito_Bucci_simple)
anova_model <- aov(vs30_m_s_from_topography ~ Lito_Bucci_simple, data = df_itaca_clean)
summary(anova_model)


#vs30_m_s_from_Mori vs ISPRA
df_itaca_clean <- df_itaca[!is.na(df_itaca$Vs30_MORI_simple) & df_itaca$Vs30_MORI_simple > 5, ]


ggplot(df_itaca_clean, aes(x = as.factor(Lito_ISPRA_simple), y = Vs30_MORI_simple, fill = as.factor(Lito_ISPRA_simple))) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.6) +
  theme_minimal() +
  labs(title = "Boxplot della proxy vs30_m_s_Mori per le 50 classi litologiche ISPRA",
       x = "Classe litologica ISPRA",
       y = "Vs30_MORI_simple") +
  theme(axis.text.x = element_text(angle = 90, size = 8),
        legend.position = "none")

# ANOVA vera e propria
df_itaca_clean$Lito_ISPRA_simple <- as.factor(df_itaca_clean$Lito_ISPRA_simple)
anova_model_ispra <- aov(Vs30_MORI_simple ~ Lito_ISPRA_simple, data = df_itaca_clean)
summary(anova_model_ispra)



#vs30_m_s_from_MOri vs Bucci
df_itaca_clean <- df_itaca[!is.na(df_itaca$Vs30_MORI_simple) &
                             df_itaca$Vs30_MORI_simple > 5 &
                             !is.na(df_itaca$Lito_Bucci_simple)&
                             df_itaca$Lito_Bucci_simple != 5&
                             df_itaca$Lito_Bucci_simple != 11, ]




ggplot(df_itaca_clean, aes(x = as.factor(Lito_Bucci_simple), y = Vs30_MORI_simple, fill = as.factor(Lito_Bucci_simple))) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.6) +
  theme_minimal() +
  labs(title = "Boxplot della proxy vs30_mori per le 11 classi litologiche Bucci",
       x = "Classe litologica Bucci",
       y = "Vs30_MORI_simple") +
  theme(axis.text.x = element_text(angle = 90, size = 8),
        legend.position = "none")

# ANOVA vera e propria
df_itaca_clean$Lito_Bucci_simple <- as.factor(df_itaca_clean$Lito_Bucci_simple)
anova_model<- aov(Vs30_MORI_simple ~ Lito_Bucci_simple, data = df_itaca_clean)
summary(anova_model)

