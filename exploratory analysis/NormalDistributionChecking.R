#Date of last modification : 30 march
#Your name ; Mattia De Bartolomeis
#which dataset you used : ITACA
#Topics you tackled :  Analysis of the normality of dS2S_SA over different time periods.
#Summary of main results : Some distributions conform to normality, while others do not.



# Carica il pacchetto necessario (se non già installato, decommenta la riga seguente)
# install.packages("readr")
#install.packages("readr")

library(readr)

# INSERT THE RIGHT PATH TO THE DATASET
file_path <- "/Users/mattiadebartolomeis/Desktop/APPLIED_STASTIC_PROJECT/ITACAs2s_2_flatfile/ITACAs2s_SA_2_0.csv"
data <- read.csv2(file_path, header = TRUE, stringsAsFactors = FALSE, dec = '.')

# Visualizza le prime righe del dataset
head(data)

# Visualizza la struttura del dataset
str(data)

# Statistiche descrittive per le variabili numeriche
summary(data)





#CHECK NORMALITY OF S2S WITH SHAPIRO TEST AND Q-Q PLOT
#Shapiro-Wilk Test:
# Hypothesis: The null hypothesis (H₀) of the Shapiro-Wilk test is that the data come from a normal distribution.
#p-value WITH 95% confidence level:
# If p-value > 0.05, you do not reject the null hypothesis, which suggests that the data do not show significant deviations from normality.
# If p-value < 0.05, you reject the null hypothesis, indicating that the data likely do not follow a normal distribution.

#QQ-Plot (Quantile-Quantile Plot):
#The QQ-plot compares the quantiles of your data with those of a normal distribution.
#If the points align along the straight line (in this case, the red line), the data distribution is compatible with a normal distribution.
#If the points deviate noticeably from the line (especially at the extremes), this indicates a departure from normality.

# List of columns for the periods
period_cols <- c("dS2S_SA_T0_010", "dS2S_SA_T0_025", "dS2S_SA_T0_040", "dS2S_SA_T0_050",
                 "dS2S_SA_T0_070", "dS2S_SA_T0_100", "dS2S_SA_T0_150", "dS2S_SA_T0_200",
                 "dS2S_SA_T0_250", "dS2S_SA_T0_300", "dS2S_SA_T0_350", "dS2S_SA_T0_400",
                 "dS2S_SA_T0_450", "dS2S_SA_T0_500", "dS2S_SA_T0_600", "dS2S_SA_T0_700",
                 "dS2S_SA_T0_750", "dS2S_SA_T0_800", "dS2S_SA_T0_900", "dS2S_SA_T1_000",
                 "dS2S_SA_T1_200", "dS2S_SA_T1_400", "dS2S_SA_T1_600", "dS2S_SA_T1_800",
                 "dS2S_SA_T2_000", "dS2S_SA_T2_500", "dS2S_SA_T3_000", "dS2S_SA_T3_500",
                 "dS2S_SA_T4_000", "dS2S_SA_T4_500", "dS2S_SA_T5_000", "dS2S_SA_T6_000",
                 "dS2S_SA_T7_000", "dS2S_SA_T8_000", "dS2S_SA_T9_000", "dS2S_SA_T10_000")

# Create an empty data frame to store p-values
results <- data.frame(period = character(), p_value = numeric(), stringsAsFactors = FALSE)

# Loop over each period column to perform the Shapiro-Wilk test and generate QQ-plots
for (col in period_cols) {
  # Extract the data vector
  vec <- na.omit(data[[col]])
  
  # Check that the sample is large enough (at least 3 values)
  if (length(vec) >= 3) {
    # Perform the Shapiro-Wilk test
    sh_test <- shapiro.test(vec)
    cat("Results of the Shapiro-Wilk test for", col, ":\n")
    print(sh_test)
    cat("\n")
    
    # Store the p-value in the results data frame
    results <- rbind(results, data.frame(period = col, p_value = sh_test$p.value))
    
    # Generate the QQ-plot
    qqnorm(vec, main = paste("QQ-Plot for", col))
    qqline(vec, col = "red")
    
    # Pause to view each plot if running interactively
    readline(prompt = "Press [Enter] to continue...")
  } else {
    cat("The variable", col, "has too few observations to run the test.\n\n")
  }
}




#PLOT EACH S_S TIME PERIOD TO IT'S P-VALUE IN THE SHAPIRO TEST

# Extract numeric time period from the column names by removing the prefix 'dS2S_SA_T'
# Replace the underscore with a decimal point, then convert to numeric
results$time_period <- as.numeric(gsub("_", ".", gsub("dS2S_SA_T", "", results$period)))

# Load ggplot2 for plotting
library(ggplot2)

# Plot the Shapiro-Wilk p-values as a function of the time period -> 
#values above this line suggest normality, while values below suggest non-normality.
ggplot(results, aes(x = time_period, y = p_value)) +
  geom_line() +
  geom_point() +
  labs(title = "Shapiro-Wilk p-values vs. Time Period",
       x = "Time Period",
       y = "p-value") +
  theme_minimal() +
  # Add a horizontal dashed red line at the 0.05 significance level
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")






