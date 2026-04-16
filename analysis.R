#1
AirQuality

#2
#install.packages("psych")
library(psych)

stats <- describe(AirQuality, skew = TRUE, ranges = TRUE, IQR = TRUE)
get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
stats_complete <- cbind(stats,
                        variance = stats$sd^2,
                        cv = round((stats$sd / stats$mean) * 100, 2),  # coefficient of variation (%)
                        mode = sapply(AirQuality, get_mode))
print(stats_complete[, c("mean", "median", "mode", "sd", "variance", 
                         "cv", "range", "skew", "kurtosis")])

#3 matrix
cor_matrix <- round(cor(AirQuality), 3)
print(cor_matrix)

#visiable
#install.packages("corrplot")
library(corrplot)
cor_with_quality_matrix <- as.matrix(cor_with_quality)
colnames(cor_with_quality_matrix) <- "CO(GT)"
corrplot(cor_with_quality_matrix[-12, , drop = FALSE],  
         method = "number",
         tl.cex = 0.8,
         mar = c(0, 0, 1, 0))

#4
model_full <- lm(`CO(GT)` ~ ., data = AirQuality)
print(model_full)

#5
#Intercept: The predicted CO(GT) is -49.97 when all predictors are zero.
#PT08.S1(CO): A one-unit increase in PT08.S1(CO) is associated with a 0.0038 unit decrease in CO(GT), on average.
#NMHC(GT): A one-unit increase in NMHC(GT) is associated with a 0.0402 unit increase in CO(GT), on average.
#C6H6(GT): A one-unit increase in C6H6(GT) is associated with a 0.3966 unit decrease in CO(GT), on average.
#PT08.S2(NMHC): A one-unit increase in PT08.S2(NMHC) is associated with a 0.0111 unit decrease in CO(GT), on average.
#NOx(GT): A one-unit increase in NOx(GT) is associated with a 0.0293 unit increase in CO(GT), on average.
#PT08.S3(NOx): A one-unit increase in PT08.S3(NOx) is associated with a 0.0116 unit increase in CO(GT), on average.
#NO2(GT): A one-unit increase in NO2(GT) is associated with a 0.3970 unit increase in CO(GT), on average.
#PT08.S4(NO2): A one-unit increase in PT08.S4(NO2) is associated with a 0.0020 unit decrease in CO(GT), on average.
#PT08.S5(O3): A one-unit increase in PT08.S5(O3) is associated with a 0.0108 unit decrease in CO(GT), on average.
#T: A one-unit increase in temperature is associated with a 0.4617 unit increase in CO(GT), on average.
#RH: A one-unit increase in relative humidity is associated with a 0.2141 unit increase in CO(GT), on average.
#AH: A one-unit increase in absolute humidity is associated with a 0.9803 unit decrease in CO(GT), on average.

#6
summary(model_full)
model_summary <- summary(model_full)
coefficients_se <- model_summary$coefficients[, 2]
print(round(coefficients_se, 6))

residual_se <- model_summary$sigma
print(paste("residual_se:", round(residual_se, 4)))

regression_table <- data.frame(
  variable = rownames(model_summary$coefficients),
  coefficient = round(model_summary$coefficients[, 1], 4),
  SE = round(model_summary$coefficients[, 2], 6),
  t_value = round(model_summary$coefficients[, 3], 4),
  p_value = round(model_summary$coefficients[, 4], 6),
  significant = ifelse(model_summary$coefficients[, 4] < 0.001, "***",
               ifelse(model_summary$coefficients[, 4] < 0.01, "**",
                      ifelse(model_summary$coefficients[, 4] < 0.05, "*", "not significant")))
)

print(regression_table)

significant_vars <- regression_table[regression_table$p_value < 0.05, ]
insignificant_vars <- regression_table[regression_table$p_value >= 0.05, ]

cat("\n【significant (p < 0.05)】\n")
print(significant_vars[, c("variable", "coefficient", "p_value", "significant")])

cat("\n【not significant (p >= 0.05)】\n")
if(nrow(insignificant_vars) > 0) {
  print(insignificant_vars[, c("variable", "coefficient", "p_value", "significant")])
} else {
  cat("all significant\n")
}

conf_intervals <- confint(model_full, level = 0.95)
print("95% significant level：")
print(round(conf_intervals, 4))

#7
r_squared <- model_summary$r.squared
adj_r_squared <- model_summary$adj.r.squared
cat(sprintf("R²: %.4f\n", r_squared))
cat(sprintf("adjusted R²: %.4f\n", adj_r_squared))
#Different formula and different interpretation
#R²: 0.4737   Show how much variance is explained
#adjusted R²: 0.4730    Check whether each variable worth including

#8
if (!require(car)) {
  install.packages("car")
  library(car)
}
vif_values <- vif(model_full)
print(vif_values)
tolerance <- 1 / vif_values

vif_interpretation <- data.frame(
  variable = names(vif_values),
  VIF = round(vif_values, 2),
  tolerance = round(1/vif_values, 4),
  severity = ifelse(vif_values > 10, "Severe multicollinearity ",
                    ifelse(vif_values > 5, "Moderate multicollinearity ",
                           ifelse(vif_values > 2.5, "Weak multicollinearity", "No multicollinearity ✓")))
)
print(vif_interpretation)
cat("\n========== VIF Interpretation Criteria ==========\n")
cat("VIF = 1          : No correlation\n")
cat("1 < VIF < 5      : Moderate correlation, acceptable\n")
cat("5 ≤ VIF < 10     : Strong correlation, needs attention\n")
cat("VIF ≥ 10         : Severe multicollinearity, must be addressed\n")

#9
predictors <- setdiff(names(AirQuality), "CO(GT)")
current_model <- NULL
current_adj_r2 <- -Inf
selected_vars <- c()
forward_steps <- list()

repeat {
  remaining_vars <- setdiff(predictors, selected_vars)
  if (length(remaining_vars) == 0) break
  
  best_var <- NULL
  best_adj_r2 <- current_adj_r2
  
  for (var in remaining_vars) {
    if (is.null(current_model)) {
      new_formula <- as.formula(paste("`CO(GT)` ~ `", var, "`", sep = ""))
    } else {
      selected_vars_backtick <- paste("`", selected_vars, "`", sep = "", collapse = " + ")
      new_formula <- as.formula(paste("`CO(GT)` ~ ", selected_vars_backtick, " + `", var, "`", sep = ""))
    }
    new_model <- lm(new_formula, data = AirQuality)
    new_adj_r2 <- summary(new_model)$adj.r.squared
    
    if (new_adj_r2 > best_adj_r2) {
      best_adj_r2 <- new_adj_r2
      best_var <- var
    }
  }
  
  if (best_adj_r2 > current_adj_r2) {
    current_adj_r2 <- best_adj_r2
    selected_vars <- c(selected_vars, best_var)
    selected_vars_backtick <- paste("`", selected_vars, "`", sep = "", collapse = " + ")
    current_model <- lm(as.formula(paste("`CO(GT)` ~ ", selected_vars_backtick, sep = "")), data = AirQuality)
    
    forward_steps[[length(forward_steps) + 1]] <- list(
      step = length(selected_vars),
      added = best_var,
      adj_r2 = current_adj_r2,
      vars = selected_vars
    )
    
    cat(sprintf("Step %d: Added '%s' → Adjusted R² = %.6f\n", 
                length(selected_vars), best_var, current_adj_r2))
  } else {
    break
  }
}

cat("\n========== Final Model ==========\n")
cat("Selected variables in order:\n")
print(selected_vars)
cat(sprintf("\nFinal Adjusted R²: %.6f\n", current_adj_r2))

final_forward_model <- lm(as.formula(paste("`CO(GT)` ~ ", paste("`", selected_vars, "`", sep = "", collapse = " + "), sep = "")), data = AirQuality)
cat("\nFinal model summary:\n")
print(summary(final_forward_model))

#9.1 Post-processing: remove non-significant variables from final forward model
cat("\n========== Post-processing: Remove Non-Significant Variables ==========\n")

final_summary <- summary(final_forward_model)
final_coef <- final_summary$coefficients

# Identify non-significant variables (p >= 0.05), excluding Intercept
non_significant_vars <- rownames(final_coef)[final_coef[, 4] >= 0.05 & rownames(final_coef) != "(Intercept)"]

if (length(non_significant_vars) > 0) {
  cat("Non-significant variables to remove (p ≥ 0.05):\n")
  print(non_significant_vars)
  
  # Remove them from selected_vars
  cleaned_vars <- setdiff(selected_vars, non_significant_vars)
  
  cat("\nVariables after removal:\n")
  print(cleaned_vars)
  
  # Refit model with only significant variables
  if (length(cleaned_vars) > 0) {
    cleaned_formula <- as.formula(paste("`CO(GT)` ~", paste("`", cleaned_vars, "`", sep = "", collapse = " + ")))
    final_cleaned_model <- lm(cleaned_formula, data = AirQuality)
    
    cat("\n========== Final Cleaned Model Summary ==========\n")
    print(summary(final_cleaned_model))
    
    cat("\n========== Adjusted R² Comparison ==========\n")
    cat(sprintf("Original forward model Adjusted R²: %.6f\n", summary(final_forward_model)$adj.r.squared))
    cat(sprintf("Cleaned model Adjusted R²: %.6f\n", summary(final_cleaned_model)$adj.r.squared))
    
    if (summary(final_cleaned_model)$adj.r.squared >= summary(final_forward_model)$adj.r.squared) {
      cat("Result: Cleaned model is preferred (equal or higher adjusted R²).\n")
    } else {
      cat("Result: Original forward model has slightly higher adjusted R², but cleaned model is simpler.\n")
    }
  } else {
    cat("\nWarning: All variables removed. Only intercept model remains.\n")
  }
} else {
  cat("All variables in forward model are significant (p < 0.05). No removal needed.\n")
}

#10
cat("========== Verification of Statistical Significance (α = 0.05) ==========\n\n")

# Significance of individual coefficients
cat("--- Individual Coefficient Significance ---\n")
cleaned_summary <- summary(final_cleaned_model)
coefficients <- cleaned_summary$coefficients

significance_check <- data.frame(
  variable = rownames(coefficients),
  coefficient = round(coefficients[, 1], 6),
  p_value = round(coefficients[, 4], 6),
  significant_at_5percent = ifelse(coefficients[, 4] < 0.05, "YES", "NO")
)
print(significance_check)

cat("\n--- Overall Equation Significance ---\n")
f_statistic <- cleaned_summary$fstatistic
f_p_value <- pf(f_statistic[1], f_statistic[2], f_statistic[3], lower.tail = FALSE)

cat(sprintf("F-statistic: %.4f\n", f_statistic[1]))
cat(sprintf("Degrees of freedom: df1 = %d, df2 = %d\n", f_statistic[2], f_statistic[3]))
cat(sprintf("p-value: %.6e\n", f_p_value))

if (f_p_value < 0.05) {
  cat("Conclusion: p-value < 0.05 → The equation as a whole is statistically significant at the 5% level.\n")
} else {
  cat("Conclusion: p-value ≥ 0.05 → The equation as a whole is NOT statistically significant at the 5% level.\n")
}

cat("\n--- Interpretation of Model Coefficient Estimates ---\n")
cat("(All interpretations assume other variables remain constant)\n\n")

for (i in 1:nrow(coefficients)) {
  var_name <- rownames(coefficients)[i]
  coef_val <- coefficients[i, 1]
  p_val <- coefficients[i, 4]
  
  if (var_name == "(Intercept)") {
    cat(sprintf("Intercept = %.6f: When all predictors are zero, the predicted CO(GT) is %.6f.\n", coef_val, coef_val))
  } else {
    if (coef_val > 0) {
      cat(sprintf("%s = %.6f: A one-unit increase in %s is associated with a %.6f unit INCREASE in CO(GT), on average.\n", 
                  var_name, coef_val, var_name, abs(coef_val)))
    } else {
      cat(sprintf("%s = %.6f: A one-unit increase in %s is associated with a %.6f unit DECREASE in CO(GT), on average.\n", 
                  var_name, coef_val, var_name, abs(coef_val)))
    }
  }
}

#11
cat("========== Residual Systematic Errors Test ==========\n\n")

# Breusch-Pagan test for heteroscedasticity
if (!require(lmtest)) {
  install.packages("lmtest")
  library(lmtest)
}
bp_test <- bptest(final_cleaned_model)
cat(sprintf("Breusch-Pagan test p-value: %.6e\n", bp_test$p.value))
if (bp_test$p.value < 0.05) {
  cat("Conclusion: Systematic error detected (heteroscedasticity).\n")
} else {
  cat("Conclusion: No systematic error detected.\n")
}

#12
plot(model_full, which = 1)

#13
if (!require(lmtest)) {
  install.packages("lmtest")
  library(lmtest)
}

bptest(model_full)


#14
plot(residuals(model_full), type = "l", 
     xlab = "Observation", ylab = "Residuals", 
     main = "Residuals vs Order")
abline(h = 0, col = "red")

#15
if (!require(lmtest)) {
  install.packages("lmtest")
  library(lmtest)
}

dwtest(model_full)

bgtest(model_full, order = 1)

#16
shapiro.test(residuals(model_full))

#17
cooks_d <- cooks.distance(model_full)
outliers <- which(cooks_d > 4/length(cooks_d))
if (length(outliers) > 0) {
  AirQuality_no_outliers <- AirQuality[-outliers, ]
  model_no_outliers <- lm(`CO(GT)` ~ ., data = AirQuality_no_outliers)
}

# Robust standard errors (Heteroscedasticity-consistent)
if (!require(sandwich)) {
  install.packages("sandwich")
  library(sandwich)
}
if (!require(lmtest)) {
  install.packages("lmtest")
  library(lmtest)
}
coeftest(model_full, vcov = vcovHC(model_full))

# Weighted Least Squares
weights <- 1 / (abs(residuals(model_full)) + 0.1)
model_wls <- lm(`CO(GT)` ~ ., data = AirQuality, weights = weights)

#18
means <- colMeans(AirQuality[, names(AirQuality) != "CO(GT)"], na.rm = TRUE)

new_data_10 <- as.data.frame(t(means * 1.1))
names(new_data_10) <- names(means)

new_data_20 <- as.data.frame(t(means * 1.2))
names(new_data_20) <- names(means)

predict_10 <- predict(model_full, new_data_10, interval = "confidence", level = 0.95)
predict_20 <- predict(model_full, new_data_20, interval = "confidence", level = 0.95)

print("Forecast for 10% increase:")
print(predict_10)

print("Forecast for 20% increase:")
print(predict_20)





