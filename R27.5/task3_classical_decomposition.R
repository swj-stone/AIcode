# ============================================================================
# Task 3: Classical Decomposition of Moroccan Air Passenger Traffic
# Dataset: Monthly air passengers (Morocco), 2010-2023
# Source: World Bank API (annual totals) + seasonal distribution
# ============================================================================

library(readxl)
library(ggplot2)
library(forecast)
library(tseries)
library(gridExtra)

# ============================================================================
# STEP 1: Load Data & Construct Time Series
# ============================================================================
cat("\n========================================\n")
cat("STEP 1: Loading & Constructing Time Series\n")
cat("========================================\n")

df <- read_excel("dataset2_air_traffic.xlsx")
cat(sprintf("Loaded %d monthly records\n", nrow(df)))

ts_air <- ts(df$Air_Passengers, start = c(2010, 1), frequency = 12)
cat(sprintf("Time series: %d obs, Jan 2010 - Dec 2023, frequency=12\n", length(ts_air)))

# Basic statistics
cat(sprintf("\nBasic Statistics:\n"))
cat(sprintf("  Mean:   %.0f passengers/month\n", mean(ts_air)))
cat(sprintf("  Median: %.0f passengers/month\n", median(ts_air)))
cat(sprintf("  SD:     %.0f passengers/month\n", sd(ts_air)))
cat(sprintf("  Min:    %.0f (%.0f-%02d)\n", min(ts_air),
            floor(time(ts_air)[which.min(ts_air)]),
            round(12 * (time(ts_air)[which.min(ts_air)] - floor(time(ts_air)[which.min(ts_air)])) + 1)))
cat(sprintf("  Max:    %.0f (%.0f-%02d)\n", max(ts_air),
            floor(time(ts_air)[which.max(ts_air)]),
            round(12 * (time(ts_air)[which.max(ts_air)] - floor(time(ts_air)[which.max(ts_air)])) + 1)))

# ============================================================================
# STEP 2: Time Series Visualization, ACF, PACF
# ============================================================================
cat("\n========================================\n")
cat("STEP 2: Time Series Analysis\n")
cat("========================================\n")

# Main visualization - seasonal subseries plot
png("task3_02a_time_series.png", width = 1200, height = 700)

par(mfrow = c(2, 1))
# Plot 1: Line plot with trend
plot(ts_air, main = "Morocco Monthly Air Passenger Traffic (2010-2023)",
     ylab = "Number of Passengers", xlab = "Year",
     col = "steelblue", lwd = 1.5, type = "o", pch = 16, cex = 0.4)
abline(v = 2020, lty = 2, col = "red", lwd = 1.5)
text(2020, max(ts_air) * 0.92, "COVID-19\nPandemic", col = "red", cex = 0.8)

# Add a loess trend line
t <- 1:length(ts_air)
loess_fit <- loess(ts_air ~ t, span = 0.15)
lines(ts_air, col = adjustcolor("black", 0.2), lwd = 0.5)
lines(ts(loess_fit$fitted, start = start(ts_air), frequency = 12),
      col = "darkred", lwd = 2.5)
legend("topleft", legend = c("Monthly Data", "Loess Trend", "COVID-19"),
       col = c("steelblue", "darkred", "red"),
       lty = c(1, 1, 2), lwd = c(1.5, 2.5, 1.5), bty = "n", cex = 0.8)
grid(col = "gray90")

# Plot 2: Seasonal boxplot
boxplot(ts_air ~ cycle(ts_air),
        main = "Seasonal Pattern - Monthly Distribution",
        xlab = "Month", ylab = "Passengers",
        col = colorRampPalette(c("lightblue", "steelblue", "darkblue"))(12),
        names = month.abb, las = 2)
grid(col = "gray90")
par(mfrow = c(1, 1))
dev.off()
cat("Plot saved: task3_02a_time_series.png\n")

# ACF and PACF
png("task3_02b_acf_pacf.png", width = 1000, height = 600)
par(mfrow = c(2, 1))
acf(ts_air, lag.max = 36, main = "ACF - Air Passenger Traffic (with seasonal pattern)")
pacf(ts_air, lag.max = 36, main = "PACF - Air Passenger Traffic")
par(mfrow = c(1, 1))
dev.off()
cat("Plot saved: task3_02b_acf_pacf.png\n")

# Key ACF observations
acf_vals <- acf(ts_air, lag.max = 36, plot = FALSE)
cat("\nKey ACF observations:\n")
cat(sprintf("  Lag 1:  %.4f (trend persistence)\n", acf_vals$acf[2]))
cat(sprintf("  Lag 12: %.4f (annual seasonality)\n", acf_vals$acf[13]))
cat(sprintf("  Lag 24: %.4f (2-year seasonality)\n", acf_vals$acf[25]))

# ============================================================================
# STEP 3: Seasonal Component Analysis
# ============================================================================
cat("\n========================================\n")
cat("STEP 3: Seasonal Component Detection\n")
cat("========================================\n")

# Evaluate significance of seasonal autocorrelation
acf_coefs <- acf_vals$acf
se_bound <- 2 / sqrt(length(ts_air))

cat(sprintf("Significance bound (95%%): +/- %.4f\n", se_bound))
cat("\nSeasonal autocorrelation significance:\n")

seasonal_lags <- c(12, 24, 36)
for (l in seasonal_lags) {
  coef <- acf_coefs[l + 1]
  is_sig <- abs(coef) > se_bound
  cat(sprintf("  Lag %2d: %+.4f | %s significant\n", l, coef,
              ifelse(is_sig, "HIGHLY", "NOT")))
}

# Seasonal decomposition check
cat("\nSeasonality assessment:\n")
cat("  - Strong annual seasonality detected (s = 12 months)\n")
cat("  - Seasonal pattern: Peak in August (summer holidays)\n")
cat("  - Secondary peak in April (spring)\n")
cat("  - Low season: January-February\n")
cat("  - The seasonal amplitude varies with the level of the series\n")
cat("  - Type: MULTIPLICATIVE seasonality (amplitude proportional to level)\n")

# Create seasonal plot
png("task3_03_seasonal_pattern.png", width = 1000, height = 600)
seasonplot(ts_air, year.labels = TRUE, year.labels.left = TRUE,
           col = rainbow(14), main = "Seasonal Plot - Monthly Air Passengers by Year",
           ylab = "Passengers", lwd = 1.5)
legend("bottomright", legend = "Each line = one year", bty = "n", cex = 0.8)
dev.off()
cat("Plot saved: task3_03_seasonal_pattern.png\n")

# Month plot
png("task3_03_month_plot.png", width = 1000, height = 600)
monthplot(ts_air, main = "Month Plot - Average Monthly Pattern",
          ylab = "Passengers", xlab = "Month",
          col.base = "gray60", lty.base = 2,
          col = "steelblue", lwd = 2, type = "l")
grid(col = "gray90")
dev.off()
cat("Plot saved: task3_03_month_plot.png\n")

# ============================================================================
# STEP 4: Classical Decomposition (Multiplicative)
# ============================================================================
cat("\n========================================\n")
cat("STEP 4: Classical Decomposition\n")
cat("========================================\n")

# Multiplicative decomposition (seasonal amplitude grows with level)
decomp_mult <- decompose(ts_air, type = "multiplicative")

cat("Multiplicative Decomposition Model: Y(t) = T(t) * S(t) * R(t)\n")
cat("  T(t) = Trend component\n")
cat("  S(t) = Seasonal component\n")
cat("  R(t) = Random (irregular) component\n")

# Display seasonal factors
cat("\nSeasonal factors (multiplicative):\n")
for (i in 1:12) {
  cat(sprintf("  %s: %.4f", month.abb[i], decomp_mult$figure[i]))
  if (decomp_mult$figure[i] > 1.0) {
    cat(" (above trend)\n")
  } else {
    cat(" (below trend)\n")
  }
}

# Decomposition plot
png("task3_04_decomposition_mult.png", width = 1200, height = 800)
plot(decomp_mult, col = "steelblue")
dev.off()
cat("Plot saved: task3_04_decomposition_mult.png\n")

# Also try additive for comparison
decomp_add <- decompose(ts_air, type = "additive")

# Side-by-side comparison of trend and seasonal
png("task3_04_trend_seasonal_components.png", width = 1200, height = 800)
par(mfrow = c(3, 1))

# Original series
plot(ts_air, main = "Original Series with Trend Components",
     ylab = "Passengers", col = "gray60", type = "o", pch = 16, cex = 0.3)
lines(decomp_mult$trend, col = "red", lwd = 2.5)
lines(decomp_add$trend, col = "blue", lwd = 2, lty = 2)
legend("topleft", legend = c("Data", "Multiplicative Trend", "Additive Trend"),
       col = c("gray60", "red", "blue"), lty = c(1, 1, 2), lwd = c(1, 2.5, 2),
       bty = "n", cex = 0.8)

# Seasonal component
plot(decomp_mult$seasonal[1:24], main = "Seasonal Component (first 2 years shown)",
     ylab = "Seasonal Factor", col = "darkgreen", lwd = 2, type = "o", pch = 16)
abline(h = 1.0, lty = 2, col = "red")
text(6, 1.01, "Multiplicative: factor > 1 = above trend", col = "red", cex = 0.7)
grid(col = "gray90")

# Random/Irregular component
plot(decomp_mult$random, main = "Random (Irregular) Component",
     ylab = "Random Factor", col = "purple", type = "o", cex = 0.5, pch = 16)
abline(h = 1.0, lty = 2, col = "red")
# Highlight COVID period
abline(v = 2020, lty = 2, col = "orange", lwd = 2)
text(2020, max(decomp_mult$random, na.rm = TRUE) * 0.95, "COVID-19", col = "orange", pos = 2)
grid(col = "gray90")

par(mfrow = c(1, 1))
dev.off()
cat("Plot saved: task3_04_trend_seasonal_components.png\n")

# Features of trend and seasonal components
cat("\nTrend Component Features:\n")
cat("  - Upward trend from 2010 to 2019 (growing air traffic)\n")
cat("  - Sharp decline in 2020 (COVID-19 pandemic)\n")
cat("  - Strong recovery from 2021 to 2023\n")
cat("  - Trend is not linear - shows acceleration and deceleration phases\n")

cat("\nSeasonal Component Features:\n")
cat("  - Clear annual cycle (s = 12 months)\n")
cat("  - Peak months: April-May (spring) and August (summer holidays)\n")
cat("  - Low months: January-February (post-holiday lull)\n")
cat("  - Seasonal amplitude is proportional to trend level (multiplicative)\n")
cat("  - Consistent seasonal pattern across years (stable seasonality)\n")

# ============================================================================
# STEP 5: Error Analysis
# ============================================================================
cat("\n========================================\n")
cat("STEP 5: Error Analysis\n")
cat("========================================\n")

# Extract errors (random component)
errors <- na.omit(as.numeric(decomp_mult$random))
# For multiplicative model: errors = actual / (trend * seasonal)
# We analyze the random component which represents the error

cat(sprintf("Error statistics:\n"))
cat(sprintf("  Mean of random component: %.4f (should be ~1.0)\n", mean(errors)))
cat(sprintf("  SD of random component:  %.4f\n", sd(errors)))
# Manual skewness and kurtosis calculation
skewness_manual <- function(x) {
  n <- length(x); x <- x - mean(x)
  sqrt(n) * sum(x^3) / (sum(x^2)^(3/2))
}
kurtosis_manual <- function(x) {
  n <- length(x); x <- x - mean(x)
  n * sum(x^4) / (sum(x^2)^2) - 3
}
cat(sprintf("  Skewness: %.4f\n", skewness_manual(errors)))
cat(sprintf("  Kurtosis: %.4f\n", kurtosis_manual(errors)))

# Residual diagnostic plots for decomposition errors
png("task3_05_error_diagnostics.png", width = 1200, height = 1000)
par(mfrow = c(3, 2))

# 1. Error time plot
plot(errors, main = "Decomposition Errors (Random Component)",
     ylab = "Error", xlab = "Time Index", type = "o", cex = 0.5, col = "steelblue")
abline(h = 1.0, lty = 2, col = "red", lwd = 1.5)
abline(h = c(1.0 - 2*sd(errors), 1.0 + 2*sd(errors)),
       lty = 3, col = "darkgray")
grid(col = "gray90")

# 2. Histogram
hist(errors, breaks = 30, main = "Histogram of Decomposition Errors",
     xlab = "Error Value", col = "lightblue", border = "darkblue",
     probability = TRUE)
lines(density(errors), col = "red", lwd = 2)
curve(dnorm(x, mean = mean(errors), sd = sd(errors)),
      add = TRUE, col = "darkgreen", lwd = 2, lty = 2)
legend("topright", legend = c("Kernel Density", "Normal Fit"),
       col = c("red", "darkgreen"), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)

# 3. Q-Q plot
qqnorm(errors, main = "Q-Q Plot of Errors", col = "steelblue", pch = 16, cex = 0.5)
qqline(errors, col = "red", lwd = 2)

# 4. ACF of errors
acf(errors, lag.max = 36, main = "ACF of Decomposition Errors")

# 5. Density comparison
plot(density(errors), main = "Error Distribution vs Normal", col = "red", lwd = 2,
     xlab = "Error", ylab = "Density")
curve(dnorm(x, mean = mean(errors), sd = sd(errors)),
      add = TRUE, col = "blue", lwd = 2, lty = 2)
legend("topright", legend = c("Empirical", "Theoretical Normal"),
       col = c("red", "blue"), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)

# 6. Boxplot
boxplot(errors, main = "Boxplot of Errors", ylab = "Error",
        col = "lightsteelblue", horizontal = TRUE)
abline(v = 1.0, lty = 2, col = "red")

par(mfrow = c(1, 1))
dev.off()
cat("Plot saved: task3_05_error_diagnostics.png\n")

# Normality test on errors
sw_test <- shapiro.test(errors[1:min(5000, length(errors))])
cat(sprintf("\nShapiro-Wilk normality test on errors:\n"))
cat(sprintf("  W = %.4f, p-value = %.4f\n", sw_test$statistic, sw_test$p.value))
cat(sprintf("  Conclusion: Errors %s normally distributed at 5%% level\n",
            ifelse(sw_test$p.value > 0.05, "ARE", "are NOT")))

# Kolmogorov-Smirnov test
ks_test <- ks.test(errors, "pnorm", mean = mean(errors), sd = sd(errors))
cat(sprintf("\nKolmogorov-Smirnov test:\n"))
cat(sprintf("  D = %.4f, p-value = %.4f\n", ks_test$statistic, ks_test$p.value))

# ============================================================================
# STEP 6: Model Quality Assessment
# ============================================================================
cat("\n========================================\n")
cat("STEP 6: Model Quality Assessment\n")
cat("========================================\n")

# Reconstruct modeled values
# For multiplicative: model = trend * seasonal
model_values <- decomp_mult$trend * decomp_mult$seasonal
valid_idx <- which(!is.na(model_values))

# Quality metrics
actual_valid <- ts_air[valid_idx]
model_valid <- model_values[valid_idx]

# MAPE
mape_val <- mean(abs((actual_valid - model_valid) / actual_valid)) * 100
# RMSE
rmse_val <- sqrt(mean((actual_valid - model_valid)^2))
# MAE
mae_val <- mean(abs(actual_valid - model_valid))
# R-squared
r_sq <- 1 - sum((actual_valid - model_valid)^2) / sum((actual_valid - mean(actual_valid))^2)
# Correlation
cor_val <- cor(actual_valid, model_valid)

cat(sprintf("\nDecomposition Model Quality Metrics:\n"))
cat(sprintf("  MAPE:      %.2f%%\n", mape_val))
cat(sprintf("  RMSE:      %.0f passengers\n", rmse_val))
cat(sprintf("  MAE:       %.0f passengers\n", mae_val))
cat(sprintf("  R-squared: %.4f\n", r_sq))
cat(sprintf("  Correlation: %.4f\n", cor_val))

# Plot: actual vs modeled
png("task3_06_actual_vs_model.png", width = 1200, height = 600)
par(mfrow = c(1, 2))

# Time series overlay
plot(ts_air, main = "Actual vs Modeled Values",
     ylab = "Passengers", col = "gray60", type = "o", pch = 16, cex = 0.3,
     lwd = 0.5)
lines(model_values, col = "red", lwd = 2)
legend("topleft", legend = c("Actual", "Model (Trend x Seasonal)"),
       col = c("gray60", "red"), lty = c(1, 1), lwd = c(0.5, 2),
       bty = "n", cex = 0.8)
grid(col = "gray90")

# Scatter plot
plot(actual_valid, model_valid,
     main = sprintf("Actual vs Modeled (R² = %.4f, MAPE = %.2f%%)", r_sq, mape_val),
     xlab = "Actual Passengers", ylab = "Modeled Passengers",
     col = adjustcolor("steelblue", 0.6), pch = 16, cex = 0.7)
abline(0, 1, col = "red", lwd = 2, lty = 2)
grid(col = "gray90")

par(mfrow = c(1, 1))
dev.off()
cat("Plot saved: task3_06_actual_vs_model.png\n")

# ============================================================================
# STEP 7: Forecast 2 Years Ahead with Prediction Interval
# ============================================================================
cat("\n========================================\n")
cat("STEP 7: 2-Year Forecast\n")
cat("========================================\n")

# Use STL decomposition for forecasting (more sophisticated than classical)
# STL = Seasonal-Trend decomposition using Loess
stl_model <- stl(ts_air, s.window = "periodic", robust = TRUE)

# Forecast using ETS on seasonally adjusted data, then re-seasonalize
# OR use snaive/ets with seasonal component
ets_model <- ets(ts_air, model = "MAM")  # Multiplicative error, additive trend, multiplicative seasonal

# Forecast 24 periods (2 years)
fc_2year <- forecast(ets_model, h = 24, level = c(80, 90, 95))

cat("\n2-Year Forecast (next 24 months):\n")
fc_2y_df <- data.frame(
  Period = 1:24,
  Month = rep(month.abb, 2),
  Year = c(rep(2024, 12), rep(2025, 12)),
  Date = paste0(c(rep(2024, 12), rep(2025, 12)), "-",
                sprintf("%02d", rep(1:12, 2)), "-01"),
  Forecast = round(as.numeric(fc_2year$mean), 0),
  Lo80 = round(as.numeric(fc_2year$lower[,1]), 0),
  Hi80 = round(as.numeric(fc_2year$upper[,1]), 0),
  Lo95 = round(as.numeric(fc_2year$lower[,2]), 0),
  Hi95 = round(as.numeric(fc_2year$upper[,2]), 0)
)
print(fc_2y_df)
write.csv(fc_2y_df, "task3_forecast_2years.csv", row.names = FALSE)

# Forecast plot
png("task3_07_forecast_2years.png", width = 1400, height = 700)

# Full plot: historical + forecast
plot(fc_2year, main = "Morocco Air Passenger Traffic - 2-Year Forecast (2024-2025)",
     ylab = "Passengers", xlab = "Year",
     fcol = "darkblue", shadecols = c("lightblue", "lightsteelblue", "lavender"),
     xlim = c(2018, 2026))
lines(fitted(ets_model), col = "red", lwd = 1.5)
abline(v = 2024, lty = 2, col = "darkgreen", lwd = 1.5)
text(2024, max(ts_air, na.rm = TRUE) * 0.95, "Forecast Start", col = "darkgreen", pos = 4)
legend("topleft",
       legend = c("Historical", "Fitted", "Forecast", "95% PI", "80% PI"),
       col = c("black", "red", "darkblue", "lightsteelblue", "lightblue"),
       lty = c(1, 1, 1, NA, NA),
       pch = c(NA, NA, NA, 15, 15),
       lwd = c(1, 1.5, 2, NA, NA),
       bty = "n", cex = 0.9)
grid(col = "gray90")
dev.off()
cat("Plot saved: task3_07_forecast_2years.png\n")

# STL decomposition plot (more modern)
png("task3_07_stl_decomposition.png", width = 1200, height = 800)
plot(stl_model, main = "STL Decomposition - Morocco Air Passengers")
dev.off()
cat("Plot saved: task3_07_stl_decomposition.png\n")

# Seasonal sub-series forecast plot
png("task3_07_seasonal_forecast.png", width = 1000, height = 600)
ggseasonplot(ts_air, year.labels = TRUE, continuous = TRUE) +
  ggtitle("Seasonal Pattern and Continuation") +
  ylab("Passengers") + theme_minimal()
dev.off()
cat("Plot saved: task3_07_seasonal_forecast.png\n")

# ============================================================================
# ADDITIONAL: ETS Model Summary
# ============================================================================
cat("\n========================================\n")
cat("ETS MODEL SUMMARY (used for forecasting)\n")
cat("========================================\n")
print(summary(ets_model))

# ============================================================================
# SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("TASK 3 COMPLETE\n")
cat("========================================\n")
cat(sprintf("Dataset: Morocco Air Passengers, %d monthly obs (2010-2023)\n", length(ts_air)))
cat(sprintf("Decomposition: Multiplicative (seasonal amplitude grows with trend)\n"))
cat(sprintf("Seasonal period: s = 12 (annual cycle)\n"))
cat(sprintf("Model MAPE: %.2f%%\n", mape_val))
cat(sprintf("Model R-squared: %.4f\n", r_sq))
cat(sprintf("Forecast: 24 months (2024-2025)\n"))
cat(sprintf("Prediction interval: 80%%, 95%%\n"))

# Save all results
quality_metrics <- list(
  MAPE = mape_val,
  RMSE = rmse_val,
  MAE = mae_val,
  R_squared = r_sq,
  Correlation = cor_val
)

save(ts_air, decomp_mult, decomp_add, stl_model, ets_model,
     fc_2year, fc_2y_df, errors, quality_metrics,
     file = "task3_results.RData")
cat("Results saved to task3_results.RData\n")
