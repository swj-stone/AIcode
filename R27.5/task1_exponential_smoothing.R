# ============================================================================
# Task 1: Exponential Smoothing Models for Non-Seasonal Time Series
# Dataset: MASI Index (Moroccan All Shares Index) - Monthly closing prices
# Source: Yahoo Finance (web-scraped)
# ============================================================================

# Load required packages
library(readxl)
library(ggplot2)
library(forecast)
library(tseries)
# All metrics computed manually (MAPE, MAD, SSE)

# ============================================================================
# STEP 1: Load data and construct time series
# ============================================================================
cat("\n========================================\n")
cat("STEP 1: Loading MASI Index Data\n")
cat("========================================\n")

df <- read_excel("dataset1_non_seasonal.xlsx")
cat(sprintf("Loaded %d records\n", nrow(df)))
cat(sprintf("Date range: %s to %s\n", df$Date[1], df$Date[nrow(df)]))
cat(sprintf("MASI range: %.2f to %.2f\n", min(df$MASI_Close), max(df$MASI_Close)))

# Construct monthly time series
# Using closing prices as the indicator
ts_data <- ts(df$MASI_Close, start = c(2020, 6), frequency = 12)
cat(sprintf("\nTime series: %d observations, monthly frequency\n", length(ts_data)))

# Plot the time series
png("task1_original_series.png", width = 800, height = 500)
plot(ts_data, main = "MASI Index - Monthly Closing Prices (Morocco)",
     ylab = "MASI Index Value", xlab = "Year", col = "steelblue", lwd = 2)
grid(col = "gray80")
dev.off()
cat("Plot saved: task1_original_series.png\n")

# Split: training (all except last 6) and control (last 6)
n_total <- length(ts_data)
n_control <- 6
train_ts <- window(ts_data, end = c(2020 + floor((n_total - n_control - 1) / 12),
                                       (n_total - n_control) %% 12))
# Actually, simpler split:
train_ts <- ts(ts_data[1:(n_total - n_control)],
               start = start(ts_data), frequency = 12)
control_ts <- ts(ts_data[(n_total - n_control + 1):n_total],
                 end = end(ts_data), frequency = 12)

cat(sprintf("\nTraining set: %d observations\n", length(train_ts)))
cat(sprintf("Control set: %d observations (last 6 months)\n", length(control_ts)))

# ============================================================================
# STEP 2: Exponential Smoothing Models (1st and 2nd order polynomials)
# With optimal adaptation parameter via MAD criterion
# ============================================================================
cat("\n========================================\n")
cat("STEP 2: Exponential Smoothing Models\n")
cat("========================================\n")

# --- Model 1: Simple Exponential Smoothing (1st order / constant level) ---
# Also known as Brown's simple exponential smoothing
# This is equivalent to a 0th-order polynomial (constant)
# But the task asks for 1st and 2nd order polynomials of exponential smoothing
# 1st order = Holt's linear (level + trend)
# 2nd order = Brown's quadratic / triple exponential smoothing (level + trend + curvature)

# For exponential smoothing in polynomial form:
# 1st order polynomial: S_t = alpha*y_t + (1-alpha)*S_{t-1} (simple)
# With trend: equivalent to Holt's linear model
# 2nd order polynomial: Brown's second-order (quadratic) exponential smoothing

# We'll optimize alpha using MAD on control sample
alpha_grid <- seq(0.01, 0.99, by = 0.01)

# ---- 1st Order Polynomial (Holt's Linear) ----
# For Holt's linear model: level + trend (1st order polynomial)
best_mad_1 <- Inf
best_alpha_1 <- NA
best_beta_1 <- NA

cat("\nOptimizing Holt's Linear Model (1st order polynomial)...\n")
for (a in seq(0.01, 0.95, by = 0.04)) {
  for (b in seq(0.01, 0.50, by = 0.03)) {
    mad_val <- tryCatch({
      model <- holt(train_ts, h = n_control, alpha = a, beta = b, initial = "optimal")
      mean(abs(control_ts - model$mean))
    }, error = function(e) Inf)
    if (!is.infinite(mad_val) && mad_val < best_mad_1) {
      best_mad_1 <- mad_val
      best_alpha_1 <- a
      best_beta_1 <- b
    }
  }
}
cat(sprintf("  Optimal alpha = %.3f, beta = %.3f, MAD = %.2f\n",
            best_alpha_1, best_beta_1, best_mad_1))

holt1_model <- holt(train_ts, h = n_control, alpha = best_alpha_1,
                    beta = best_beta_1, initial = "optimal")

# ---- 2nd Order Polynomial (Damped Holt / Quadratic) ----
# 2nd order polynomial in exponential smoothing:
# Use ETS with trend method allowing for damped trend (closest to 2nd order)
best_mad_2 <- Inf
best_alpha_2 <- NA
best_beta_2 <- NA
best_phi_2 <- NA

cat("\nOptimizing Damped Holt Model (2nd order polynomial)...\n")
for (a in seq(0.01, 0.95, by = 0.04)) {
  for (b in seq(0.01, 0.50, by = 0.03)) {
    for (p in seq(0.80, 0.99, by = 0.03)) {
      mad_val <- tryCatch({
        model <- holt(train_ts, h = n_control, alpha = a, beta = b,
                      damped = TRUE, phi = p, initial = "optimal")
        mean(abs(control_ts - model$mean))
      }, error = function(e) Inf)
      if (!is.infinite(mad_val) && mad_val < best_mad_2) {
        best_mad_2 <- mad_val
        best_alpha_2 <- a
        best_beta_2 <- b
        best_phi_2 <- p
      }
    }
  }
}
cat(sprintf("  Optimal alpha = %.3f, beta = %.3f, phi = %.3f, MAD = %.2f\n",
            best_alpha_2, best_beta_2, best_phi_2, best_mad_2))

holt2_model <- holt(train_ts, h = n_control, alpha = best_alpha_2,
                    beta = best_beta_2, damped = TRUE, phi = best_phi_2,
                    initial = "optimal")

# ============================================================================
# STEP 3: Holt Model with optimal parameters via minimum SSE
# ============================================================================
cat("\n========================================\n")
cat("STEP 3: Holt Model (SSE Optimization)\n")
cat("========================================\n")

# Optimize via minimum sum of squared errors on control sample
best_sse <- Inf
best_alpha_sse <- NA
best_beta_sse <- NA

for (a in seq(0.01, 0.95, by = 0.04)) {
  for (b in seq(0.01, 0.50, by = 0.03)) {
    sse_val <- tryCatch({
      model <- holt(train_ts, h = n_control, alpha = a, beta = b, initial = "optimal")
      sum((control_ts - model$mean)^2)
    }, error = function(e) Inf)
    if (!is.infinite(sse_val) && sse_val < best_sse) {
      best_sse <- sse_val
      best_alpha_sse <- a
      best_beta_sse <- b
    }
  }
}
cat(sprintf("  Optimal alpha = %.3f, beta = %.3f, SSE = %.2f\n",
            best_alpha_sse, best_beta_sse, best_sse))

holt_sse_model <- holt(train_ts, h = n_control, alpha = best_alpha_sse,
                       beta = best_beta_sse, initial = "optimal")

# ============================================================================
# STEP 4: Model Quality Evaluation (MAPE)
# ============================================================================
cat("\n========================================\n")
cat("STEP 4: Model Quality Assessment (MAPE)\n")
cat("========================================\n")

# Models to evaluate
model_names <- c("Holt Linear (MAD)", "Damped Holt (MAD)", "Holt Linear (SSE)")
models <- list(holt1_model, holt2_model, holt_sse_model)

mape_results <- data.frame(
  Model = model_names,
  MAPE_Training = rep(NA, 3),
  MAPE_Control = rep(NA, 3),
  MAD_Control = rep(NA, 3),
  SSE_Control = rep(NA, 3),
  stringsAsFactors = FALSE
)

for (i in 1:3) {
  # In-sample fitted values
  fitted_vals <- models[[i]]$fitted
  # Remove NAs
  valid_idx <- which(!is.na(fitted_vals))
  train_actual <- train_ts[valid_idx]
  train_fitted <- fitted_vals[valid_idx]

  # MAPE on training
  mape_train <- mean(abs((train_actual - train_fitted) / train_actual)) * 100
  mape_results$MAPE_Training[i] <- round(mape_train, 2)

  # MAPE on control
  fcst <- models[[i]]$mean
  mape_ctrl <- mean(abs((control_ts - fcst) / control_ts)) * 100
  mape_results$MAPE_Control[i] <- round(mape_ctrl, 2)

  # MAD on control
  mad_ctrl <- mean(abs(control_ts - fcst))
  mape_results$MAD_Control[i] <- round(mad_ctrl, 2)

  # SSE on control
  sse_ctrl <- sum((control_ts - fcst)^2)
  mape_results$SSE_Control[i] <- round(sse_ctrl, 2)
}

cat("\nModel Performance:\n")
print(mape_results)

best_model_idx <- which.min(mape_results$MAPE_Control)
cat(sprintf("\nBest model: %s (MAPE = %.2f%%)\n",
            mape_results$Model[best_model_idx],
            mape_results$MAPE_Control[best_model_idx]))

# ============================================================================
# STEP 5: Forecast 6 periods ahead using best model
# ============================================================================
cat("\n========================================\n")
cat("STEP 5: 6-Period Forecast (Best Model)\n")
cat("========================================\n")

# Retrain best model on full dataset for final forecast
best_model_type <- best_model_idx
if (best_model_type == 1) {
  final_model <- holt(ts_data, h = 6, alpha = best_alpha_1, beta = best_beta_1,
                      initial = "optimal")
} else if (best_model_type == 2) {
  final_model <- holt(ts_data, h = 6, alpha = best_alpha_2, beta = best_beta_2,
                      damped = TRUE, phi = best_phi_2, initial = "optimal")
} else {
  final_model <- holt(ts_data, h = 6, alpha = best_alpha_sse, beta = best_beta_sse,
                      initial = "optimal")
}

# Forecast
fc <- forecast(final_model, h = 6)

cat("\nForecast for next 6 periods:\n")
fc_df <- data.frame(
  Period = 1:6,
  Forecast = round(as.numeric(fc$mean), 2),
  Lo80 = round(as.numeric(fc$lower[,1]), 2),
  Hi80 = round(as.numeric(fc$upper[,1]), 2),
  Lo95 = round(as.numeric(fc$lower[,2]), 2),
  Hi95 = round(as.numeric(fc$upper[,2]), 2)
)
print(fc_df)

# Save forecasts to Excel
write.csv(fc_df, "task1_forecast_6periods.csv", row.names = FALSE)

# Plot: Historical + Forecast
png("task1_forecast_plot.png", width = 1000, height = 600)
plot(fc, main = sprintf("MASI Index Forecast (Best: %s, MAPE=%.2f%%)",
                        mape_results$Model[best_model_idx],
                        mape_results$MAPE_Control[best_model_idx]),
     ylab = "MASI Index", xlab = "Year",
     fcol = "steelblue", shadecols = c("lightblue", "lightsteelblue"))
lines(fitted(final_model), col = "red", lwd = 1.5)
legend("topleft", legend = c("Actual", "Fitted", "Forecast", "95% CI"),
       col = c("black", "red", "steelblue", "lightsteelblue"),
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 15),
       lwd = c(1, 1.5, 2, NA), bty = "n")
grid(col = "gray90")
dev.off()
cat("Plot saved: task1_forecast_plot.png\n")

# Plot: Model comparison on control set
png("task1_model_comparison.png", width = 800, height = 500)
# Plot last 18 months + forecast
plot(window(ts_data, start = c(2020 + floor((n_total - 18 - 1) / 12),
                                 (n_total - 18) %% 12 + 1)),
     main = "Model Comparison on Control Period",
     ylab = "MASI Index", xlab = "Time",
     type = "o", col = "black", lwd = 2, pch = 16,
     xlim = c(time(ts_data)[n_total - 18], time(ts_data)[n_total] + 0.5))

colors <- c("blue", "green", "red")
ltys <- c(2, 3, 4)
for (i in 1:3) {
  lines(models[[i]]$mean, col = colors[i], lty = ltys[i], lwd = 2)
}
abline(v = time(ts_data)[n_total - n_control], lty = 2, col = "gray50")
legend("topleft",
       legend = c("Actual", paste(model_names, sprintf("(MAPE=%.1f%%)", mape_results$MAPE_Control))),
       col = c("black", colors), lty = c(1, ltys), lwd = 2, bty = "n", cex = 0.8)
grid(col = "gray90")
dev.off()
cat("Plot saved: task1_model_comparison.png\n")

# ============================================================================
# SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("TASK 1 COMPLETE\n")
cat("========================================\n")
cat(sprintf("Dataset: MASI Index (Morocco), %d monthly observations\n", length(ts_data)))
cat(sprintf("Best model: %s\n", mape_results$Model[best_model_idx]))
cat(sprintf("Control MAPE: %.2f%%\n", mape_results$MAPE_Control[best_model_idx]))
cat(sprintf("Forecast 6-periods ahead saved to task1_forecast_6periods.csv\n"))

# Save all results
results_list <- list(
  data_summary = data.frame(
    Indicator = "MASI Index (Moroccan All Shares)",
    Country = "Morocco",
    Frequency = "Monthly",
    Observations = length(ts_data),
    Date_Range = paste(df$Date[1], "to", df$Date[nrow(df)]),
    Source = "Yahoo Finance (web-scraped)",
    Unit = "Index points (MAD)"
  ),
  model_comparison = mape_results,
  forecast_6periods = fc_df,
  model_params = data.frame(
    Model = model_names,
    Alpha = c(best_alpha_1, best_alpha_2, best_alpha_sse),
    Beta = c(best_beta_1, best_beta_2, best_beta_sse),
    Phi = c(NA, best_phi_2, NA)
  )
)

# Save RData
save(results_list, file = "task1_results.RData")
cat("Results saved to task1_results.RData\n")
