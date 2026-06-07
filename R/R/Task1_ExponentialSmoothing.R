###############################################################################
# Practical Work №2 - Task 1: Exponential Smoothing & Holt Model
# Dataset: S&P 500 Monthly Close (non-seasonal)
###############################################################################

# ---- Load libraries ----
library(readxl)
library(ggplot2)
library(forecast)
library(tseries)
library(Metrics)

# ---- Load data ----
raw <- read_excel("Dataset1_SP500_monthly.xlsx")
names(raw) <- c("Month", "SP500")
cat("Loaded", nrow(raw), "observations from", as.character(raw$Month[1]),
    "to", as.character(raw$Month[nrow(raw)]), "\n")

# Create time series object (monthly frequency, start from Jan 2010)
sp500_ts <- ts(raw$SP500, start = c(2010, 1), frequency = 12)
cat("\nTime series summary:\n")
print(summary(sp500_ts))

# ---- Plot the series ----
png("Task1_plots/01_original_series.png", width = 800, height = 500)
plot(sp500_ts, main = "S&P 500 Monthly Close (2010-2024)",
     ylab = "Index Value (USD)", xlab = "Year", col = "steelblue", lwd = 2)
grid(col = "gray80")
dev.off()

# ---- Split into training and control sample ----
# Control sample: last 12 values (1 year)
n_total  <- length(sp500_ts)
n_control <- 12
train_ts  <- window(sp500_ts, end = c(2024 - (n_control/12), 12 - n_control))
control_ts <- window(sp500_ts, start = c(2024, 1))

cat("\nTraining set:", start(train_ts), "-", end(train_ts), "(", length(train_ts), "obs )\n")
cat("Control set:",  start(control_ts), "-", end(control_ts), "(", length(control_ts), "obs )\n")

###############################################################################
# 2. Exponential smoothing - 1st and 2nd order polynomials
###############################################################################

# ----- 2a. Simple Exponential Smoothing (1st order) -----
# Optimize alpha on training set by minimizing MAD on control set
alpha_grid <- seq(0.01, 0.99, by = 0.01)
mad_ses <- sapply(alpha_grid, function(a) {
  fit  <- ses(train_ts, alpha = a, h = n_control, initial = "optimal")
  fcst <- fit$mean
  mean(abs(control_ts - fcst))
})

opt_alpha_ses <- alpha_grid[which.min(mad_ses)]
cat("\n--- Simple Exponential Smoothing (1st order) ---\n")
cat(sprintf("Optimal alpha (min MAD on control): %.3f (MAD = %.2f)\n",
            opt_alpha_ses, min(mad_ses)))

# Fit final SES with optimal alpha
ses_final <- ses(train_ts, alpha = opt_alpha_ses, h = n_control, initial = "optimal")

# Calculate MAPE on control set
ses_fcst <- ses_final$mean
ses_mape <- 100 * mean(abs((control_ts - ses_fcst) / control_ts))
cat(sprintf("MAPE on control set: %.2f%%\n", ses_mape))

# ----- 2b. Second-order exponential smoothing (Holt's linear, but not damped) -----
# We test second-order (quadratic) smoothing via Brown's method
# SES with second smoothing (double exponential smoothing / Brown's method)

mad_brown <- sapply(alpha_grid, function(a) {
  fit  <- holt(train_ts, alpha = a, beta = a, h = n_control, initial = "optimal")
  fcst <- fit$mean
  mean(abs(control_ts - fcst))
})

opt_alpha_brown <- alpha_grid[which.min(mad_brown)]
cat("\n--- Double Exponential Smoothing (2nd order / Brown) ---\n")
cat(sprintf("Optimal alpha=beta (min MAD on control): %.3f (MAD = %.2f)\n",
            opt_alpha_brown, min(mad_brown)))

brown_final <- holt(train_ts, alpha = opt_alpha_brown, beta = opt_alpha_brown,
                     h = n_control, initial = "optimal")
brown_fcst <- brown_final$mean
brown_mape <- 100 * mean(abs((control_ts - brown_fcst) / control_ts))
cat(sprintf("MAPE on control set: %.2f%%\n", brown_mape))

###############################################################################
# 3. Holt model (1st order polynomial with separate alpha and beta)
###############################################################################

# Grid search for optimal alpha and beta using SSE on control set
alpha_vals <- seq(0.01, 0.99, by = 0.05)
beta_vals  <- seq(0.01, 0.99, by = 0.05)

best_sse <- Inf
best_alpha <- NA
best_beta  <- NA

for (a in alpha_vals) {
  for (b in beta_vals) {
    fit <- tryCatch({
      holt(train_ts, alpha = a, beta = b, h = n_control, initial = "optimal")
    }, error = function(e) NULL)
    if (!is.null(fit)) {
      sse <- sum((control_ts - fit$mean)^2)
      if (sse < best_sse) {
        best_sse   <- sse
        best_alpha <- a
        best_beta  <- b
      }
    }
  }
}

cat("\n--- Holt Model (1st order polynomial) ---\n")
cat(sprintf("Optimal alpha: %.3f (min SSE on control)\n", best_alpha))
cat(sprintf("Optimal beta:  %.3f (min SSE on control)\n", best_beta))
cat(sprintf("Minimum SSE on control set: %.2f\n", best_sse))

# Fit final Holt model
holt_final <- holt(train_ts, alpha = best_alpha, beta = best_beta,
                    h = n_control, initial = "optimal")
holt_fcst <- holt_final$mean
holt_mape <- 100 * mean(abs((control_ts - holt_fcst) / control_ts))
cat(sprintf("MAPE on control set: %.2f%%\n", holt_mape))

###############################################################################
# 4. Model quality comparison (MAPE)
###############################################################################

cat("\n", paste(rep("=", 55), collapse = ""), "\n", sep = "")
cat("MODEL COMPARISON (MAPE on control set)\n")
cat(paste(rep("=", 55), collapse = ""), "\n", sep = "")
cat(sprintf("  Simple Exp Smoothing (1st order, alpha=%.3f):  MAPE = %.2f%%\n",
            opt_alpha_ses, ses_mape))
cat(sprintf("  Double Exp Smoothing (2nd order, alpha=%.3f):  MAPE = %.2f%%\n",
            opt_alpha_brown, brown_mape))
cat(sprintf("  Holt Model (alpha=%.3f, beta=%.3f):           MAPE = %.2f%%\n",
            best_alpha, best_beta, holt_mape))

# Determine best model
models <- data.frame(
  Model = c("SES (1st order)", "Double ES (2nd order)", "Holt (1st order poly)"),
  Alpha = c(opt_alpha_ses, opt_alpha_brown, best_alpha),
  Beta  = c(NA, opt_alpha_brown, best_beta),
  MAPE  = c(ses_mape, brown_mape, holt_mape)
)

best_idx <- which.min(models$MAPE)
cat(sprintf("\n>>> BEST MODEL: %s with MAPE = %.2f%%\n",
            models$Model[best_idx], models$MAPE[best_idx]))

###############################################################################
# 5. Forecast 6 periods ahead using best model
###############################################################################

# Re-fit best model on ALL data and forecast 6 months
cat("\n--- Forecasting 6 periods ahead ---\n")

if (best_idx == 1) {
  best_model <- ses(sp500_ts, alpha = opt_alpha_ses, h = 6, initial = "optimal")
  model_name <- "Simple Exponential Smoothing"
} else if (best_idx == 2) {
  best_model <- holt(sp500_ts, alpha = opt_alpha_brown, beta = opt_alpha_brown,
                      h = 6, initial = "optimal")
  model_name <- "Double Exponential Smoothing (Brown)"
} else {
  best_model <- holt(sp500_ts, alpha = best_alpha, beta = best_beta,
                      h = 6, initial = "optimal")
  model_name <- "Holt Linear Trend"
}

cat(sprintf("Forecasting with: %s\n", model_name))
cat(sprintf("Model parameters: alpha = %.3f", best_model$model$par[1]))
if (best_idx > 1) cat(sprintf(", beta = %.3f", best_model$model$par[2]))
cat("\n\n6-Month Forecast (S&P 500 Close):\n")
fcst_df <- data.frame(
  Month = as.character(time(best_model$mean)),
  Forecast = as.numeric(best_model$mean),
  Lo80 = as.numeric(best_model$lower[, 1]),
  Hi80 = as.numeric(best_model$upper[, 1]),
  Lo95 = as.numeric(best_model$lower[, 2]),
  Hi95 = as.numeric(best_model$upper[, 2])
)
print(fcst_df, row.names = FALSE)

# ---- Plot final forecast ----
dir.create("Task1_plots", showWarnings = FALSE)
png("Task1_plots/02_forecast_6months.png", width = 900, height = 550)
plot(best_model, main = paste("S&P 500 6-Month Forecast -", model_name),
     ylab = "Index Value (USD)", xlab = "Year",
     fcol = "tomato", flwd = 2, shaded = TRUE)
grid(col = "gray80")
dev.off()

# ---- Save model details ----
sink("Task1_models_summary.txt")
cat("TASK 1 - MODEL SUMMARY\n")
cat(paste(rep("=", 55), collapse = ""), "\n\n")
cat(sprintf("Dataset: S&P 500 Monthly Close (%d obs)\n", n_total))
cat(sprintf("Training: %d obs, Control: %d obs\n\n", length(train_ts), n_control))

cat("1. Simple Exponential Smoothing (1st order):\n")
cat(sprintf("   Optimal alpha (MAD): %.3f\n", opt_alpha_ses))
cat(sprintf("   MAPE: %.2f%%\n\n", ses_mape))

cat("2. Double Exponential Smoothing (2nd order / Brown):\n")
cat(sprintf("   Optimal alpha=beta (MAD): %.3f\n", opt_alpha_brown))
cat(sprintf("   MAPE: %.2f%%\n\n", brown_mape))

cat("3. Holt Model (1st order polynomial):\n")
cat(sprintf("   Optimal alpha (SSE): %.3f\n", best_alpha))
cat(sprintf("   Optimal beta (SSE):  %.3f\n", best_beta))
cat(sprintf("   MAPE: %.2f%%\n\n", holt_mape))

cat("BEST MODEL:", model_name, "\n")
cat(sprintf("   MAPE: %.2f%%\n", models$MAPE[best_idx]))
cat("\nFORECAST (6 months ahead):\n")
print(fcst_df)
sink()

cat("\nTask 1 complete. Results saved to Task1_plots/ and Task1_models_summary.txt\n")
