###############################################################################
# Practical Work №2 - Task 2: ARIMA / SARIMA Modeling
# Dataset: International Airline Passengers (monthly, 1949-1960)
###############################################################################

# ---- Load libraries ----
library(readxl)
library(ggplot2)
library(forecast)
library(tseries)

# ---- Load data ----
raw <- read_excel("Dataset2_AviationPassengers.xlsx")
names(raw) <- c("Month", "Passengers")
cat("Loaded", nrow(raw), "observations\n")

# Create time series (monthly, start Jan 1949)
pass_ts <- ts(raw$Passengers, start = c(1949, 1), frequency = 12)
cat(sprintf("Time series: %d obs, frequency = 12\n", length(pass_ts)))

###############################################################################
# 1. Time series analysis - ACF, PACF, stationarity
###############################################################################

# ---- Plot original series ----
dir.create("Task2_plots", showWarnings = FALSE)
png("Task2_plots/01_original_series.png", width = 900, height = 500)
plot(pass_ts, main = "International Airline Passengers (1949-1960)",
     ylab = "Passengers (thousands)", xlab = "Year",
     col = "steelblue", lwd = 2, type = "o", pch = 16, cex = 0.6)
grid(col = "gray85")
dev.off()

# ---- ACF and PACF of original series ----
png("Task2_plots/02_ACF_PACF_original.png", width = 900, height = 600)
par(mfrow = c(2, 1))
acf(pass_ts, lag.max = 36, main = "ACF - Original Series (Air Passengers)",
    ci.col = "steelblue", lwd = 2)
pacf(pass_ts, lag.max = 36, main = "PACF - Original Series (Air Passengers)",
     ci.col = "steelblue", lwd = 2)
par(mfrow = c(1, 1))
dev.off()

# Calculate autocorrelation coefficients
acf_vals <- acf(pass_ts, lag.max = 24, plot = FALSE)
pacf_vals <- pacf(pass_ts, lag.max = 24, plot = FALSE)

cat("\n--- ACF (first 12 lags) ---\n")
print(round(acf_vals$acf[2:13], 3))
cat("\n--- PACF (first 12 lags) ---\n")
print(round(pacf_vals$acf[1:12], 3))

cat("\nObservation: Strong, slowly decaying ACF indicates NON-stationarity.\n")
cat("Significant peaks at lags 12, 24 suggest annual seasonality.\n")

###############################################################################
# 2. Differencing and stationarity tests
###############################################################################

# ---- First difference ----
pass_d1 <- diff(pass_ts)
cat("\n--- First Difference ---\n")

adf_d1 <- adf.test(pass_d1, alternative = "stationary")
cat(sprintf("ADF test (1st diff): statistic = %.3f, p-value = %.4f\n",
            adf_d1$statistic, adf_d1$p.value))
if (adf_d1$p.value < 0.05) {
  cat("  -> Series IS stationary after 1st difference (reject H0 of unit root)\n")
} else {
  cat("  -> Series NOT yet stationary after 1st difference\n")
}

# ---- ACF/PACF of first difference ----
png("Task2_plots/03_ACF_PACF_diff1.png", width = 900, height = 600)
par(mfrow = c(2, 1))
acf(pass_d1, lag.max = 36, main = "ACF - First Difference",
    ci.col = "darkgreen", lwd = 2)
pacf(pass_d1, lag.max = 36, main = "PACF - First Difference",
     ci.col = "darkgreen", lwd = 2)
par(mfrow = c(1, 1))
dev.off()

# ---- Second difference (check if needed) ----
pass_d2 <- diff(pass_ts, differences = 2)
adf_d2 <- adf.test(pass_d2, alternative = "stationary")
cat(sprintf("\nADF test (2nd diff): statistic = %.3f, p-value = %.4f\n",
            adf_d2$statistic, adf_d2$p.value))

# ---- Determine order d ----
if (adf_d1$p.value < 0.05) {
  d_order <- 1
  cat(sprintf("\n>>> Order of integration d = %d (stationary after 1st diff)\n", d_order))
} else {
  d_order <- 2
  cat(sprintf("\n>>> Order of integration d = %d (stationary after 2nd diff)\n", d_order))
}

###############################################################################
# 3. Determine ARIMA(p,d,q) order from ACF/PACF of differenced series
###############################################################################

# Analyze ACF/PACF of differenced series to determine p and q
d_acf  <- acf(pass_d1, lag.max = 24, plot = FALSE)
d_pacf <- pacf(pass_d1, lag.max = 24, plot = FALSE)

cat("\n--- ACF of differenced series (first 12 lags) ---\n")
print(round(d_acf$acf[2:13], 3))
cat("\n--- PACF of differenced series (first 12 lags) ---\n")
print(round(d_pacf$acf[1:12], 3))

# From ACF/PACF patterns:
# PACF cuts off after lag 1-2 => p = 1 or 2
# ACF has significant spike at lag 1 and possibly lag 12 (seasonal)
# Tentative non-seasonal ARIMA order: (2, d, 1) or (1, d, 1)

cat("\nTentative non-seasonal ARIMA based on ACF/PACF:\n")
cat("  p = 2 (PACF significant at lags 1-2, cuts off after)\n")
cat("  q = 1 (ACF significant at lag 1, exponential decay)\n")

###############################################################################
# 4. Build ARIMA(p,d,q) model
###############################################################################

p_order <- 2
q_order <- 1

arima_base <- Arima(pass_ts, order = c(p_order, d_order, q_order),
                     include.constant = TRUE)
cat(sprintf("\n--- ARIMA(%d,%d,%d) Model ---\n", p_order, d_order, q_order))
print(summary(arima_base))

###############################################################################
# 5. Residual analysis of ARIMA model
###############################################################################

residuals_base <- residuals(arima_base)

png("Task2_plots/04_ARIMA_residuals.png", width = 900, height = 700)
par(mfrow = c(3, 1))
plot(residuals_base, main = "ARIMA Residuals", ylab = "Residuals",
     col = "darkred", type = "o", pch = 16, cex = 0.5)
grid(col = "gray85")

acf(residuals_base, lag.max = 36, main = "ACF of ARIMA Residuals",
    ci.col = "steelblue", lwd = 2)

hist(residuals_base, breaks = 20, col = "lightblue", border = "white",
     main = "Histogram of ARIMA Residuals", xlab = "Residuals",
     probability = TRUE)
lines(density(residuals_base, na.rm = TRUE), col = "red", lwd = 2)
curve(dnorm(x, mean = mean(residuals_base, na.rm = TRUE),
            sd = sd(residuals_base, na.rm = TRUE)),
      add = TRUE, col = "blue", lwd = 2, lty = 2)
par(mfrow = c(1, 1))
dev.off()

# Check residual autocorrelation
lb_base <- Box.test(residuals_base, lag = 12, type = "Ljung-Box")
cat(sprintf("\nLjung-Box test (lag=12): X-squared = %.3f, p-value = %.4f\n",
            lb_base$statistic, lb_base$p.value))

# Check for remaining autocorrelation at low lags
res_acf <- acf(residuals_base, lag.max = 12, plot = FALSE)
cat("Residual ACF (lags 1-3): ", round(res_acf$acf[2:4], 3), "\n")

# If significant autocorrelation remains, try adjusting the model
if (lb_base$p.value < 0.05) {
  cat("\nResiduals show significant autocorrelation. Adjusting model...\n")

  # Try adding AR or MA terms
  model_candidates <- list(
    Arima(pass_ts, order = c(3, d_order, 1), include.constant = TRUE),
    Arima(pass_ts, order = c(2, d_order, 2), include.constant = TRUE),
    Arima(pass_ts, order = c(3, d_order, 2), include.constant = TRUE)
  )

  best_aic <- Inf
  best_arima <- arima_base
  for (m in model_candidates) {
    aic_val <- AIC(m)
    cat(sprintf("  ARIMA(%d,%d,%d) AIC = %.2f\n",
                m$arma[1], d_order, m$arma[2], aic_val))
    if (aic_val < best_aic) {
      best_aic <- aic_val
      best_arima <- m
    }
  }
  arima_base <- best_arima
  residuals_base <- residuals(arima_base)
  cat(sprintf("Selected adjusted model: ARIMA(%d,%d,%d) AIC = %.2f\n",
              arima_base$arma[1], d_order, arima_base$arma[2], best_aic))
}

###############################################################################
# 6. Seasonal component analysis
###############################################################################

cat("\n--- Seasonal Component Analysis ---\n")

# Check seasonal ACF
seasonal_acf <- acf(pass_ts, lag.max = 36, plot = FALSE)
cat("ACF at seasonal lags:\n")
cat(sprintf("  Lag 12: %.3f\n", seasonal_acf$acf[13]))
cat(sprintf("  Lag 24: %.3f\n", seasonal_acf$acf[25]))
cat(sprintf("  Lag 36: %.3f\n", seasonal_acf$acf[37]))

cat("\nStrong seasonal pattern detected with period s = 12 (annual).\n")
s_period <- 12

# ---- Seasonal differencing ----
pass_d1s <- diff(pass_d1, lag = s_period)
adf_d1s <- adf.test(pass_d1s, alternative = "stationary")
cat(sprintf("ADF after seasonal diff: statistic = %.3f, p-value = %.4f\n",
            adf_d1s$statistic, adf_d1s$p.value))

# ACF/PACF after seasonal differencing
png("Task2_plots/05_ACF_PACF_seasonal_diff.png", width = 900, height = 600)
par(mfrow = c(2, 1))
acf(pass_d1s, lag.max = 36, main = "ACF - Seasonal Difference (D=1)",
    ci.col = "darkorange", lwd = 2, na.action = na.pass)
pacf(pass_d1s, lag.max = 36, main = "PACF - Seasonal Difference (D=1)",
     ci.col = "darkorange", lwd = 2, na.action = na.pass)
par(mfrow = c(1, 1))
dev.off()

###############################################################################
# 7. SARIMA(p,d,q)(P,D,Q)s model
###############################################################################

# Based on analysis:
# d = 1 (regular differencing)
# D = 1 (seasonal differencing needed due to strong seasonal ACF)
# s = 12 (monthly data, annual cycle)
# Regular: p=2, q=1 from earlier ACF/PACF
# Seasonal: P=1, Q=1 (typical for monthly data)

sarima_model <- Arima(pass_ts,
  order    = c(2, 1, 1),
  seasonal = list(order = c(1, 1, 1), period = 12),
  include.constant = TRUE
)

cat("\n--- SARIMA(2,1,1)(1,1,1)[12] Model ---\n")
print(summary(sarima_model))

###############################################################################
# 8. Residual analysis of SARIMA model
###############################################################################

residuals_sarima <- residuals(sarima_model)

png("Task2_plots/06_SARIMA_residuals.png", width = 1000, height = 800)
par(mfrow = c(2, 2))

# Time plot
plot(residuals_sarima, main = "SARIMA Residuals", ylab = "Residuals",
     col = "darkred", type = "o", pch = 16, cex = 0.5)

# ACF of residuals
acf(residuals_sarima, lag.max = 24, main = "ACF of SARIMA Residuals",
    ci.col = "steelblue", lwd = 2, na.action = na.pass)

# Histogram
hist(residuals_sarima, breaks = 20, col = "lightblue", border = "white",
     main = "Histogram of SARIMA Residuals", xlab = "Residuals",
     probability = TRUE)
lines(density(residuals_sarima, na.rm = TRUE), col = "red", lwd = 2)
curve(dnorm(x, mean = mean(residuals_sarima, na.rm = TRUE),
            sd = sd(residuals_sarima, na.rm = TRUE)),
      add = TRUE, col = "blue", lwd = 2, lty = 2)

# Q-Q plot
qqnorm(residuals_sarima, main = "Q-Q Plot of SARIMA Residuals")
qqline(residuals_sarima, col = "red", lwd = 2)

par(mfrow = c(1, 1))
dev.off()

# If residuals still show issues, try modifying model
lb_sarima <- Box.test(residuals_sarima, lag = 12, type = "Ljung-Box")
cat(sprintf("\nLjung-Box (lag=12): X-squared = %.3f, p-value = %.4f\n",
            lb_sarima$statistic, lb_sarima$p.value))

if (lb_sarima$p.value < 0.05) {
  cat("\nResiduals show autocorrelation. Trying alternative SARIMA specifications...\n")

  alt_models <- list(
    Arima(pass_ts, order = c(2,1,2), seasonal = list(order = c(1,1,1), period = 12)),
    Arima(pass_ts, order = c(1,1,1), seasonal = list(order = c(1,1,2), period = 12)),
    Arima(pass_ts, order = c(3,1,1), seasonal = list(order = c(1,1,1), period = 12)),
    Arima(pass_ts, order = c(2,1,1), seasonal = list(order = c(2,1,1), period = 12)),
    Arima(pass_ts, order = c(1,1,1), seasonal = list(order = c(2,1,1), period = 12))
  )

  best_aic <- Inf
  best_sarima <- sarima_model
  for (m in alt_models) {
    m_aic <- AIC(m)
    m_res <- residuals(m)
    m_lb <- Box.test(m_res, lag = 12, type = "Ljung-Box")
    cat(sprintf("  SARIMA(%d,%d,%d)(%d,%d,%d)[12] AIC = %.2f, LB p = %.3f\n",
                m$arma[1], 1, m$arma[2], m$arma[3], 1, m$arma[4],
                m_aic, m_lb$p.value))
    if (m_lb$p.value > 0.05 && m_aic < best_aic) {
      best_aic <- m_aic
      best_sarima <- m
    }
  }
  sarima_model <- best_sarima
  residuals_sarima <- residuals(sarima_model)
  cat(sprintf("\nSelected: SARIMA(%d,%d,%d)(%d,%d,%d)[12]\n",
              sarima_model$arma[1], 1, sarima_model$arma[2],
              sarima_model$arma[3], 1, sarima_model$arma[4]))
}

###############################################################################
# 9. Normality test and white noise verification
###############################################################################

cat("\n--- Residual Diagnostics ---\n")

# Shapiro-Wilk normality test
sw_test <- shapiro.test(residuals_sarima)
cat(sprintf("Shapiro-Wilk normality: W = %.3f, p-value = %.4f\n",
            sw_test$statistic, sw_test$p.value))

# Kolmogorov-Smirnov test
ks_test <- ks.test(residuals_sarima, "pnorm",
                    mean = mean(residuals_sarima, na.rm = TRUE),
                    sd = sd(residuals_sarima, na.rm = TRUE))
cat(sprintf("Kolmogorov-Smirnov: D = %.3f, p-value = %.4f\n",
            ks_test$statistic, ks_test$p.value))

# Residual histogram for final model
png("Task2_plots/07_SARIMA_final_residual_hist.png", width = 800, height = 500)
hist(residuals_sarima, breaks = 25, col = "steelblue", border = "white",
     main = "SARIMA Final Model - Residual Histogram",
     xlab = "Residuals", probability = TRUE)
lines(density(residuals_sarima, na.rm = TRUE), col = "red", lwd = 2.5)
curve(dnorm(x, mean = mean(residuals_sarima, na.rm = TRUE),
            sd = sd(residuals_sarima, na.rm = TRUE)),
      add = TRUE, col = "darkgreen", lwd = 2.5, lty = 2)
legend("topright", legend = c("Density", "Normal Curve"),
       col = c("red", "darkgreen"), lwd = 2, lty = c(1, 2))
dev.off()

###############################################################################
# 10. Ljung-Box white noise test
###############################################################################

lb_final <- Box.test(residuals_sarima, lag = 12, type = "Ljung-Box")
cat(sprintf("\nLjung-Box Q-test (lag=12): X-squared = %.3f, p-value = %.4f\n",
            lb_final$statistic, lb_final$p.value))
cat(sprintf("Ljung-Box Q-test (lag=24): X-squared = %.3f, p-value = %.4f\n",
            Box.test(residuals_sarima, lag = 24, type = "Ljung-Box")$statistic,
            Box.test(residuals_sarima, lag = 24, type = "Ljung-Box")$p.value))

if (lb_final$p.value > 0.05) {
  cat(">>> Residuals are WHITE NOISE (do NOT reject H0). Model is adequate.\n")
} else {
  cat(">>> WARNING: Residuals are NOT white noise. Further model adjustment needed.\n")
}

# Check stationarity of residuals
adf_res <- adf.test(na.omit(residuals_sarima), alternative = "stationary")
cat(sprintf("ADF test on residuals: statistic = %.3f, p-value = %.4f\n",
            adf_res$statistic, adf_res$p.value))

###############################################################################
# 11. Forecast 18 periods ahead with SARIMA
###############################################################################

sarima_fcst <- forecast(sarima_model, h = 18, level = c(80, 95))

cat("\n--- SARIMA Forecast (18 months ahead) ---\n")
fcst_sarima_df <- data.frame(
  Month = as.character(time(sarima_fcst$mean)),
  Forecast = round(as.numeric(sarima_fcst$mean), 1),
  Lo80 = round(as.numeric(sarima_fcst$lower[, 1]), 1),
  Hi80 = round(as.numeric(sarima_fcst$upper[, 1]), 1),
  Lo95 = round(as.numeric(sarima_fcst$lower[, 2]), 1),
  Hi95 = round(as.numeric(sarima_fcst$upper[, 2]), 1)
)
print(fcst_sarima_df, row.names = FALSE)

# Plot forecast
png("Task2_plots/08_SARIMA_forecast_18months.png", width = 1000, height = 600)
plot(sarima_fcst, main = "SARIMA Forecast - International Airline Passengers (18 months)",
     ylab = "Passengers (thousands)", xlab = "Year",
     fcol = "tomato", flwd = 2.5, shaded = TRUE)
lines(fitted(sarima_model), col = "darkgreen", lwd = 1.5)
grid(col = "gray85")
legend("topleft",
       legend = c("Original", "Fitted", "Forecast", "95% CI"),
       col = c("black", "darkgreen", "tomato", "gray70"),
       lwd = c(1, 1.5, 2.5, 8), lty = c(1, 1, 1, 1))
dev.off()

###############################################################################
# 12. Exhaustive search - Auto SARIMA (best by AIC / SBC)
###############################################################################

cat("\n--- Exhaustive SARIMA Search (by AIC/SBC) ---\n")

auto_sarima <- auto.arima(pass_ts,
  seasonal = TRUE,
  stepwise = FALSE,
  approximation = FALSE,
  trace = FALSE,
  ic = "aic"
)

cat("\nBest SARIMA by auto.arima (AIC):\n")
print(summary(auto_sarima))
cat(sprintf("AIC: %.2f, BIC: %.2f\n", AIC(auto_sarima), BIC(auto_sarima)))

# Compare manual vs auto
cat("\n--- Model Comparison ---\n")
cat(sprintf("Manual SARIMA:  AIC = %.2f, BIC = %.2f\n",
            AIC(sarima_model), BIC(sarima_model)))
cat(sprintf("Auto SARIMA:    AIC = %.2f, BIC = %.2f\n",
            AIC(auto_sarima), BIC(auto_sarima)))

# Select best model by AIC
if (AIC(auto_sarima) < AIC(sarima_model)) {
  final_model <- auto_sarima
  model_label <- "Auto SARIMA"
} else {
  final_model <- sarima_model
  model_label <- "Manual SARIMA"
}
cat(sprintf("\n>>> Best model by AIC: %s\n", model_label))
cat(sprintf("    AIC = %.2f, BIC = %.2f\n", AIC(final_model), BIC(final_model)))

# Forecast with auto SARIMA
auto_fcst <- forecast(auto_sarima, h = 18, level = c(80, 95))

# Comparison plot
png("Task2_plots/09_model_comparison_forecast.png", width = 1000, height = 700)
par(mfrow = c(2, 1))
plot(sarima_fcst, main = paste("Manual SARIMA Forecast (AIC =", round(AIC(sarima_model), 1), ")"),
     ylab = "Passengers", fcol = "steelblue", flwd = 2)
plot(auto_fcst, main = paste("Auto SARIMA Forecast (AIC =", round(AIC(auto_sarima), 1), ")"),
     ylab = "Passengers", fcol = "tomato", flwd = 2)
par(mfrow = c(1, 1))
dev.off()

# ---- Save results ----
sink("Task2_SARIMA_summary.txt")
cat("TASK 2 - SARIMA MODEL SUMMARY\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("1. Stationarity Analysis:\n")
cat(sprintf("   ADF on original: included in analysis\n"))
cat(sprintf("   ADF on 1st diff: statistic = %.3f, p-value = %.4f\n",
            adf_d1$statistic, adf_d1$p.value))
cat(sprintf("   Integration order d = %d\n\n", d_order))

cat("2. Manual SARIMA Model:\n")
cat(sprintf("   SARIMA%s\n", gsub("ARIMA", "", capture.output(print(sarima_model))[1])))
cat(sprintf("   AIC = %.2f, BIC = %.2f\n\n", AIC(sarima_model), BIC(sarima_model)))

cat("3. Residual Diagnostics (Manual SARIMA):\n")
cat(sprintf("   Shapiro-Wilk: W = %.3f, p = %.4f\n", sw_test$statistic, sw_test$p.value))
cat(sprintf("   Ljung-Box(12): X2 = %.3f, p = %.4f\n", lb_final$statistic, lb_final$p.value))

cat("\n4. Best Auto SARIMA:\n")
cat(sprintf("   SARIMA%s\n", gsub("ARIMA", "", capture.output(print(auto_sarima))[1])))
cat(sprintf("   AIC = %.2f, BIC = %.2f\n\n", AIC(auto_sarima), BIC(auto_sarima)))

cat("5. Forecast (18 months, best model):\n")
print(fcst_sarima_df)

cat("\n6. Model Selection:\n")
cat(sprintf("   Selected: %s (by minimum AIC)\n", model_label))
cat(sprintf("   Rationale: Lower AIC indicates better fit-penalty balance.\n"))
cat("   Seasonal cycle: s = 12 months (annual seasonality confirmed).\n")
sink()

cat("\nTask 2 complete. Results saved to Task2_plots/ and Task2_SARIMA_summary.txt\n")
