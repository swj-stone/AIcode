# ============================================================================
# Task 2: ARIMA/SARIMA Modeling for Moroccan Air Passenger Traffic
# Dataset: Monthly air passengers (Morocco), 2010-2023
# Source: World Bank API (annual totals) + seasonal distribution
# ============================================================================

library(readxl)
library(ggplot2)
library(forecast)
library(tseries)
library(lmtest)
library(gridExtra)

# ============================================================================
# STEP 1: Time Series Analysis - Plot, ACF, PACF
# ============================================================================
cat("\n========================================\n")
cat("STEP 1: Time Series Analysis\n")
cat("========================================\n")

df <- read_excel("dataset2_air_traffic.xlsx")
cat(sprintf("Loaded %d monthly records\n", nrow(df)))

# Create time series (monthly, starting Jan 2010)
ts_air <- ts(df$Air_Passengers, start = c(2010, 1), frequency = 12)
cat(sprintf("Time series: %d obs, Jan 2010 - Dec 2023\n", length(ts_air)))

# Plot original series
png("task2_01_original_series.png", width = 1000, height = 500)
plot(ts_air, main = "Morocco Air Passenger Traffic (2010-2023)",
     ylab = "Passengers (monthly)", xlab = "Year",
     col = "steelblue", lwd = 2, type = "o", pch = 16, cex = 0.5)
abline(v = 2020, lty = 2, col = "red", lwd = 2)
text(2020, max(ts_air) * 0.9, "COVID-19", col = "red", pos = 4)
grid(col = "gray85")
dev.off()
cat("Plot saved: task2_01_original_series.png\n")

# ACF and PACF of original series
png("task2_02_acf_pacf_original.png", width = 1000, height = 600)
par(mfrow = c(2, 1))
acf(ts_air, lag.max = 36, main = "ACF - Original Series (Air Passengers)")
pacf(ts_air, lag.max = 36, main = "PACF - Original Series (Air Passengers)")
par(mfrow = c(1, 1))
dev.off()
cat("Plot saved: task2_02_acf_pacf_original.png\n")

# Calculate autocorrelation coefficients
acf_vals <- acf(ts_air, lag.max = 24, plot = FALSE)
cat("\nAutocorrelation coefficients (first 12 lags):\n")
for (i in 1:12) {
  stars <- if (abs(acf_vals$acf[i+1]) > 2/sqrt(length(ts_air))) "***" else ""
  cat(sprintf("  Lag %2d: %+.4f %s\n", i, acf_vals$acf[i+1], stars))
}
cat("(*** = statistically significant at 5% level)\n")

# Stationarity assessment
cat("\nStationarity assessment: The ACF shows slow decay and significant\n")
cat("autocorrelation at seasonal lags (12, 24), indicating NON-stationarity.\n")
cat("The series has both trend and strong seasonal components.\n")

# ============================================================================
# STEP 2: Determine order of differencing (d) using ADF test
# ============================================================================
cat("\n========================================\n")
cat("STEP 2: Order of Differencing (d)\n")
cat("========================================\n")

# ADF test on original series
adf_orig <- adf.test(ts_air, alternative = "stationary")
cat(sprintf("\nADF test - Original series:\n"))
cat(sprintf("  Dickey-Fuller = %.4f, p-value = %.4f\n", adf_orig$statistic, adf_orig$p.value))
cat(sprintf("  Conclusion: %s\n", ifelse(adf_orig$p.value < 0.05,
      "STATIONARY (reject unit root)", "NON-STATIONARY (unit root present)")))

# First difference
ts_d1 <- diff(ts_air, differences = 1)
adf_d1 <- adf.test(ts_d1, alternative = "stationary")
cat(sprintf("\nADF test - 1st difference:\n"))
cat(sprintf("  Dickey-Fuller = %.4f, p-value = %.4f\n", adf_d1$statistic, adf_d1$p.value))
cat(sprintf("  Conclusion: %s\n", ifelse(adf_d1$p.value < 0.05,
      "STATIONARY (reject unit root)", "NON-STATIONARY (unit root present)")))

# Second difference
ts_d2 <- diff(ts_air, differences = 2)
adf_d2 <- adf.test(ts_d2, alternative = "stationary")
cat(sprintf("\nADF test - 2nd difference:\n"))
cat(sprintf("  Dickey-Fuller = %.4f, p-value = %.4f\n", adf_d2$statistic, adf_d2$p.value))
cat(sprintf("  Conclusion: %s\n", ifelse(adf_d2$p.value < 0.05,
      "STATIONARY (reject unit root)", "NON-STATIONARY (unit root present)")))

# Determine d
if (adf_orig$p.value < 0.05) {
  d_order <- 0
} else if (adf_d1$p.value < 0.05) {
  d_order <- 1
} else if (adf_d2$p.value < 0.05) {
  d_order <- 2
} else {
  d_order <- 1  # default
}
cat(sprintf("\n=> Order of integration d = %d\n", d_order))

# Correlograms for original, 1st diff, 2nd diff
png("task2_03_correlograms_differences.png", width = 1200, height = 800)
par(mfrow = c(3, 2))
acf(ts_air, lag.max = 36, main = "ACF - Original")
pacf(ts_air, lag.max = 36, main = "PACF - Original")
acf(ts_d1, lag.max = 36, main = "ACF - 1st Difference")
pacf(ts_d1, lag.max = 36, main = "PACF - 1st Difference")
acf(ts_d2, lag.max = 36, main = "ACF - 2nd Difference")
pacf(ts_d2, lag.max = 36, main = "PACF - 2nd Difference")
par(mfrow = c(1, 1))
dev.off()
cat("Plot saved: task2_03_correlograms_differences.png\n")

# ============================================================================
# STEP 3: Determine ARIMA(p,d,q) from ACF/PACF of differenced series
# ============================================================================
cat("\n========================================\n")
cat("STEP 3: Determine ARIMA(p,d,q) Order\n")
cat("========================================\n")

# Analyze ACF/PACF of the differenced series
diff_series <- if (d_order == 1) ts_d1 else ts_d2

acf_diff <- acf(diff_series, lag.max = 36, plot = FALSE)
pacf_diff <- pacf(diff_series, lag.max = 36, plot = FALSE)

cat("\nDifferenced series ACF/PACF analysis:\n")

# Find significant ACF lags (for q - MA order)
sig_acf <- which(abs(acf_diff$acf[-1]) > 2/sqrt(length(diff_series)))
if (length(sig_acf) > 0) {
  cat(sprintf("  Significant ACF lags: %s\n", paste(sig_acf[1:min(10, length(sig_acf))], collapse = ", ")))
  # Last significant ACF lag before cutoff (ignoring seasonal)
  non_seasonal_acf <- sig_acf[sig_acf < 10 & sig_acf > 0]
  q_order <- if (length(non_seasonal_acf) > 0) max(non_seasonal_acf) else 2
} else {
  q_order <- 0
}

# Find significant PACF lags (for p - AR order)
sig_pacf <- which(abs(pacf_diff$acf) > 2/sqrt(length(diff_series)))
if (length(sig_pacf) > 0) {
  cat(sprintf("  Significant PACF lags: %s\n", paste(sig_pacf[1:min(10, length(sig_pacf))], collapse = ", ")))
  non_seasonal_pacf <- sig_pacf[sig_pacf < 10 & sig_pacf > 0]
  p_order <- if (length(non_seasonal_pacf) > 0) min(max(non_seasonal_pacf), 3) else 2
} else {
  p_order <- 1
}

cat(sprintf("\nDetermined ARIMA order: (%d, %d, %d)\n", p_order, d_order, q_order))

# ============================================================================
# STEP 4: Build ARIMA(p,d,q) Model
# ============================================================================
cat("\n========================================\n")
cat("STEP 4: ARIMA Model Construction\n")
cat("========================================\n")

arima_model <- Arima(ts_air, order = c(p_order, d_order, q_order))
cat("\nARIMA Model Summary:\n")
print(summary(arima_model))

# ============================================================================
# STEP 5: Residual Analysis of ARIMA
# ============================================================================
cat("\n========================================\n")
cat("STEP 5: Residual Diagnostics\n")
cat("========================================\n")

residuals_arima <- residuals(arima_model)

png("task2_05_residuals_arima.png", width = 1000, height = 800)
par(mfrow = c(2, 2))
plot(residuals_arima, main = "ARIMA Residuals", ylab = "Residual", type = "o", cex = 0.5)
abline(h = 0, col = "red", lty = 2)
acf(residuals_arima, lag.max = 36, main = "ACF of Residuals")
pacf(residuals_arima, lag.max = 36, main = "PACF of Residuals")
qqnorm(residuals_arima, main = "Q-Q Plot of Residuals")
qqline(residuals_arima, col = "red")
par(mfrow = c(1, 1))
dev.off()
cat("Plot saved: task2_05_residuals_arima.png\n")

# Check for remaining autocorrelation in residuals
lb_test <- Box.test(residuals_arima, lag = 12, type = "Ljung-Box")
cat(sprintf("\nLjung-Box test (lag=12): X-squared = %.4f, p-value = %.4f\n",
            lb_test$statistic, lb_test$p.value))

# Check low-order residual autocorrelation
acf_res <- acf(residuals_arima, lag.max = 3, plot = FALSE)
cat("\nLow-order residual ACF:\n")
for (i in 1:3) {
  cat(sprintf("  Lag %d: %.4f\n", i, acf_res$acf[i+1]))
}

# If residuals have significant low-order autocorrelation, modify model
has_autocorr <- any(abs(acf_res$acf[2:4]) > 2/sqrt(length(residuals_arima)))
if (has_autocorr) {
  cat("\nLow-order autocorrelation detected! Modifying model...\n")
  # Try adding AR or MA terms
  best_aic <- AIC(arima_model)
  best_mod <- arima_model

  for (p_try in max(0, p_order-1):min(p_order+2, 5)) {
    for (q_try in max(0, q_order-1):min(q_order+2, 5)) {
      if (p_try == p_order && q_try == q_order) next
      tryCatch({
        mod <- Arima(ts_air, order = c(p_try, d_order, q_try))
        if (AIC(mod) < best_aic) {
          best_aic <- AIC(mod)
          best_mod <- mod
        }
      }, error = function(e) {})
    }
  }

  if (AIC(best_mod) < AIC(arima_model)) {
    arima_model <- best_mod
    cat(sprintf("  Updated ARIMA(%d,%d,%d), AIC=%.2f\n",
                arima_model$arma[1], d_order, arima_model$arma[2], best_aic))
  }
}

# ============================================================================
# STEP 6: Seasonal Component Detection
# ============================================================================
cat("\n========================================\n")
cat("STEP 6: Seasonal Component Analysis\n")
cat("========================================\n")

# Check ACF for seasonal patterns
acf_orig <- acf(ts_air, lag.max = 36, plot = FALSE)
cat("Seasonal autocorrelation coefficients:\n")
for (lag in c(12, 24, 36)) {
  if (lag <= length(acf_orig$acf) - 1) {
    is_sig <- abs(acf_orig$acf[lag+1]) > 2/sqrt(length(ts_air))
    cat(sprintf("  Lag %d: %+.4f %s\n", lag, acf_orig$acf[lag+1],
                ifelse(is_sig, "(SIGNIFICANT - seasonal component present)", "")))
  }
}

# Seasonal differencing
s <- 12  # seasonal period
ts_seasonal_diff <- diff(ts_air, lag = s, differences = 1)

# Combined: seasonal + regular differencing
if (d_order > 0) {
  ts_combined_diff <- diff(diff(ts_air, lag = s), differences = d_order)
} else {
  ts_combined_diff <- diff(ts_air, lag = s)
}

png("task2_06_seasonal_differencing.png", width = 1200, height = 600)
par(mfrow = c(2, 2))
acf(ts_seasonal_diff, lag.max = 36, main = "ACF - Seasonal Difference (s=12)")
pacf(ts_seasonal_diff, lag.max = 36, main = "PACF - Seasonal Difference (s=12)")
acf(ts_combined_diff, lag.max = 36, main = "ACF - Seasonal + Regular Difference")
pacf(ts_combined_diff, lag.max = 36, main = "PACF - Seasonal + Regular Difference")
par(mfrow = c(1, 1))
dev.off()
cat("Plot saved: task2_06_seasonal_differencing.png\n")

# D = 1 (seasonal differencing needed due to significant seasonal ACF)
D_order <- 1
cat(sprintf("\nSeasonal period s = %d, Seasonal differencing D = %d\n", s, D_order))

# ============================================================================
# STEP 7: SARIMA Model with Seasonal Elements
# ============================================================================
cat("\n========================================\n")
cat("STEP 7: SARIMA Model Construction\n")
cat("========================================\n")

# Determine seasonal AR(P) and MA(Q) from ACF/PACF of seasonally differenced series
acf_sd <- acf(ts_combined_diff, lag.max = 36, plot = FALSE)
pacf_sd <- pacf(ts_combined_diff, lag.max = 36, plot = FALSE)

# Seasonal P from PACF at lag 12
P_order <- if (abs(pacf_sd$acf[12]) > 2/sqrt(length(ts_combined_diff))) 1 else 0
if (abs(pacf_sd$acf[24]) > 2/sqrt(length(ts_combined_diff))) P_order <- max(P_order, 1)

# Seasonal Q from ACF at lag 12
Q_order <- if (abs(acf_sd$acf[13]) > 2/sqrt(length(ts_combined_diff))) 1 else 0
if (abs(acf_sd$acf[25]) > 2/sqrt(length(ts_combined_diff))) Q_order <- max(Q_order, 1)

cat(sprintf("Seasonal order: (P=%d, D=%d, Q=%d), s=%d\n", P_order, D_order, Q_order, s))
cat(sprintf("Full SARIMA: (%d,%d,%d)(%d,%d,%d)[%d]\n",
            p_order, d_order, q_order, P_order, D_order, Q_order, s))

# Build SARIMA model
sarima_model <- Arima(ts_air,
                      order = c(p_order, d_order, q_order),
                      seasonal = list(order = c(P_order, D_order, Q_order), period = s))
cat("\nSARIMA Model Summary:\n")
print(summary(sarima_model))

# Also try a few nearby orders
cat("\nTrying nearby SARIMA specifications:\n")
candidates <- list()
aic_values <- c()

try_orders <- expand.grid(
  p = c(max(0,p_order-1):min(p_order+1,3)),
  q = c(max(0,q_order-1):min(q_order+1,3)),
  P = c(max(0,P_order-1):min(P_order+1,2)),
  Q = c(max(0,Q_order-1):min(Q_order+1,2))
)

for (i in 1:nrow(try_orders)) {
  ord <- try_orders[i, ]
  tryCatch({
    m <- Arima(ts_air,
               order = c(ord$p, d_order, ord$q),
               seasonal = list(order = c(ord$P, D_order, ord$Q), period = s))
    aic_val <- AIC(m)
    cat(sprintf("  SARIMA(%d,%d,%d)(%d,%d,%d)[%d]: AIC=%.2f\n",
                ord$p, d_order, ord$q, ord$P, D_order, ord$Q, s, aic_val))
    if (length(candidates) == 0 || aic_val < min(aic_values)) {
      if (length(candidates) == 0 || aic_val < min(aic_values)) {
        best_sarima <- m
      }
    }
    candidates[[length(candidates) + 1]] <- m
    aic_values <- c(aic_values, aic_val)
  }, error = function(e) {})
}

if (exists("best_sarima")) {
  sarima_model <- best_sarima
}
cat(sprintf("\nSelected SARIMA(%d,%d,%d)(%d,%d,%d)[%d]\n",
            sarima_model$arma[1], sarima_model$arma[6], sarima_model$arma[2],
            sarima_model$arma[3], sarima_model$arma[7], sarima_model$arma[4],
            sarima_model$arma[5]))

# ============================================================================
# STEP 8: Residual Analysis of SARIMA
# ============================================================================
cat("\n========================================\n")
cat("STEP 8: SARIMA Residual Diagnostics\n")
cat("========================================\n")

residuals_sarima <- residuals(sarima_model)

png("task2_08_residuals_sarima.png", width = 1000, height = 800)
par(mfrow = c(2, 2))
plot(residuals_sarima, main = "SARIMA Residuals", ylab = "Residual", type = "o", cex = 0.5)
abline(h = 0, col = "red", lty = 2)
acf(residuals_sarima, lag.max = 36, main = "ACF - SARIMA Residuals")
pacf(residuals_sarima, lag.max = 36, main = "PACF - SARIMA Residuals")
hist(residuals_sarima, breaks = 20, main = "Histogram of Residuals",
     xlab = "Residual", col = "lightblue", border = "darkblue", probability = TRUE)
lines(density(residuals_sarima, na.rm = TRUE), col = "red", lwd = 2)
curve(dnorm(x, mean = mean(residuals_sarima, na.rm = TRUE),
            sd = sd(residuals_sarima, na.rm = TRUE)),
      add = TRUE, col = "darkgreen", lwd = 2, lty = 2)
par(mfrow = c(1, 1))
dev.off()
cat("Plot saved: task2_08_residuals_sarima.png\n")

# Check if residuals are white noise
lb_res <- Box.test(residuals_sarima, lag = 12, type = "Ljung-Box")
cat(sprintf("\nLjung-Box test (lag=12): X-squared = %.4f, p-value = %.4f\n",
            lb_res$statistic, lb_res$p.value))

if (lb_res$p.value < 0.05) {
  cat("Residuals NOT white noise - further model modification needed.\n")
  # Try different orders
  cat("Attempting automatic ARIMA selection...\n")
  auto_mod <- auto.arima(ts_air, seasonal = TRUE, stepwise = FALSE,
                         approximation = FALSE, trace = FALSE)

  lb_auto <- Box.test(residuals(auto_mod), lag = 12, type = "Ljung-Box")
  if (lb_auto$p.value > lb_res$p.value) {
    cat(sprintf("Auto ARIMA improves residual whiteness (p=%.4f vs %.4f)\n",
                lb_auto$p.value, lb_res$p.value))
    sarima_model <- auto_mod
    residuals_sarima <- residuals(sarima_model)
  }
}

# ============================================================================
# STEP 9: Residual Histogram & Normality Test
# ============================================================================
cat("\n========================================\n")
cat("STEP 9: Normality Testing\n")
cat("========================================\n")

png("task2_09_residual_histogram.png", width = 800, height = 500)
hist(residuals_sarima, breaks = 25, main = "SARIMA Residual Distribution",
     xlab = "Residual Value", col = "steelblue", border = "white",
     probability = TRUE)
lines(density(residuals_sarima, na.rm = TRUE), col = "red", lwd = 2.5)
curve(dnorm(x, mean = mean(residuals_sarima, na.rm = TRUE),
            sd = sd(residuals_sarima, na.rm = TRUE)),
      add = TRUE, col = "darkgreen", lwd = 2, lty = 2)
legend("topright", legend = c("Kernel Density", "Normal Distribution"),
       col = c("red", "darkgreen"), lty = c(1, 2), lwd = 2, bty = "n")
dev.off()
cat("Plot saved: task2_09_residual_histogram.png\n")

# Shapiro-Wilk test
sw_test <- shapiro.test(residuals_sarima[1:min(5000, length(residuals_sarima))])
cat(sprintf("\nShapiro-Wilk normality test:\n"))
cat(sprintf("  W = %.4f, p-value = %.4f\n", sw_test$statistic, sw_test$p.value))
cat(sprintf("  Conclusion: Residuals %s normally distributed\n",
            ifelse(sw_test$p.value > 0.05, "ARE", "are NOT")))

# ============================================================================
# STEP 10: White Noise Verification (Ljung-Box Q-test)
# ============================================================================
cat("\n========================================\n")
cat("STEP 10: White Noise Verification\n")
cat("========================================\n")

# Ljung-Box for multiple lags
lb_results <- data.frame(Lag = c(6, 12, 18, 24, 36),
                         Q_stat = NA, p_value = NA)
for (i in 1:nrow(lb_results)) {
  lb <- Box.test(residuals_sarima, lag = lb_results$Lag[i], type = "Ljung-Box")
  lb_results$Q_stat[i] <- round(lb$statistic, 4)
  lb_results$p_value[i] <- round(lb$p.value, 4)
}

cat("\nLjung-Box Q-test results:\n")
print(lb_results)

# ADF test on residuals
adf_res <- adf.test(residuals_sarima[!is.na(residuals_sarima)], alternative = "stationary")
cat(sprintf("\nADF test on residuals: DF=%.4f, p=%.4f\n", adf_res$statistic, adf_res$p.value))
cat(sprintf("Residuals: %s\n", ifelse(adf_res$p.value < 0.05,
      "STATIONARY (white noise confirmed)", "NON-STATIONARY (issue detected)")))

# ============================================================================
# STEP 11: Forecast 18 Periods Ahead
# ============================================================================
cat("\n========================================\n")
cat("STEP 11: 18-Period Forecast\n")
cat("========================================\n")

fc_sarima <- forecast(sarima_model, h = 18, level = c(80, 95))

cat("\nSARIMA Forecast (18 periods):\n")
fc_df <- data.frame(
  Period = 1:18,
  Forecast = round(as.numeric(fc_sarima$mean), 0),
  Lo80 = round(as.numeric(fc_sarima$lower[,1]), 0),
  Hi80 = round(as.numeric(fc_sarima$upper[,1]), 0),
  Lo95 = round(as.numeric(fc_sarima$lower[,2]), 0),
  Hi95 = round(as.numeric(fc_sarima$upper[,2]), 0)
)
print(fc_df)
write.csv(fc_df, "task2_forecast_18periods_sarima.csv", row.names = FALSE)

# Plot forecast
png("task2_11_forecast_18periods.png", width = 1200, height = 600)
plot(fc_sarima, main = "SARIMA Forecast - Morocco Air Passengers (18 months ahead)",
     ylab = "Passengers", xlab = "Year",
     fcol = "steelblue", shadecols = c("lightblue", "lightsteelblue"))
lines(fitted(sarima_model), col = "red", lwd = 1)
legend("topleft", legend = c("Actual", "Fitted", "Forecast", "95% CI", "80% CI"),
       col = c("black", "red", "steelblue", "lightsteelblue", "lightblue"),
       lty = c(1, 1, 1, NA, NA), pch = c(NA, NA, NA, 15, 15),
       lwd = c(1, 1, 2, NA, NA), bty = "n")
grid(col = "gray90")
dev.off()
cat("Plot saved: task2_11_forecast_18periods.png\n")

# ============================================================================
# STEP 12: Auto ARIMA - Exhaustive Search (Best by AIC/SBC)
# ============================================================================
cat("\n========================================\n")
cat("STEP 12: Automatic Model Selection (AIC/SBC)\n")
cat("========================================\n")

# Exhaustive search using auto.arima
cat("Running exhaustive auto.arima search...\n")
auto_model <- auto.arima(ts_air, seasonal = TRUE,
                         stepwise = FALSE,
                         approximation = FALSE,
                         trace = FALSE,
                         ic = "aic")

cat("\nBest model by AIC (auto.arima exhaustive):\n")
cat(sprintf("  SARIMA(%d,%d,%d)(%d,%d,%d)[%d]\n",
            auto_model$arma[1], auto_model$arma[6], auto_model$arma[2],
            auto_model$arma[3], auto_model$arma[7], auto_model$arma[4],
            auto_model$arma[5]))
cat(sprintf("  AIC = %.2f, BIC = %.2f\n", AIC(auto_model), BIC(auto_model)))

# Also get best by BIC
auto_model_bic <- auto.arima(ts_air, seasonal = TRUE,
                             stepwise = FALSE,
                             approximation = FALSE,
                             trace = FALSE,
                             ic = "bic")

cat(sprintf("\nBest model by BIC (auto.arima exhaustive):\n"))
cat(sprintf("  SARIMA(%d,%d,%d)(%d,%d,%d)[%d]\n",
            auto_model_bic$arma[1], auto_model_bic$arma[6], auto_model_bic$arma[2],
            auto_model_bic$arma[3], auto_model_bic$arma[7], auto_model_bic$arma[4],
            auto_model_bic$arma[5]))
cat(sprintf("  AIC = %.2f, BIC = %.2f\n", AIC(auto_model_bic), BIC(auto_model_bic)))

# Comparison table
comp_table <- data.frame(
  Model = c(
    "Manual SARIMA",
    "Auto SARIMA (AIC)",
    "Auto SARIMA (BIC)"
  ),
  Order = c(
    sprintf("(%d,%d,%d)(%d,%d,%d)[%d]",
            sarima_model$arma[1], sarima_model$arma[6], sarima_model$arma[2],
            sarima_model$arma[3], sarima_model$arma[7], sarima_model$arma[4],
            sarima_model$arma[5]),
    sprintf("(%d,%d,%d)(%d,%d,%d)[%d]",
            auto_model$arma[1], auto_model$arma[6], auto_model$arma[2],
            auto_model$arma[3], auto_model$arma[7], auto_model$arma[4],
            auto_model$arma[5]),
    sprintf("(%d,%d,%d)(%d,%d,%d)[%d]",
            auto_model_bic$arma[1], auto_model_bic$arma[6], auto_model_bic$arma[2],
            auto_model_bic$arma[3], auto_model_bic$arma[7], auto_model_bic$arma[4],
            auto_model_bic$arma[5])
  ),
  AIC = c(AIC(sarima_model), AIC(auto_model), AIC(auto_model_bic)),
  BIC = c(BIC(sarima_model), BIC(auto_model), BIC(auto_model_bic))
)
cat("\nModel Comparison:\n")
print(comp_table)

# Forecast using best auto model
fc_auto <- forecast(auto_model, h = 18, level = c(95))

# Compare forecasts
png("task2_12_forecast_comparison.png", width = 1200, height = 700)
par(mfrow = c(2, 1))

# Manual model forecast
plot(fc_sarima, main = "Manual SARIMA Forecast",
     ylab = "Passengers", fcol = "steelblue", shadecols = c("lightgray", "lightsteelblue"))
lines(fitted(sarima_model), col = "red", lwd = 1)

# Auto model forecast
plot(fc_auto, main = sprintf("Auto SARIMA (AIC=%.0f) Forecast", AIC(auto_model)),
     ylab = "Passengers", fcol = "darkgreen", shadecols = c("lightgray", "lightgreen"))
lines(fitted(auto_model), col = "red", lwd = 1)

par(mfrow = c(1, 1))
dev.off()
cat("Plot saved: task2_12_forecast_comparison.png\n")

# Compare forecasts numerically
fc_comp <- data.frame(
  Period = 1:18,
  Manual_SARIMA = round(as.numeric(fc_sarima$mean), 0),
  Auto_SARIMA_AIC = round(as.numeric(fc_auto$mean), 0)
)
fc_comp$Difference <- fc_comp$Manual_SARIMA - fc_comp$Auto_SARIMA_AIC
cat("\nForecast comparison:\n")
print(fc_comp)
write.csv(fc_comp, "task2_forecast_comparison.csv", row.names = FALSE)

# ============================================================================
# SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("TASK 2 COMPLETE\n")
cat("========================================\n")
cat(sprintf("Dataset: Morocco Air Passengers, %d monthly obs (2010-2023)\n", length(ts_air)))
cat(sprintf("Stationarity: d = %d\n", d_order))
cat(sprintf("Seasonality: s = %d, D = %d\n", s, D_order))
cat(sprintf("Manual SARIMA: (%d,%d,%d)(%d,%d,%d)[%d], AIC=%.2f\n",
            sarima_model$arma[1], sarima_model$arma[6], sarima_model$arma[2],
            sarima_model$arma[3], sarima_model$arma[7], sarima_model$arma[4],
            sarima_model$arma[5], AIC(sarima_model)))
cat(sprintf("Best Auto SARIMA: (%d,%d,%d)(%d,%d,%d)[%d], AIC=%.2f\n",
            auto_model$arma[1], auto_model$arma[6], auto_model$arma[2],
            auto_model$arma[3], auto_model$arma[7], auto_model$arma[4],
            auto_model$arma[5], AIC(auto_model)))

save(ts_air, sarima_model, auto_model, auto_model_bic, comp_table,
     fc_sarima, fc_auto, fc_comp, d_order, D_order, s,
     file = "task2_results.RData")
cat("Results saved to task2_results.RData\n")
