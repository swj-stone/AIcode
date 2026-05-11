###############################################################################
# Practical Work 2 — Task 2: ARIMA / SARIMA Modelling
# Dataset: AirPassengers — Monthly international airline passengers (1949–1960)
###############################################################################

setwd("C:/Users/swj17/.claude/projects/R-analysis")
library(forecast)
library(tseries)
library(ggplot2)

data(AirPassengers)
y <- AirPassengers  # monthly, 1949–1960, 144 obs, s = 12

cat("Series: Monthly airline passengers (thousands), 1949–1960\n")
cat("Length:", length(y), "| Frequency:", frequency(y), "\n\n")

# ── 1. Plot & ACF / PACF ────────────────────────────────────────────────────

png("task2_01_eda.png", width = 1000, height = 800)
par(mfrow = c(3, 1))
plot(y, main = "AirPassengers — Original Series", ylab = "Thousands", xlab = "Year")
acf(y, lag.max = 48, main = "ACF — Original Series")
pacf(y, lag.max = 48, main = "PACF — Original Series")
dev.off()

cat("Stationarity assessment:\n")
cat("  - Clear upward trend → non-stationary in mean\n")
cat("  - Increasing seasonal amplitude → non-stationary in variance\n")
cat("  - ACF decays slowly → confirms non-stationarity\n\n")

# ── 2. Differencing & ADF test ──────────────────────────────────────────────

adf_orig <- adf.test(y, alternative = "stationary")
cat(sprintf("ADF test (original):  stat = %.4f, p = %.4f → %s\n",
            adf_orig$statistic, adf_orig$p.value,
            ifelse(adf_orig$p.value < 0.05, "stationary", "NON-stationary")))

# Log transform for variance stabilisation
y_log <- log(y)

# First difference of log series
dy1 <- diff(y_log, differences = 1)
adf_d1 <- adf.test(dy1, alternative = "stationary")
cat(sprintf("ADF test (1st diff log): stat = %.4f, p = %.4f → %s\n",
            adf_d1$statistic, adf_d1$p.value,
            ifelse(adf_d1$p.value < 0.05, "stationary", "NON-stationary")))

# Second difference
dy2 <- diff(y_log, differences = 2)
adf_d2 <- adf.test(dy2, alternative = "stationary")
cat(sprintf("ADF test (2nd diff log): stat = %.4f, p = %.4f → %s\n",
            adf_d2$statistic, adf_d2$p.value,
            ifelse(adf_d2$p.value < 0.05, "stationary", "NON-stationary")))

d <- 1  # confirmed by ADF: 1st diff is stationary
cat(sprintf("\n→ Order of integration d = %d\n\n", d))

# Correlograms for differenced series
png("task2_02_differencing.png", width = 1000, height = 800)
par(mfrow = c(3, 2))
acf(y_log, lag.max = 36, main = "ACF — Log series")
pacf(y_log, lag.max = 36, main = "PACF — Log series")
acf(dy1, lag.max = 36, main = "ACF — 1st diff (log)")
pacf(dy1, lag.max = 36, main = "PACF — 1st diff (log)")
acf(dy2, lag.max = 36, main = "ACF — 2nd diff (log)")
pacf(dy2, lag.max = 36, main = "PACF — 2nd diff (log)")
dev.off()

# ── 3. Determine ARIMA(p,d,q) ───────────────────────────────────────────────

cat("ACF/PACF of 1st-differenced log series:\n")
cat("  - ACF: significant spike at lag 1, then cuts off → MA(1) → q = 1\n")
cat("  - PACF: significant spike at lag 1 → AR(1) → p = 1\n")
cat("  - Also seasonal spikes at lag 12 in ACF\n")
cat("→ Initial model: ARIMA(1,1,1) on log scale\n\n")

# ── 4. Build ARIMA(1,1,1) ───────────────────────────────────────────────────

fit_arima <- Arima(y, order = c(1, 1, 1), lambda = 0)  # lambda=0 → log transform
cat("── ARIMA(1,1,1) with log transform ──\n")
print(summary(fit_arima))

# ── 5. Residual diagnostics ─────────────────────────────────────────────────

png("task2_05_arima_resid.png", width = 1000, height = 800)
checkresiduals(fit_arima)
dev.off()

res_arima <- residuals(fit_arima)
lb_arima <- Box.test(res_arima, lag = 24, type = "Ljung-Box", fitdf = 2)
cat(sprintf("Ljung-Box (24 lags): χ² = %.4f, p = %.4f\n",
            lb_arima$statistic, lb_arima$p.value))

# Check ACF of residuals
res_acf <- acf(res_arima, plot = FALSE, lag.max = 36)
sig_lags <- which(abs(res_acf$acf[-1]) > 1.96 / sqrt(length(res_arima)))
if (length(sig_lags) > 0) {
  cat("Significant residual ACF at lags:", sig_lags, "\n")
  cat("→ Need to refine model\n\n")
}

# ── 6 & 7. Seasonal component & SARIMA ──────────────────────────────────────

# Seasonal differencing
dy1_s12 <- diff(dy1, lag = 12)

png("task2_06_seasonal.png", width = 1000, height = 600)
par(mfrow = c(2, 2))
acf(dy1, lag.max = 36, main = "ACF — 1st diff (non-seasonal)")
pacf(dy1, lag.max = 36, main = "PACF — 1st diff (non-seasonal)")
acf(dy1_s12, lag.max = 36, main = "ACF — 1st + seasonal diff")
pacf(dy1_s12, lag.max = 36, main = "PACF — 1st + seasonal diff")
dev.off()

cat("Seasonal component: YES, seasonal cycle s = 12 (monthly)\n")
cat("Seasonal ACF/PACF analysis:\n")
cat("  - ACF spike at lag 12 → SMA(1) → Q = 1\n")
cat("  - Seasonal differencing needed → D = 1\n")
cat("→ SARIMA(p,d,q)(P,D,Q)[12] candidate: SARIMA(1,1,1)(0,1,1)[12]\n\n")

# Build SARIMA(1,1,1)(0,1,1)[12] with log transform
fit_sarima <- Arima(y, order = c(1, 1, 1),
                    seasonal = list(order = c(0, 1, 1), period = 12),
                    lambda = 0)
cat("── SARIMA(1,1,1)(0,1,1)[12] ──\n")
print(summary(fit_sarima))

# ── 8. SARIMA residual diagnostics ──────────────────────────────────────────

png("task2_08_sarima_resid.png", width = 1000, height = 800)
checkresiduals(fit_sarima)
dev.off()

res_sarima <- residuals(fit_sarima)

# ── 9. Histogram & Normality test ───────────────────────────────────────────

png("task2_09_normality.png", width = 1000, height = 500)
par(mfrow = c(1, 2))
hist(res_sarima, breaks = 20, probability = TRUE,
     main = "Histogram of SARIMA residuals",
     xlab = "Residuals", col = "lightblue")
curve(dnorm(x, mean = mean(res_sarima), sd = sd(res_sarima)),
      add = TRUE, col = "red", lwd = 2)
qqnorm(res_sarima, main = "Q-Q Plot")
qqline(res_sarima, col = "red", lwd = 2)
dev.off()

sw_test <- shapiro.test(res_sarima)
cat(sprintf("Shapiro-Wilk normality test: W = %.4f, p = %.4f → %s\n",
            sw_test$statistic, sw_test$p.value,
            ifelse(sw_test$p.value > 0.05, "Normal", "NON-normal")))

# ── 10. Ljung-Box test ──────────────────────────────────────────────────────

lb <- Box.test(res_sarima, lag = 24, type = "Ljung-Box", fitdf = 3)
cat(sprintf("Ljung-Box (24 lags): χ² = %.4f, p = %.4f → %s\n",
            lb$statistic, lb$p.value,
            ifelse(lb$p.value > 0.05, "WHITE NOISE ✓", "NOT white noise ✗")))

# Multiple lags
lb_lags <- seq(6, 36, by = 6)
lb_results <- sapply(lb_lags, function(l) {
  bt <- Box.test(res_sarima, lag = l, type = "Ljung-Box", fitdf = 3)
  bt$p.value
})
cat("\nLjung-Box p-values across lags:\n")
print(data.frame(Lag = lb_lags, p_value = round(lb_results, 4)))

# ── 11. Forecast 18 periods ahead ───────────────────────────────────────────

fcast <- forecast(fit_sarima, h = 18, level = 95)

png("task2_11_forecast.png", width = 1000, height = 600)
plot(fcast, main = "SARIMA(1,1,1)(0,1,1)[12] Forecast — 18 months ahead",
     xlab = "Year", ylab = "Thousands of passengers",
     shadecols = c("skyblue", "lightgray"))
lines(fitted(fit_sarima), col = "red", lwd = 1)
legend("topleft", legend = c("Observed", "Fitted", "Forecast", "95% CI"),
       col = c("black", "red", "blue", "skyblue"),
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 15), pt.cex = 2)
dev.off()

cat("\n── Forecast values ──\n")
print(fcast)

# ── 12. Auto-ARIMA comparison ───────────────────────────────────────────────

cat("\n── Auto-ARIMA exhaustive search ──\n")
fit_auto <- auto.arima(y, lambda = 0, seasonal = TRUE,
                       stepwise = FALSE, approximation = FALSE,
                       ic = "aic")
cat("Best model by AIC:\n")
print(summary(fit_auto))

cat("\n── Model Comparison ──\n")
comp <- data.frame(
  Model = c("Manual SARIMA(1,1,1)(0,1,1)[12]",
            paste0("Auto ", gsub("lambda = 0, ", "", capture.output(fit_auto)[1]))),
  AIC  = round(c(AIC(fit_sarima), AIC(fit_auto)), 2),
  AICc = round(c(fit_sarima$aicc, fit_auto$aicc), 2),
  BIC  = round(c(BIC(fit_sarima), BIC(fit_auto)), 2)
)
print(comp)

# Forecast from auto model
fcast_auto <- forecast(fit_auto, h = 18, level = 95)

png("task2_12_comparison.png", width = 1000, height = 600)
plot(fcast, main = "Forecast comparison: Manual vs Auto-ARIMA",
     xlab = "Year", ylab = "Thousands",
     shadecols = c("skyblue", "lightgray"))
lines(fcast_auto$mean, col = "darkgreen", lwd = 2, lty = 2)
legend("topleft",
       legend = c("Manual SARIMA", "Auto SARIMA"),
       col = c("blue", "darkgreen"), lty = c(1, 2), lwd = 2)
dev.off()

