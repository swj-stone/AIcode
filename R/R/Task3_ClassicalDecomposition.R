###############################################################################
# Practical Work №2 - Task 3: Classical Decomposition
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

# Create time series
pass_ts <- ts(raw$Passengers, start = c(1949, 1), frequency = 12)
cat(sprintf("Loaded time series: %d months (1949-1960), frequency = 12\n", length(pass_ts)))

###############################################################################
# 1-2. Time series analysis, ACF and PACF
###############################################################################

dir.create("Task3_plots", showWarnings = FALSE)

# ---- Plot original series ----
png("Task3_plots/01_original_series.png", width = 900, height = 500)
plot(pass_ts, main = "International Airline Passengers (1949-1960)",
     ylab = "Passengers (thousands)", xlab = "Year",
     col = "steelblue", lwd = 2, type = "o", pch = 16, cex = 0.7)
grid(col = "gray85")
dev.off()

# ---- ACF and PACF ----
acf_vals <- acf(pass_ts, lag.max = 36, plot = FALSE)
pacf_vals <- pacf(pass_ts, lag.max = 36, plot = FALSE)

png("Task3_plots/02_ACF_PACF.png", width = 900, height = 600)
par(mfrow = c(2, 1))
acf(pass_ts, lag.max = 36, main = "ACF - Airline Passengers",
    ci.col = "steelblue", lwd = 2)
pacf(pass_ts, lag.max = 36, main = "PACF - Airline Passengers",
     ci.col = "steelblue", lwd = 2)
par(mfrow = c(1, 1))
dev.off()

cat("\n--- ACF Analysis ---\n")
cat(sprintf("Lag 1:  %.3f\n", acf_vals$acf[2]))
cat(sprintf("Lag 6:  %.3f\n", acf_vals$acf[7]))
cat(sprintf("Lag 12: %.3f (seasonal peak)\n", acf_vals$acf[13]))
cat(sprintf("Lag 24: %.3f (seasonal peak)\n", acf_vals$acf[25]))

cat("\n--- PACF Analysis ---\n")
cat(sprintf("Lag 1:  %.3f\n", pacf_vals$acf[1]))
cat(sprintf("Lag 12: %.3f (seasonal)\n", pacf_vals$acf[12]))
cat(sprintf("Lag 24: %.3f (seasonal)\n", pacf_vals$acf[24]))

###############################################################################
# 3. Seasonal component analysis
###############################################################################

cat("\n--- Seasonal Component Analysis ---\n")

# Check seasonal lags for significance
# ACF at lag 12 >> 2/sqrt(n) threshold => significant seasonality
threshold <- 2 / sqrt(length(pass_ts))
cat(sprintf("Significance threshold (2/sqrt(n)): %.3f\n", threshold))

seasonal_lags <- c(12, 24, 36)
for (lag in seasonal_lags) {
  acf_val <- acf_vals$acf[lag + 1]
  is_sig <- abs(acf_val) > threshold
  cat(sprintf("  ACF at lag %2d: %.3f %s\n", lag, acf_val,
              ifelse(is_sig, "(SIGNIFICANT - seasonality confirmed)", "")))
}

cat("\nSeasonality type analysis:\n")
cat("  - Period: s = 12 months (annual cycle)\n")
cat("  - Behavior: Multiplicative (seasonal amplitude grows with trend level)\n")
cat("  - Evidence: The spread between seasonal peaks and troughs widens over time.\n")
cat("             This is characteristic of multiplicative seasonality.\n")

# ---- Visual evidence of multiplicative seasonality ----
# Plot seasonal subseries
png("Task3_plots/03_seasonal_subseries.png", width = 900, height = 500)
monthplot(pass_ts, main = "Seasonal Subseries Plot - Airline Passengers",
          ylab = "Passengers (thousands)", xlab = "Month",
          col = "steelblue", lwd = 2, type = "o", pch = 16)
grid(col = "gray85")
dev.off()

# ---- Boxplot by month ----
png("Task3_plots/04_seasonal_boxplot.png", width = 800, height = 500)
boxplot(split(pass_ts, cycle(pass_ts)),
        main = "Monthly Distribution of Airline Passengers",
        xlab = "Month", ylab = "Passengers (thousands)",
        col = c("steelblue", "tomato", "darkgreen", "goldenrod",
                "purple", "orange", "cyan", "pink",
                "brown", "darkblue", "darkred", "darkgreen"),
        names = month.abb)
grid(col = "gray85")
dev.off()

###############################################################################
# 4. Classical decomposition
###############################################################################

# Try both additive and multiplicative decomposition
# Visual inspection suggests multiplicative (increasing seasonal amplitude)
decomp_mult <- decompose(pass_ts, type = "multiplicative")
decomp_add  <- decompose(pass_ts, type = "additive")

# ---- Plot multiplicative decomposition ----
png("Task3_plots/05_decomposition_multiplicative.png", width = 1000, height = 800)
plot(decomp_mult, col = "steelblue")
dev.off()

# ---- Plot additive decomposition (for comparison) ----
png("Task3_plots/06_decomposition_additive.png", width = 1000, height = 800)
plot(decomp_add, col = "darkgreen")
dev.off()

# ---- Extracted components ----
trend_mult  <- decomp_mult$trend
seasonal_mult <- decomp_mult$seasonal
random_mult  <- decomp_mult$random

cat("\n--- Classical Decomposition (Multiplicative) ---\n")
cat("Model: Y(t) = Trend(t) * Seasonal(t) * Random(t)\n\n")

cat("Trend component features:\n")
cat(sprintf("  Range: %.1f to %.1f\n", min(trend_mult, na.rm = TRUE),
            max(trend_mult, na.rm = TRUE)))
cat("  Pattern: Monotonically increasing, approximately linear growth\n\n")

cat("Seasonal component features:\n")
cat(sprintf("  Range: %.3f to %.3f\n", min(seasonal_mult, na.rm = TRUE),
            max(seasonal_mult, na.rm = TRUE)))
cat("  Pattern: Regular annual cycle, peak in Jul-Aug, trough in Nov\n")

# ---- Plot trend and seasonal separately ----
png("Task3_plots/07_trend_and_seasonal.png", width = 1000, height = 600)
par(mfrow = c(2, 1))
plot(trend_mult, main = "Trend Component (Multiplicative Decomposition)",
     ylab = "Trend", xlab = "Year", col = "steelblue", lwd = 2)
grid(col = "gray85")

plot(seasonal_mult, main = "Seasonal Component (Multiplicative Decomposition)",
     ylab = "Seasonal Factor", xlab = "Year",
     col = "tomato", lwd = 2)
grid(col = "gray85")
par(mfrow = c(1, 1))
dev.off()

###############################################################################
# 5. Error analysis of decomposition model
###############################################################################

# Multiplicative model: random = Y / (trend * seasonal) = random component
# Ideally random component ~ 1.0 for multiplicative model
errors_mult <- na.omit(as.numeric(random_mult))

cat("\n--- Error Analysis (Multiplicative Model) ---\n")
cat(sprintf("Mean of random component: %.4f (ideal: 1.0000)\n", mean(errors_mult)))
cat(sprintf("SD of random component:   %.4f\n", sd(errors_mult)))
cat(sprintf("Min: %.4f, Max: %.4f\n", min(errors_mult), max(errors_mult)))

# ---- Error plots ----
png("Task3_plots/08_error_analysis.png", width = 1000, height = 800)
par(mfrow = c(2, 2))

# Time plot of random component
plot(errors_mult, type = "o", pch = 16, cex = 0.5, col = "darkred",
     main = "Random Component (Multiplicative Model Errors)",
     ylab = "Random Factor", xlab = "Time Index")
abline(h = 1, col = "blue", lwd = 2, lty = 2)
grid(col = "gray85")

# Histogram of errors
hist(errors_mult, breaks = 25, col = "lightblue", border = "white",
     main = "Histogram of Random Component",
     xlab = "Random Factor", probability = TRUE)
lines(density(errors_mult), col = "red", lwd = 2)
curve(dnorm(x, mean = mean(errors_mult), sd = sd(errors_mult)),
      add = TRUE, col = "blue", lwd = 2, lty = 2)
abline(v = 1, col = "darkgreen", lwd = 2, lty = 2)

# ACF of errors
acf(errors_mult, lag.max = 24, main = "ACF of Random Component",
    ci.col = "steelblue", lwd = 2)

# Q-Q plot
qqnorm(errors_mult, main = "Q-Q Plot of Random Component")
qqline(errors_mult, col = "red", lwd = 2)

par(mfrow = c(1, 1))
dev.off()

# Normality tests on errors
sw_error <- shapiro.test(errors_mult)
cat(sprintf("\nShapiro-Wilk normality test: W = %.3f, p-value = %.4f\n",
            sw_error$statistic, sw_error$p.value))

ks_error <- ks.test(errors_mult, "pnorm",
                     mean = mean(errors_mult), sd = sd(errors_mult))
cat(sprintf("Kolmogorov-Smirnov test: D = %.3f, p-value = %.4f\n",
            ks_error$statistic, ks_error$p.value))

if (sw_error$p.value > 0.05) {
  cat(">>> Errors are normally distributed (do NOT reject H0).\n")
} else {
  cat(">>> Errors may NOT be normally distributed. Consider robust methods.\n")
}

# Check for remaining autocorrelation in errors
lb_error <- Box.test(errors_mult, lag = 12, type = "Ljung-Box")
cat(sprintf("Ljung-Box test (lag=12): X-squared = %.3f, p-value = %.4f\n",
            lb_error$statistic, lb_error$p.value))

###############################################################################
# 6. Model quality assessment
###############################################################################

cat("\n--- Model Quality Assessment ---\n")

# Reconstruct fitted values from multiplicative decomposition
# Fitted = Trend * Seasonal (ignoring random component)
fitted_mult <- trend_mult * seasonal_mult
fitted_mult_clean <- na.omit(fitted_mult)
actual_clean <- window(pass_ts, start = time(fitted_mult_clean)[1],
                       end = time(fitted_mult_clean)[length(fitted_mult_clean)])

actual_vals <- as.numeric(actual_clean)
fitted_vals <- as.numeric(fitted_mult_clean)

# Calculate quality metrics
mae  <- mean(abs(actual_vals - fitted_vals))
rmse <- sqrt(mean((actual_vals - fitted_vals)^2))
mape <- 100 * mean(abs((actual_vals - fitted_vals) / actual_vals))

cat(sprintf("MAE:  %.2f thousand passengers\n", mae))
cat(sprintf("RMSE: %.2f thousand passengers\n", rmse))
cat(sprintf("MAPE: %.2f%%\n", mape))

# ---- Plot original vs fitted ----
png("Task3_plots/09_original_vs_fitted.png", width = 1000, height = 550)
plot(pass_ts, main = "Original vs Fitted Values - Multiplicative Decomposition",
     ylab = "Passengers (thousands)", xlab = "Year",
     col = "black", lwd = 2, type = "o", pch = 16, cex = 0.5)
lines(fitted_mult, col = "red", lwd = 2.5, lty = 2)
grid(col = "gray85")
legend("topleft", legend = c("Original", "Fitted (Trend x Seasonal)"),
       col = c("black", "red"), lwd = 2, lty = c(1, 2))
dev.off()

# ---- Plot with error band ----
png("Task3_plots/10_residuals_vs_fitted.png", width = 900, height = 500)
plot(fitted_vals, actual_vals - fitted_vals,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values", ylab = "Residuals",
     col = "steelblue", pch = 16, cex = 0.8)
abline(h = 0, col = "red", lwd = 2, lty = 2)
grid(col = "gray85")
dev.off()

###############################################################################
# 7. Forecast 2 years ahead using decomposition
###############################################################################

cat("\n--- Forecasting 2 Years Ahead (24 months) ---\n")

# For multiplicative decomposition, we need to forecast trend and reuse seasonal factors
# Method: STL + ETS (more robust than classical decompose)

stl_fit <- stl(pass_ts, s.window = "periodic", t.window = 13)
stl_fcst <- forecast(stl_fit, method = "ets", h = 24, level = c(80, 95))

cat("\n2-Year Forecast (24 months):\n")
fcst_df <- data.frame(
  Month = as.character(time(stl_fcst$mean)),
  Forecast = round(as.numeric(stl_fcst$mean), 1),
  Lo80 = round(as.numeric(stl_fcst$lower[, 1]), 1),
  Hi80 = round(as.numeric(stl_fcst$upper[, 1]), 1),
  Lo95 = round(as.numeric(stl_fcst$lower[, 2]), 1),
  Hi95 = round(as.numeric(stl_fcst$upper[, 2]), 1)
)
print(fcst_df, row.names = FALSE)

# ---- Plot forecast ----
png("Task3_plots/11_forecast_2years.png", width = 1000, height = 600)
plot(stl_fcst, main = "2-Year Forecast - International Airline Passengers (Classical Decomposition + ETS)",
     ylab = "Passengers (thousands)", xlab = "Year",
     fcol = "tomato", flwd = 2.5, shaded = TRUE)
lines(fitted(stl_fit), col = "darkgreen", lwd = 1.5)
grid(col = "gray85")
legend("topleft",
       legend = c("Original", "Fitted (STL)", "Forecast", "95% CI"),
       col = c("black", "darkgreen", "tomato", "gray70"),
       lwd = c(1, 1.5, 2.5, 8))
dev.off()

# ---- Also show classical decomposition forecast approach ----
# Use the seasonal factors from decompose for the forecast
# Forecast trend using simple linear regression
t <- 1:length(trend_mult)
t_valid <- which(!is.na(trend_mult))
trend_lm <- lm(trend_mult[t_valid] ~ t[t_valid])
future_t <- (length(pass_ts) + 1):(length(pass_ts) + 24)
trend_fcst <- predict(trend_lm, newdata = data.frame(t = future_t))

# Seasonal factors repeat every 12 months
seasonal_idx <- (length(pass_ts) + 1 - 1) %% 12 + 1  # next month index
seasonal_fcst <- rep(seasonal_mult[1:12], length.out = 24)
# Align to correct starting month
start_month <- (length(pass_ts)) %% 12 + 1
if (start_month <= 12) {
  seasonal_fcst <- rep(seasonal_mult[start_month:(start_month + 11)], length.out = 24)
  # Handle wrap-around
}

# Actually, let's use a simpler approach for classical decomp forecast
seasonal_factors <- as.numeric(seasonal_mult[1:12])
start_idx <- (end(pass_ts)[2] %% 12) + 1
if (start_idx > 12) start_idx <- 1

# Rotate seasonal factors to align with forecast start
rotated_seasonal <- seasonal_factors[c(start_idx:12, 1:(start_idx - 1))]
seasonal_fcst_24 <- rep(rotated_seasonal, length.out = 24)
trend_fcst_24 <- predict(trend_lm, newdata = data.frame(t = future_t))
decomp_fcst_24 <- trend_fcst_24 * seasonal_fcst_24

cat("\nClassical Decomposition Forecast (24 months):\n")
decomp_fcst_df <- data.frame(
  Month = format(seq(as.Date("1961-01-01"), by = "month", length.out = 24), "%Y-%m"),
  Forecast = round(decomp_fcst_24, 1)
)
print(decomp_fcst_df, row.names = FALSE)

# ---- Plot classical decomposition forecast ----
png("Task3_plots/12_classical_decomp_forecast.png", width = 1000, height = 600)
plot(pass_ts, main = "Classical Decomposition Forecast - 2 Years Ahead",
     ylab = "Passengers (thousands)", xlab = "Year",
     col = "black", lwd = 2, type = "o", pch = 16, cex = 0.5,
     xlim = c(1949, 1963), ylim = c(100, max(decomp_fcst_24) * 1.1))
lines(fitted_mult, col = "darkgreen", lwd = 2, lty = 2)
lines(ts(decomp_fcst_24, start = c(1961, 1), frequency = 12),
      col = "tomato", lwd = 2.5, lty = 1, type = "o", pch = 17, cex = 0.8)
grid(col = "gray85")
legend("topleft",
       legend = c("Original", "Fitted (Trend x Seasonal)", "Forecast (Classical)"),
       col = c("black", "darkgreen", "tomato"),
       lwd = c(2, 2, 2.5), lty = c(1, 2, 1))
dev.off()

###############################################################################
# Save summary
###############################################################################

sink("Task3_Decomposition_summary.txt")
cat("TASK 3 - CLASSICAL DECOMPOSITION SUMMARY\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("1. Time Series Characteristics:\n")
cat("   Variable: International airline passengers (thousands/month)\n")
cat("   Period: Jan 1949 - Dec 1960 (144 months)\n")
cat("   Frequency: Monthly (s = 12)\n\n")

cat("2. Seasonality:\n")
cat("   Type: MULTIPLICATIVE\n")
cat("   Evidence: Increasing seasonal amplitude proportional to trend level\n")
cat("   Seasonal cycle: 12 months (annual)\n")
cat(sprintf("   Peak month: July-August (seasonal factor: %.3f)\n",
            max(seasonal_mult[1:12])))
cat(sprintf("   Trough month: November (seasonal factor: %.3f)\n\n",
            min(seasonal_mult[1:12])))

cat("3. Decomposition Model:\n")
cat("   Method: Classical multiplicative decomposition\n")
cat("   Y(t) = Trend(t) * Seasonal(t) * Random(t)\n\n")

cat("4. Model Quality:\n")
cat(sprintf("   MAE:  %.2f thousand passengers\n", mae))
cat(sprintf("   RMSE: %.2f thousand passengers\n", rmse))
cat(sprintf("   MAPE: %.2f%%\n\n", mape))

cat("5. Error Analysis:\n")
cat(sprintf("   Mean of random component: %.4f (target: 1.0)\n", mean(errors_mult)))
cat(sprintf("   Shapiro-Wilk normality: W = %.3f, p = %.4f\n",
            sw_error$statistic, sw_error$p.value))
cat(sprintf("   Ljung-Box (lag=12): X2 = %.3f, p = %.4f\n",
            lb_error$statistic, lb_error$p.value))

cat("\n6. Forecast (2 years / 24 months):\n")
cat("   First 6 forecast values:\n")
for (i in 1:6) {
  cat(sprintf("   %s: %.0f thousand passengers\n",
              decomp_fcst_df$Month[i], decomp_fcst_df$Forecast[i]))
}
sink()

cat("\nTask 3 complete. Results saved to Task3_plots/ and Task3_Decomposition_summary.txt\n")
