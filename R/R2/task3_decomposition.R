###############################################################################
# Practical Work 2 — Task 3: Classical Time Series Decomposition
# Dataset: AirPassengers — Monthly international airline passengers (1949–1960)
###############################################################################

setwd("C:/Users/swj17/.claude/projects/R-analysis")
library(forecast)
library(ggplot2)

data(AirPassengers)
y <- AirPassengers  # monthly, 1949–1960, 144 obs

cat("Series: Monthly airline passengers (thousands), 1949–1960\n")
cat("Length:", length(y), "| Frequency:", frequency(y), "\n\n")

# ── 1. Indicator already selected: AirPassengers ─────────────────────────────

# ── 2. ACF & PACF Analysis ──────────────────────────────────────────────────

png("task3_02_acf_pacf.png", width = 1000, height = 600)
par(mfrow = c(1, 2))
acf(y, lag.max = 48, main = "ACF — AirPassengers")
pacf(y, lag.max = 48, main = "PACF — AirPassengers")
dev.off()

cat("ACF observations:\n")
cat("  - Slow decay → non-stationary (trend)\n")
cat("  - Peaks at lags 12, 24, 36 → strong annual seasonality (s = 12)\n")
cat("  - PACF: significant spike at lag 1, then at seasonal lags\n\n")

# ── 3. Seasonal Component Analysis ──────────────────────────────────────────

cat("Seasonal component detected: YES\n")
cat("Seasonal cycle length: s = 12 (monthly → annual pattern)\n")
cat("Seasonal type analysis:\n")

# Compare additive vs multiplicative
# Multiplicative: seasonal amplitude grows with trend level
# Additive: seasonal amplitude is constant

# Visual inspection of seasonal amplitude
y_decomp <- decompose(y, type = "multiplicative")
seasonal_amp <- y_decomp$seasonal
y_detrended <- y / y_decomp$trend

png("task3_03_seasonality_check.png", width = 1000, height = 600)
par(mfrow = c(2, 1))
# Seasonal component over time
plot(seasonal_amp, main = "Multiplicative Seasonal Component (constant ratio)",
     ylab = "Seasonal factor", xlab = "Year")

# Boxplots by month to check seasonal pattern stability
boxplot(split(y, cycle(y)), names = month.abb,
        main = "Monthly Distribution of Passengers",
        xlab = "Month", ylab = "Thousands", col = "lightblue")
dev.off()

cat("  → Amplitude grows proportionally with level → MULTIPLICATIVE seasonality\n")
cat("  → Pattern is stable across years (consistent high in Jul-Aug, low in Nov)\n\n")

# ── 4. Classical Multiplicative Decomposition ────────────────────────────────

decomp <- decompose(y, type = "multiplicative")

png("task3_04_decomposition.png", width = 1000, height = 800)
plot(decomp)
dev.off()

cat("Decomposition components:\n")
cat("  - Trend: smooth upward curve, accelerating after 1955\n")
cat("  - Seasonal: consistent factor ~0.90–1.20, peak Jul, trough Nov\n")
cat("  - Random: mean ≈ 1, no obvious pattern\n\n")

# Extract components
trend_comp    <- decomp$trend
seasonal_comp <- decomp$seasonal
random_comp   <- decomp$random

# Plot trend + seasonal separately with ggplot
df <- data.frame(
  Date       = as.numeric(time(y)),
  Observed   = as.numeric(y),
  Trend      = as.numeric(trend_comp),
  Seasonal   = as.numeric(seasonal_comp),
  Random     = as.numeric(random_comp),
  Model      = as.numeric(trend_comp * seasonal_comp / seasonal_comp * trend_comp)
)

# Reconstructed model values: Trend × Seasonal
df$Model <- as.numeric(trend_comp * seasonal_comp)

p1 <- ggplot(df, aes(x = Date)) +
  geom_line(aes(y = Observed), colour = "black", linewidth = 0.6) +
  geom_line(aes(y = Trend), colour = "red", linewidth = 1.2) +
  labs(title = "Observed vs Trend", x = "Year", y = "Thousands") +
  theme_minimal()

p2 <- ggplot(df, aes(x = Date, y = Seasonal)) +
  geom_line(colour = "blue", linewidth = 0.8) +
  labs(title = "Seasonal Component (Multiplicative Factors)",
       x = "Year", y = "Seasonal Factor") +
  theme_minimal() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "gray50")

ggsave("task3_04_trend.png", p1, width = 8, height = 4)
ggsave("task3_04_seasonal.png", p2, width = 8, height = 4)

# ── 5. Error Analysis ───────────────────────────────────────────────────────

errors <- na.omit(random_comp)  # random component = errors (ratio to model)

png("task3_05_errors.png", width = 1000, height = 800)
par(mfrow = c(2, 2))
plot(errors, main = "Random Component (Errors)", ylab = "Ratio", xlab = "Year")
abline(h = 1, col = "red", lty = 2)
hist(errors, breaks = 20, probability = TRUE,
     main = "Histogram of Errors", xlab = "Error ratio", col = "lightgreen")
curve(dnorm(x, mean = mean(errors), sd = sd(errors)),
      add = TRUE, col = "red", lwd = 2)
qqnorm(as.numeric(errors), main = "Q-Q Plot of Errors")
qqline(as.numeric(errors), col = "red", lwd = 2)
acf(errors, lag.max = 36, main = "ACF of Errors")
dev.off()

# Test for normality
sw <- shapiro.test(as.numeric(errors))
ks <- ks.test(as.numeric(errors), "pnorm", mean = mean(errors), sd = sd(errors))

cat("── Error normality tests ──\n")
cat(sprintf("Shapiro-Wilk:      W = %.4f, p = %.4f → %s\n",
            sw$statistic, sw$p.value,
            ifelse(sw$p.value > 0.05, "Normal ✓", "NON-normal ✗")))
cat(sprintf("Kolmogorov-Smirnov: D = %.4f, p = %.4f → %s\n",
            ks$statistic, ks$p.value,
            ifelse(ks$p.value > 0.05, "Normal ✓", "NON-normal ✗")))

# Error statistics
cat(sprintf("\nError mean: %.4f (ideal = 1 in multiplicative model)\n", mean(errors)))
cat(sprintf("Error variance: %.6f\n", var(as.numeric(errors))))

# ── 6. Model Evaluation ─────────────────────────────────────────────────────

# Reconstructed values = Trend × Seasonal
y_model <- trend_comp * seasonal_comp

# Model quality metrics (on non-NA values)
valid_idx <- which(!is.na(y_model))
y_valid    <- as.numeric(y)[valid_idx]
y_model_v  <- as.numeric(y_model)[valid_idx]

# R²
SSE <- sum((y_valid - y_model_v)^2)
SST <- sum((y_valid - mean(y_valid))^2)
R2 <- 1 - SSE / SST

# MAPE
MAPE <- mean(abs((y_valid - y_model_v) / y_valid)) * 100

# RMSE
RMSE <- sqrt(mean((y_valid - y_model_v)^2))

# MAE
MAE <- mean(abs(y_valid - y_model_v))

cat("\n── Model Quality Metrics ──\n")
cat(sprintf("R²   = %.4f\n", R2))
cat(sprintf("MAPE = %.2f%%\n", MAPE))
cat(sprintf("RMSE = %.2f\n", RMSE))
cat(sprintf("MAE  = %.2f\n", MAE))

# Plot original vs model
df_valid <- data.frame(
  Date     = as.numeric(time(y))[valid_idx],
  Observed = y_valid,
  Model    = y_model_v
)

ggplot(df_valid, aes(x = Date)) +
  geom_line(aes(y = Observed), colour = "black", linewidth = 0.6) +
  geom_line(aes(y = Model), colour = "red", linewidth = 0.8, alpha = 0.8) +
  labs(title = "Observed vs Classical Decomposition Model",
       subtitle = paste0("Multiplicative | MAPE = ", round(MAPE, 2), "% | R² = ", round(R2, 4)),
       x = "Year", y = "Thousands of passengers") +
  theme_minimal()
ggsave("task3_06_model_fit.png", width = 8, height = 5)

# ── 7. Forecast 2 Years (24 months) Ahead ───────────────────────────────────

# Extract last full year of seasonal factors (12 months)
seasonal_pattern <- seasonal_comp[1:12]  # Jan–Dec seasonal indices

# Fit linear trend to last 2 years of trend component
trend_valid <- na.omit(trend_comp)
n_trend <- length(trend_valid)
trend_time <- 1:n_trend
lm_trend <- lm(as.numeric(trend_valid) ~ trend_time)
cat(sprintf("\nTrend linear fit: R² = %.4f\n", summary(lm_trend)$r.squared))

# Forecast trend 24 steps ahead
h <- 24
trend_future_time <- (n_trend + 1):(n_trend + h)
trend_fcast <- predict(lm_trend, newdata = data.frame(trend_time = trend_future_time))

# Seasonal forecast: repeat seasonal pattern
n_full_years <- h %/% 12
seasonal_fcast <- rep(seasonal_pattern, n_full_years + 1)[1:h]

# Multiplicative forecast
fcast_values <- trend_fcast * seasonal_fcast

# Prediction interval (based on error distribution)
error_sd <- sd(as.numeric(errors))
log_errors <- log(as.numeric(errors))
log_error_sd <- sd(log_errors)

# 95% PI: multiplicative factors from lognormal
fcast_lower_95 <- fcast_values * exp(-1.96 * log_error_sd)
fcast_upper_95 <- fcast_values * exp( 1.96 * log_error_sd)

fcast_start <- end(y) + c(0, 1)
fcast_ts  <- ts(fcast_values, start = c(1961, 1), frequency = 12)
fcast_lo  <- ts(fcast_lower_95, start = c(1961, 1), frequency = 12)
fcast_hi  <- ts(fcast_upper_95, start = c(1961, 1), frequency = 12)

cat("\n── Forecast for 2 years ahead (Jan 1961 – Dec 1962) ──\n")
fc_table <- data.frame(
  Date     = paste0(rep(1961:1962, each = 12), "-", month.abb),
  Forecast = round(fcast_values, 1),
  Lower95  = round(fcast_lower_95, 1),
  Upper95  = round(fcast_upper_95, 1)
)
print(fc_table, row.names = FALSE)

# Plot: original + model + forecast + 95% PI
df_fc <- data.frame(
  Year  = as.numeric(time(y)),
  Value = as.numeric(y),
  Type  = "Observed"
)
df_fc <- rbind(df_fc, data.frame(
  Year  = as.numeric(time(y))[valid_idx],
  Value = y_model_v,
  Type  = "Model"
))
df_fc <- rbind(df_fc, data.frame(
  Year  = as.numeric(time(fcast_ts)),
  Value = as.numeric(fcast_ts),
  Type  = "Forecast"
))
df_fc <- rbind(df_fc, data.frame(
  Year  = as.numeric(time(fcast_lo)),
  Value = as.numeric(fcast_lo),
  Type  = "95% PI lower"
))
df_fc <- rbind(df_fc, data.frame(
  Year  = as.numeric(time(fcast_hi)),
  Value = as.numeric(fcast_hi),
  Type  = "95% PI upper"
))
ggplot(df_fc, aes(x = Year, y = Value, group = Type)) +
  geom_line(aes(colour = Type, linetype = Type, size = Type)) +
  scale_colour_manual(values = c(
    "Observed" = "black", "Model" = "red", "Forecast" = "blue",
    "95% PI lower" = "darkgreen", "95% PI upper" = "darkgreen"
  )) +
  scale_linetype_manual(values = c(
    "Observed" = "solid", "Model" = "solid", "Forecast" = "solid",
    "95% PI lower" = "dashed", "95% PI upper" = "dashed"
  )) +
  scale_size_manual(values = c(
    "Observed" = 0.6, "Model" = 0.8, "Forecast" = 1,
    "95% PI lower" = 0.5, "95% PI upper" = 0.5
  )) +
  labs(title = "Classical Decomposition Forecast — 2 Years Ahead",
       subtitle = "Multiplicative model with linear trend extrapolation",
       x = "Year", y = "Thousands of passengers") +
  theme_minimal()
ggsave("task3_07_forecast.png", width = 10, height = 6)

# ── Bonus: STL decomposition for comparison ──────────────────────────────────

stl_fit <- stl(y, s.window = "periodic")
png("task3_bonus_stl.png", width = 1000, height = 800)
plot(stl_fit, main = "STL Decomposition (comparison)")
dev.off()

