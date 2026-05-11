###############################################################################
# Practical Work 2 — Task 1: Exponential Smoothing
# Dataset: Johnson & Johnson Quarterly EPS (1960–1980)
# Non-seasonal quarterly data
###############################################################################

setwd("C:/Users/swj17/.claude/projects/R-analysis")
library(forecast)
library(ggplot2)

# ── 1. Data Preparation ─────────────────────────────────────────────────────

data(JohnsonJohnson)
y_all <- JohnsonJohnson  # length 84, quarterly 1960–1980

# Split: training (first 76 obs, 1960–1978 Q4) / control (last 8 obs, 1979–1980)
y_train <- window(y_all, end = c(1978, 4))
y_ctrl  <- window(y_all, start = c(1979, 1))

cat("Training set:", start(y_train), "–", end(y_train), "| n =", length(y_train), "\n")
cat("Control set: ", start(y_ctrl),  "–", end(y_ctrl),  "| n =", length(y_ctrl),  "\n\n")

# Plot
df_all <- data.frame(
  Year = as.numeric(time(y_all)),
  EPS  = as.numeric(y_all),
  Set  = ifelse(seq_along(y_all) > length(y_train), "Control", "Training")
)
ggplot(df_all, aes(x = Year, y = EPS, group = 1)) +
  geom_line(aes(colour = Set), linewidth = 1) +
  scale_colour_manual(values = c("Training" = "black", "Control" = "red")) +
  labs(title = "J&J Quarterly EPS (1960–1980)",
       subtitle = "Training (black) + Control (red)",
       x = "Year", y = "EPS ($)")
ggsave("task1_data_split.png", width = 8, height = 5)

# ── 2. Brown's Exponential Smoothing ─────────────────────────────────────────

# Brown's 1st-order (double ES) — linear trend, single parameter α
brown_linear <- function(y, alpha, h) {
  n <- length(y)
  S1 <- S2 <- numeric(n)
  S1[1] <- S2[1] <- y[1]
  for (t in 2:n) {
    S1[t] <- alpha * y[t] + (1 - alpha) * S1[t - 1]
    S2[t] <- alpha * S1[t] + (1 - alpha) * S2[t - 1]
  }
  a <- 2 * S1[n] - S2[n]
  b <- alpha / (1 - alpha) * (S1[n] - S2[n])
  fitted <- numeric(n)
  fitted[1] <- y[1]
  for (t in 2:n) {
    a_t <- 2 * S1[t - 1] - S2[t - 1]
    b_t <- alpha / (1 - alpha) * (S1[t - 1] - S2[t - 1])
    fitted[t] <- a_t + b_t
  }
  fcast <- a + b * seq_len(h)
  list(fitted = ts(fitted, start = start(y), frequency = frequency(y)),
       forecast = fcast, level = a, trend = b, alpha = alpha)
}

# Brown's 2nd-order (triple ES) — quadratic trend, single parameter α
brown_quadratic <- function(y, alpha, h) {
  n <- length(y)
  S1 <- S2 <- S3 <- numeric(n)
  S1[1] <- S2[1] <- S3[1] <- y[1]
  for (t in 2:n) {
    S1[t] <- alpha * y[t] + (1 - alpha) * S1[t - 1]
    S2[t] <- alpha * S1[t] + (1 - alpha) * S2[t - 1]
    S3[t] <- alpha * S2[t] + (1 - alpha) * S3[t - 1]
  }
  a <- 3 * S1[n] - 3 * S2[n] + S3[n]
  b <- alpha / (2 * (1 - alpha)^2) *
    ((6 - 5 * alpha) * S1[n] - 2 * (5 - 4 * alpha) * S2[n] + (4 - 3 * alpha) * S3[n])
  c <- alpha^2 / (1 - alpha)^2 * (S1[n] - 2 * S2[n] + S3[n])
  fitted <- numeric(n)
  fitted[1] <- y[1]
  for (t in 2:n) {
    at_ <- 3 * S1[t - 1] - 3 * S2[t - 1] + S3[t - 1]
    bt_ <- alpha / (2 * (1 - alpha)^2) *
      ((6 - 5 * alpha) * S1[t - 1] - 2 * (5 - 4 * alpha) * S2[t - 1] + (4 - 3 * alpha) * S3[t - 1])
    ct_ <- alpha^2 / (1 - alpha)^2 * (S1[t - 1] - 2 * S2[t - 1] + S3[t - 1])
    fitted[t] <- at_ + bt_ + 0.5 * ct_
  }
  fcast <- a + b * seq_len(h) + 0.5 * c * (seq_len(h))^2
  list(fitted = ts(fitted, start = start(y), frequency = frequency(y)),
       forecast = fcast, level = a, trend = b, quad = c, alpha = alpha)
}

# Optimise Brown's models → minimise MAD on control sample
alphas <- seq(0.01, 0.99, by = 0.01)
k <- length(y_ctrl)

# Brown linear: find best α by MAD
mad_bl <- sapply(alphas, function(a) {
  fit <- brown_linear(y_train, a, h = k)
  mean(abs(y_ctrl - fit$forecast))
})
alpha_bl <- alphas[which.min(mad_bl)]
cat(sprintf("Brown 1st-order optimal α = %.2f  (MAD = %.4f)\n", alpha_bl, min(mad_bl)))

# Brown quadratic: find best α by MAD
mad_bq <- sapply(alphas, function(a) {
  fit <- brown_quadratic(y_train, a, h = k)
  mean(abs(y_ctrl - fit$forecast))
})
alpha_bq <- alphas[which.min(mad_bq)]
cat(sprintf("Brown 2nd-order optimal α = %.2f  (MAD = %.4f)\n", alpha_bq, min(mad_bq)))

# Plot MAD curves
png("task1_brown_optimisation.png", width = 1000, height = 500)
par(mfrow = c(1, 2))
plot(alphas, mad_bl, type = "l", main = "Brown 1st-order: MAD vs α",
     xlab = "α", ylab = "MAD")
abline(v = alpha_bl, col = "red", lty = 2)
points(alpha_bl, min(mad_bl), col = "red", pch = 19)

plot(alphas, mad_bq, type = "l", main = "Brown 2nd-order: MAD vs α",
     xlab = "α", ylab = "MAD")
abline(v = alpha_bq, col = "red", lty = 2)
points(alpha_bq, min(mad_bq), col = "red", pch = 19)
dev.off()

# Fit final Brown models
fit_bl <- brown_linear(y_train, alpha_bl, h = k)
fit_bq <- brown_quadratic(y_train, alpha_bq, h = k)

# ── 3. Holt Model (1st-order polynomial) ─────────────────────────────────────

# Holt's linear trend — optimise α, β by minimum SSE on control
holt_lin <- function(y, alpha, beta, h) {
  n <- length(y)
  L <- T <- numeric(n)
  L[1] <- y[1]
  T[1] <- y[2] - y[1]
  for (t in 2:n) {
    L[t] <- alpha * y[t] + (1 - alpha) * (L[t - 1] + T[t - 1])
    T[t] <- beta * (L[t] - L[t - 1]) + (1 - beta) * T[t - 1]
  }
  fitted <- numeric(n)
  fitted[1] <- y[1]
  for (t in 2:n) fitted[t] <- L[t - 1] + T[t - 1]
  fcast <- L[n] + T[n] * seq_len(h)
  list(fitted = ts(fitted, start = start(y), frequency = frequency(y)),
       forecast = fcast, level = L[n], trend = T[n], alpha = alpha, beta = beta)
}

# Grid search for optimal (α, β)
grid <- expand.grid(alpha = seq(0.01, 0.99, by = 0.02),
                    beta  = seq(0.01, 0.99, by = 0.02))
sse_holt <- mapply(function(a, b) {
  fit <- holt_lin(y_train, a, b, h = k)
  sum((y_ctrl - fit$forecast)^2)
}, grid$alpha, grid$beta)
best <- which.min(sse_holt)
alpha_h <- grid$alpha[best]
beta_h  <- grid$beta[best]
cat(sprintf("Holt optimal — α = %.2f, β = %.2f  (SSE = %.4f)\n",
            alpha_h, beta_h, min(sse_holt)))

fit_holt <- holt_lin(y_train, alpha_h, beta_h, h = k)

# ── 4. Model Quality — MAPE ─────────────────────────────────────────────────

mape <- function(actual, pred) mean(abs((actual - pred) / actual)) * 100

cat("\n── MAPE Comparison ──\n")
mape_bl   <- mape(as.numeric(y_ctrl), fit_bl$forecast)
mape_bq   <- mape(as.numeric(y_ctrl), fit_bq$forecast)
mape_holt <- mape(as.numeric(y_ctrl), fit_holt$forecast)
cat(sprintf("Brown 1st-order: MAPE = %.2f%%\n", mape_bl))
cat(sprintf("Brown 2nd-order: MAPE = %.2f%%\n", mape_bq))
cat(sprintf("Holt:            MAPE = %.2f%%\n", mape_holt))

results <- data.frame(
  Model = c("Brown 1st-order", "Brown 2nd-order", "Holt"),
  Alpha = c(alpha_bl, alpha_bq, alpha_h),
  Beta  = c(NA, NA, beta_h),
  MAPE  = round(c(mape_bl, mape_bq, mape_holt), 2)
)
print(results)

# ── 5. Forecast 6 periods ahead using BEST model ─────────────────────────────

best_name <- results$Model[which.min(results$MAPE)]
cat(sprintf("\nBest model: %s (MAPE = %.2f%%)\n", best_name, min(results$MAPE)))

h_fcast <- 6
if (best_name == "Brown 1st-order") {
  fcast_final <- brown_linear(y_all, alpha_bl, h = h_fcast)
} else if (best_name == "Brown 2nd-order") {
  fcast_final <- brown_quadratic(y_all, alpha_bq, h = h_fcast)
} else {
  fcast_final <- holt_lin(y_all, alpha_h, beta_h, h = h_fcast)
}

fcast_ts <- ts(fcast_final$forecast,
               start = c(1981, 1), frequency = frequency(y_all))

cat("\n── Forecast for 6 periods ahead ──\n")
print(fcast_ts)

# Plot: full series + forecast
df_fcast <- data.frame(
  Year  = as.numeric(time(y_all)),
  EPS   = as.numeric(y_all),
  Type  = "Observed"
)
df_fcast <- rbind(df_fcast, data.frame(
  Year = as.numeric(time(fcast_ts)),
  EPS  = as.numeric(fcast_ts),
  Type = "Forecast"
))
ggplot(df_fcast, aes(x = Year, y = EPS, group = Type)) +
  geom_line(aes(colour = Type), linewidth = 1) +
  scale_colour_manual(values = c("Observed" = "black", "Forecast" = "blue")) +
  labs(title = paste("Forecast —", best_name),
       subtitle = "J&J Quarterly EPS with 6-period forecast",
       x = "Year", y = "EPS ($)")
ggsave("task1_forecast.png", width = 8, height = 5)

# Plot: training + control + fit
df_models <- data.frame(
  Year   = as.numeric(time(y_train)),
  Train  = as.numeric(y_train),
  Brown1 = as.numeric(fit_bl$fitted),
  Brown2 = as.numeric(fit_bq$fitted),
  Holt   = as.numeric(fit_holt$fitted)
)
df_ctrl_plot <- data.frame(
  Year    = as.numeric(time(y_ctrl)),
  Control = as.numeric(y_ctrl)
)
ggplot(df_models, aes(x = Year)) +
  geom_line(aes(y = Train,  colour = "Training"), linewidth = 0.8) +
  geom_line(data = df_ctrl_plot, aes(y = Control, colour = "Control"), linewidth = 1) +
  geom_line(aes(y = Brown1, colour = "Brown 1st"), linetype = "dashed", alpha = 0.7) +
  geom_line(aes(y = Brown2, colour = "Brown 2nd"), linetype = "dotted", alpha = 0.7) +
  geom_line(aes(y = Holt,   colour = "Holt"), linetype = "dotdash", alpha = 0.7) +
  scale_colour_manual(values = c("Training" = "black", "Control" = "darkgreen",
                                  "Brown 1st" = "blue", "Brown 2nd" = "orange",
                                  "Holt" = "red")) +
  labs(title = "Model fits vs actual data", x = "Year", y = "EPS ($)")
ggsave("task1_model_fits.png", width = 8, height = 5)

