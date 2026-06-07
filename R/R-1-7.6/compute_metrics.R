setwd('C:/Users/swj17/.claude/projects/R-1-7.6')
df <- read.csv('morocco_dataset_clean.csv')

# Interpolate any remaining NAs
var_names <- c('GDP_growth','Inflation','FDI','Pop_growth','Exports','Capital_form','Govt_consumption','Trade','Unemployment')
for(v in var_names) {
  x <- df[[v]]
  if(any(is.na(x))) df[[v]] <- approx(seq_along(x), x, xout=seq_along(x), rule=2)$y
}

# Best forward model (from step 9)
model_best <- lm(GDP_growth ~ Pop_growth + Govt_consumption, data=df)

# Corrected model with lagged DV
df$GDP_lag1 <- c(NA, df$GDP_growth[-nrow(df)])
df_lag <- df[-1,]
model_lag <- lm(GDP_growth ~ GDP_lag1 + Pop_growth + Govt_consumption, data=df_lag)

# Full model
model_full <- lm(GDP_growth ~ Inflation + FDI + Pop_growth + Exports + Capital_form + Govt_consumption + Trade + Unemployment, data=df)

# ---- METRICS for BEST MODEL ----
pred_best <- predict(model_best, df)
actual <- df$GDP_growth
resid_best <- actual - pred_best

mape_best <- mean(abs(resid_best / pmax(abs(actual), 0.01))) * 100
mae_best <- mean(abs(resid_best))
rmse_best <- sqrt(mean(resid_best^2))

cat('=== BEST MODEL (Forward Selection: Pop_growth + Govt_consumption) ===\n')
cat(sprintf('MAPE:  %.2f%%\n', mape_best))
cat(sprintf('MAE:   %.4f\n', mae_best))
cat(sprintf('RMSE:  %.4f\n', rmse_best))
s <- summary(model_best)
cat(sprintf('R-squared:        %.4f\n', s$r.squared))
cat(sprintf('Adj R-squared:    %.4f\n', s$adj.r.squared))
cat(sprintf('AIC:              %.2f\n', AIC(model_best)))
cat(sprintf('BIC:              %.2f\n', BIC(model_best)))

# ---- METRICS for LAG-CORRECTED MODEL ----
pred_lag <- predict(model_lag, df_lag)
actual_lag <- df_lag$GDP_growth
resid_lag <- actual_lag - pred_lag

mape_lag <- mean(abs(resid_lag / pmax(abs(actual_lag), 0.01))) * 100
mae_lag <- mean(abs(resid_lag))
rmse_lag <- sqrt(mean(resid_lag^2))

cat('\n=== CORRECTED MODEL (with Lagged DV) ===\n')
cat(sprintf('MAPE:  %.2f%%\n', mape_lag))
cat(sprintf('MAE:   %.4f\n', mae_lag))
cat(sprintf('RMSE:  %.4f\n', rmse_lag))
s2 <- summary(model_lag)
cat(sprintf('R-squared:        %.4f\n', s2$r.squared))
cat(sprintf('Adj R-squared:    %.4f\n', s2$adj.r.squared))
cat(sprintf('AIC:              %.2f\n', AIC(model_lag)))

# ---- METRICS for FULL MODEL ----
pred_full <- predict(model_full, df)
resid_full <- actual - pred_full
mape_full <- mean(abs(resid_full / pmax(abs(actual), 0.01))) * 100
rmse_full <- sqrt(mean(resid_full^2))

cat('\n=== FULL MODEL (All 8 factors) ===\n')
cat(sprintf('MAPE:  %.2f%%\n', mape_full))
cat(sprintf('RMSE:  %.4f\n', rmse_full))
s3 <- summary(model_full)
cat(sprintf('R-squared:        %.4f\n', s3$r.squared))
cat(sprintf('Adj R-squared:    %.4f\n', s3$adj.r.squared))
cat(sprintf('AIC:              %.2f\n', AIC(model_full)))

# ---- FORECAST VALUES ----
means <- sapply(c('Pop_growth','Govt_consumption'), function(v) mean(df[[v]]))
sc10 <- as.data.frame(t(means * 1.10))
sc20 <- as.data.frame(t(means * 1.20))
p10 <- predict(model_best, newdata=sc10, interval='confidence', level=0.95)
p20 <- predict(model_best, newdata=sc20, interval='confidence', level=0.95)
p10pi <- predict(model_best, newdata=sc10, interval='prediction', level=0.95)
p20pi <- predict(model_best, newdata=sc20, interval='prediction', level=0.95)

cat('\n=== FORECAST ===\n')
cat(sprintf('10%% Growth: Point=%.4f, 95%% CI=[%.4f, %.4f], 95%% PI=[%.4f, %.4f]\n',
            p10[1], p10[2], p10[3], p10pi[2], p10pi[3]))
cat(sprintf('20%% Growth: Point=%.4f, 95%% CI=[%.4f, %.4f], 95%% PI=[%.4f, %.4f]\n',
            p20[1], p20[2], p20[3], p20pi[2], p20pi[3]))

# ---- DIAGNOSTICS ----
library(lmtest)
dw_best <- dwtest(model_best)
bp_best <- bptest(model_best)
library(tseries)
jb_best <- jarque.bera.test(resid_best)

cat('\n=== DIAGNOSTICS (Best Model) ===\n')
cat(sprintf('Durbin-Watson:    DW=%.4f  p=%.4f\n', dw_best$statistic, dw_best$p.value))
cat(sprintf('Breusch-Pagan:    BP=%.4f  p=%.4f\n', bp_best$statistic, bp_best$p.value))
cat(sprintf('Jarque-Bera:      JB=%.4f  p=%.4f\n', jb_best$statistic, jb_best$p.value))
