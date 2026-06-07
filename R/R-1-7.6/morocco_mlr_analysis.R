###############################################################################
# Practice 1: Multiple Linear Regression — Morocco GDP Growth Forecasting
# Country: Morocco (MA) | Data: World Bank WDI | Period: 1967–2024 (58 obs)
# Dependent Variable: GDP_growth (% annual)
# Independent Variables (8): GFCF_growth, Export_growth, Import_growth,
#   Pop_growth, Agric_growth, Inflation, Unemployment, FDI_pct_GDP
###############################################################################

# ===========================================================================
# STEP 0: Setup
# ===========================================================================
required_packages <- c("lmtest", "car", "tseries", "moments", "sandwich",
                        "dplyr", "ggplot2", "corrplot", "nortest")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos="https://cran.r-project.org")
  library(pkg, character.only = TRUE)
}
setwd("C:/Users/swj17/.claude/projects/R-1-7.6")

# ===========================================================================
# STEP 1: Load & prepare data
# ===========================================================================
cat("\n", paste(rep("=", 72), collapse=""), "\n")
cat("  MLR ANALYSIS — MOROCCO GDP GROWTH (V3: Agriculture + Growth Drivers)\n")
cat(paste(rep("=", 72), collapse=""), "\n")

cat("\n========== STEP 1: DATA COLLECTION ==========\n")
df <- read.csv("morocco_dataset_v3.csv", stringsAsFactors = FALSE)
year_vec <- df$Year; df <- df[, !names(df) %in% "Year"]
for(j in 1:ncol(df)) df[[j]] <- as.numeric(df[[j]])

# Interpolate NAs
for(v in names(df)) {
  x <- df[[v]]
  nna <- sum(is.na(x))
  if(nna > 0) { df[[v]] <- approx(seq_along(x), x, xout=seq_along(x), rule=2)$y
                cat(sprintf("  %-16s %d NAs interpolated\n", v, nna)) }
}

dep_var  <- "GDP_growth"
indep_vars <- setdiff(names(df), dep_var)
n_obs <- nrow(df)

cat(sprintf("\nFinal dataset: %d obs, %d vars | Years: %d–%d\n",
            n_obs, ncol(df), min(year_vec), max(year_vec)))
cat(sprintf("DV: %s\n", dep_var))
cat(sprintf("IVs (%d): %s\n", length(indep_vars), paste(indep_vars, collapse=", ")))
write.csv(cbind(Year=year_vec, df), "morocco_dataset_final.csv", row.names=FALSE)

# ===========================================================================
# STEP 2: Descriptive Statistics
# ===========================================================================
cat("\n========== STEP 2: DESCRIPTIVE STATISTICS ==========\n")

calc_mode <- function(x) { u <- unique(x); u[which.max(tabulate(match(x, u)))] }
all_vars <- names(df)
stats_list <- lapply(all_vars, function(v) {
  x <- df[[v]]
  data.frame(Variable=v, N=length(x), Mean=round(mean(x),4), Variance=round(var(x),4),
             Std_Dev=round(sd(x),4), Mode=round(calc_mode(x),4), Median=round(median(x),4),
             Min=round(min(x),4), Max=round(max(x),4), Range=round(max(x)-min(x),4),
             Kurtosis=round(kurtosis(x),4), Skewness=round(skewness(x),4),
             CV_pct=round((sd(x)/mean(x))*100,2),
             Q1=round(quantile(x,0.25),4), Q3=round(quantile(x,0.75),4), IQR=round(IQR(x),4))
})
stats_table <- do.call(rbind, stats_list)
print(stats_table, row.names=FALSE)
write.csv(stats_table, "descriptive_statistics.csv", row.names=FALSE)

# ===========================================================================
# STEP 3: Correlation Matrix
# ===========================================================================
cat("\n========== STEP 3: CORRELATION MATRIX ==========\n")
cor_matrix <- cor(df)
cat("Correlation Matrix:\n"); print(round(cor_matrix, 4))
write.csv(round(cor_matrix,4), "correlation_matrix.csv")

png("correlation_plot.png", width=1200, height=1000, res=130)
corrplot(cor_matrix, method="color", type="upper", addCoef.col="black",
         tl.col="black", tl.cex=0.6, number.cex=0.55,
         title="Correlation Matrix — Morocco GDP Growth (V3)", mar=c(0,0,2,0))
dev.off()

cat("\nGDP_growth vs each factor:\n")
for(v in indep_vars) {
  r <- cor(df$GDP_growth, df[[v]])
  s <- if(abs(r)>0.7) "STRONG" else if(abs(r)>0.4) "moderate" else "weak"
  d <- ifelse(r>0,"+","")
  cat(sprintf("  %-18s r = %+.4f  (%s %s)\n", v, r, s, ifelse(r>0,"positive","negative")))
}

# ===========================================================================
# STEP 4: Full MLR Model (All 8 Factors)
# ===========================================================================
cat("\n========== STEP 4: FULL MLR MODEL (ALL 8 FACTORS) ==========\n")
formula_full <- as.formula(paste(dep_var, "~", paste(indep_vars, collapse=" + ")))
model_full <- lm(formula_full, data=df)
summary_full <- summary(model_full)
cat("Full Model:\n"); print(summary_full)
cat(sprintf("Residual SE: %.4f on %d df\n", summary_full$sigma, summary_full$df[2]))

# ===========================================================================
# STEP 5: Coefficient Interpretation
# ===========================================================================
cat("\n========== STEP 5: COEFFICIENT INTERPRETATION ==========\n")
coeffs <- coef(model_full)
cat(sprintf("Intercept (β₀=%.4f): Expected GDP_growth when all IVs=0.\n\n", coeffs[1]))
for(i in 2:length(coeffs)) {
  vn <- names(coeffs)[i]; cv <- coeffs[i]
  dw <- ifelse(cv>0, "increases", "decreases")
  cat(sprintf("  %-18s β=%+.4f  → 1-unit ↑ in %s %s GDP_growth by %.4f pp\n",
              vn, cv, vn, dw, abs(cv)))
}

# ===========================================================================
# STEP 6: Quality & Reliability
# ===========================================================================
cat("\n========== STEP 6: STANDARD ERRORS, SIGNIFICANCE, CI ==========\n")
coef_tab <- summary_full$coefficients
cat("Coefficients with SE, t, p:\n"); print(round(coef_tab, 4))

cat("\nSignificance (α=0.05):\n")
for(i in 1:nrow(coef_tab)) {
  pv <- coef_tab[i,"Pr(>|t|)"]
  sig <- if(pv<0.001) "*** (0.1%)" else if(pv<0.01) "** (1%)" else
         if(pv<0.05) "* (5%)" else if(pv<0.10) ". (10%)" else "— ns"
  cat(sprintf("  %-18s p=%.4f  %s\n", rownames(coef_tab)[i], pv, sig))
}
cat("\n95% Confidence Intervals:\n"); print(round(confint(model_full, level=0.95), 4))

# ===========================================================================
# STEP 7: R² and Adjusted R²
# ===========================================================================
cat("\n========== STEP 7: R² AND ADJUSTED R² ==========\n")
r2 <- summary_full$r.squared; ar2 <- summary_full$adj.r.squared
cat(sprintf("R-squared:          %.4f (%.2f%%)\n", r2, r2*100))
cat(sprintf("Adjusted R-squared: %.4f (%.2f%%)\n", ar2, ar2*100))
cat(sprintf("Difference:         %.4f\n", r2 - ar2))
cat(sprintf("\n→ %.1f%% of GDP growth variation explained by the 8 factors.\n", ar2*100))
cat("→ Adjusted R² penalizes irrelevant regressors; always ≤ R².\n")

# ===========================================================================
# STEP 8: Multicollinearity (VIF)
# ===========================================================================
cat("\n========== STEP 8: MULTICOLLINEARITY (VIF) ==========\n")
vif_vals <- vif(model_full)
for(i in seq_along(vif_vals)) {
  vv <- vif_vals[i]
  lvl <- if(vv>10) "*** SERIOUS" else if(vv>5) "** MODERATE" else "OK"
  cat(sprintf("  %-18s VIF=%7.2f  %s\n", names(vif_vals)[i], vv, lvl))
}

# ===========================================================================
# STEP 9: Forward Selection (Adjusted R² Growth)
# ===========================================================================
cat("\n========== STEP 9: FORWARD SELECTION (Δ Adj.R² > 0) ==========\n")
cat("Method: Add variable maximizing Adj.R² at each step.\n")
cat("Stop when no variable produces positive Δ Adj.R².\n\n")

selected <- character(0); remaining <- indep_vars; best_adj_r2 <- 0; step <- 0
while(length(remaining) > 0) {
  best_v <- NA; best_score <- -Inf; best_r2v <- NA; best_aicv <- NA
  for(v in remaining) {
    fmla <- as.formula(paste("GDP_growth ~", paste(c(selected, v), collapse=" + ")))
    m <- lm(fmla, data=df); s <- summary(m)
    ar2v <- s$adj.r.squared
    if(ar2v > best_score) { best_score <- ar2v; best_v <- v; best_r2v <- s$r.squared; best_aicv <- AIC(m) }
  }
  if(best_score > best_adj_r2) {
    step <- step + 1; delta <- best_score - best_adj_r2; best_adj_r2 <- best_score
    selected <- c(selected, best_v); remaining <- setdiff(remaining, best_v)
    cat(sprintf("  Step %d: + %-18s Adj.R²=%.4f (+%.4f) R²=%.4f AIC=%.1f\n",
                step, best_v, best_score, delta, best_r2v, best_aicv))
  } else { cat("\n  No improvement. Stopping.\n"); break }
}

cat(sprintf("\n=== Best model: %d variables ===\n", length(selected)))
cat(sprintf("Selected: %s\n", paste(selected, collapse=", ")))

best_formula <- as.formula(paste("GDP_growth ~", paste(selected, collapse=" + ")))
model_best <- lm(best_formula, data=df)
summary_best <- summary(model_best)
cat("\nBest Model Summary:\n"); print(summary_best)

# ===========================================================================
# STEP 10: Significance of Best Model
# ===========================================================================
cat("\n========== STEP 10: SIGNIFICANCE (α=0.05) ==========\n")
coef_best <- summary_best$coefficients
cat("Coefficients:\n")
for(i in 1:nrow(coef_best)) {
  pv <- coef_best[i,"Pr(>|t|)"]
  sig <- if(pv<0.001) "***" else if(pv<0.01) "**" else if(pv<0.05) "*" else if(pv<0.10) "." else "ns"
  cat(sprintf("  %-18s Est=%+.4f  p=%.4f  %s\n", rownames(coef_best)[i], coef_best[i,"Estimate"], pv, sig))
}
fstat <- summary_best$fstatistic
fpval <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
cat(sprintf("\nF-stat = %.4f on %d,%d df; p = %.2e\n", fstat[1], fstat[2], fstat[3], fpval))
cat(ifelse(fpval<0.05, "→ Model IS statistically significant.\n", "→ NOT significant.\n"))
cat("\n95% CI:\n"); print(round(confint(model_best, level=0.95), 4))

# ===========================================================================
# STEP 11: RESET Test
# ===========================================================================
cat("\n========== STEP 11: RESET TEST ==========\n")
reset_r <- resettest(model_best, power=2:3, type="fitted")
cat(sprintf("Ramsey RESET: F=%.4f p=%.4f  %s\n",
            reset_r$statistic, reset_r$p.value,
            ifelse(reset_r$p.value>0.05, "→ OK, linear adequate.", "→ WARNING: misspecification.")))

# ===========================================================================
# STEP 12: Graphical Homoskedasticity
# ===========================================================================
cat("\n========== STEP 12: GRAPHICAL HOMOSKEDASTICITY ==========\n")
resid_b <- residuals(model_best); fitted_b <- fitted(model_best); std_b <- rstandard(model_best)
png("homoskedasticity_plots.png", width=1500, height=500, res=130)
par(mfrow=c(1,3))
plot(fitted_b, resid_b, pch=19, col=rgb(0,0,0.5,0.4), main="Residuals vs Fitted",
     xlab="Fitted", ylab="Residuals"); abline(h=0, col="red", lwd=2, lty=2)
lines(lowess(fitted_b, resid_b), col="blue", lwd=2)
plot(fitted_b, sqrt(abs(std_b)), pch=19, col=rgb(0,0,0.5,0.4), main="Scale-Location",
     xlab="Fitted", ylab=expression(sqrt("|Std Resid|")))
lines(lowess(fitted_b, sqrt(abs(std_b))), col="blue", lwd=2)
plot(fitted_b, std_b, pch=19, col=rgb(0,0,0.5,0.4), main="Std Residuals vs Fitted",
     xlab="Fitted", ylab="Std Residuals"); abline(h=c(-2,0,2), col=c("orange","red","orange"), lty=c(3,2,3))
par(mfrow=c(1,1)); dev.off(); cat("Saved: homoskedasticity_plots.png\n")

# ===========================================================================
# STEP 13: Homoskedasticity Tests
# ===========================================================================
cat("\n========== STEP 13: HOMOSKEDASTICITY TESTS ==========\n")
bp <- bptest(model_best, studentize=TRUE)
gq <- gqtest(model_best, fraction=0.2, order.by=fitted(model_best))
cat(sprintf("Breusch-Pagan: BP=%.4f df=%d p=%.4f  %s\n",
            bp$statistic, bp$parameter, bp$p.value,
            ifelse(bp$p.value>0.05, "→ Homoskedasticity ✓", "→ HETEROSKEDASTICITY ✗")))
cat(sprintf("Goldfeld-Quandt: GQ=%.4f p=%.4f  %s\n",
            gq$statistic, gq$p.value,
            ifelse(gq$p.value>0.05, "→ Homoskedasticity ✓", "→ HETEROSKEDASTICITY ✗")))

# ===========================================================================
# STEP 14: Graphical Autocorrelation
# ===========================================================================
cat("\n========== STEP 14: GRAPHICAL AUTOCORRELATION ==========\n")
n_res <- length(resid_b)
png("autocorrelation_plots.png", width=1500, height=1000, res=130)
par(mfrow=c(2,2))
plot(year_vec, resid_b, type="b", pch=19, col="darkblue", main="Residuals Over Time",
     xlab="Year", ylab="Residuals"); abline(h=0, col="red", lwd=2, lty=2)
plot(resid_b[-n_res], resid_b[-1], pch=19, col=rgb(0,0,0.5,0.4),
     xlab=expression(e[t-1]), ylab=expression(e[t]), main="Lag-1 Residuals")
abline(h=0, v=0, col="red", lty=2)
acf(resid_b, main="ACF", lag.max=20)
pacf(resid_b, main="PACF", lag.max=20)
par(mfrow=c(1,1)); dev.off(); cat("Saved: autocorrelation_plots.png\n")

# ===========================================================================
# STEP 15: Durbin-Watson & Breusch-Godfrey
# ===========================================================================
cat("\n========== STEP 15: AUTOCORRELATION TESTS ==========\n")
dw <- dwtest(model_best, alternative="two.sided")
cat(sprintf("Durbin-Watson: DW=%.4f p=%.4f  %s\n",
            dw$statistic, dw$p.value,
            ifelse(dw$p.value>0.05, "→ No autocorrelation ✓", "→ AUTOCORRELATION ✗")))
cat("Breusch-Godfrey LM:\n")
for(k in 1:3) {
  bg <- bgtest(model_best, order=k)
  cat(sprintf("  Order %d: LM=%.4f df=%d p=%.4f  %s\n",
              k, bg$statistic, bg$parameter, bg$p.value,
              ifelse(bg$p.value>0.05, "OK", "AUTOCORR!")))
}

# ===========================================================================
# STEP 16: Normality of Residuals
# ===========================================================================
cat("\n========== STEP 16: NORMALITY OF RESIDUALS ==========\n")
png("normality_plots.png", width=1500, height=550, res=130)
par(mfrow=c(1,2))
hist(resid_b, breaks=14, prob=TRUE, col="lightblue", main="Histogram of Residuals",
     xlab="Residuals"); curve(dnorm(x, mean(resid_b), sd(resid_b)), add=TRUE, col="darkred", lwd=2)
qqnorm(resid_b, main="Normal Q-Q Plot", pch=19, col="darkblue"); qqline(resid_b, col="red", lwd=2)
par(mfrow=c(1,1)); dev.off(); cat("Saved: normality_plots.png\n")

sw <- shapiro.test(resid_b)
jb <- jarque.bera.test(resid_b)
ad <- ad.test(resid_b)
cat(sprintf("Shapiro-Wilk:     W=%.4f p=%.4f\n", sw$statistic, sw$p.value))
cat(sprintf("Jarque-Bera:      JB=%.4f p=%.4f\n", jb$statistic, jb$p.value))
cat(sprintf("Anderson-Darling: A=%.4f p=%.4f\n", ad$statistic, ad$p.value))
cat(sprintf("Skewness=%.4f  Kurtosis=%.4f\n", skewness(resid_b), kurtosis(resid_b)))
cat(ifelse(sw$p.value>0.05, "→ Residuals ARE normal ✓\n", "→ Residuals NOT normal ✗\n"))

# ===========================================================================
# STEP 17: Gauss-Markov Corrections
# ===========================================================================
cat("\n========== STEP 17: GAUSS-MARKOV CORRECTIONS ==========\n")
n_hetero <- bp$p.value < 0.05; n_auto <- dw$p.value < 0.05; n_spec <- reset_r$p.value < 0.05
cat("Diagnosis:\n")
cat(sprintf("  Heteroskedasticity: %s\n", ifelse(n_hetero,"VIOLATED ✗","OK ✓")))
cat(sprintf("  Autocorrelation:    %s\n", ifelse(n_auto,"VIOLATED ✗","OK ✓")))
cat(sprintf("  Misspecification:   %s\n", ifelse(n_spec,"VIOLATED ✗","OK ✓")))

corrections <- c(); model_final <- model_best

if(n_hetero || n_auto) {
  cat("\n→ HAC robust SE (Newey-West):\n")
  print(round(coeftest(model_final, vcov=vcovHAC(model_final)), 4))
  corrections <- c(corrections, "HAC robust standard errors (Newey-West)")
}
if(n_auto) {
  cat("\n→ Adding lagged GDP_growth to fix autocorrelation:\n")
  df$GDP_lag1 <- c(NA, df$GDP_growth[-nrow(df)]); df_lag <- df[-1,]
  lag_fmla <- update(formula(model_final), . ~ GDP_lag1 + .)
  m_lag <- lm(lag_fmla, data=df_lag)
  dw_lag <- dwtest(m_lag)
  cat(sprintf("  DW after lag: %.4f p=%.4f\n", dw_lag$statistic, dw_lag$p.value))
  if(dw_lag$p.value > 0.05) {
    cat("  Autocorrelation resolved.\n")
    corrections <- c(corrections, "Lagged dependent variable added")
    model_final <- m_lag; df_final <- df_lag
  }
}

# Outlier check
cooks <- cooks.distance(model_final); thresh <- 4/length(cooks); bad <- which(cooks > thresh)
if(length(bad) > 0) {
  cat(sprintf("\n→ %d influential obs (Cook's D > %.4f):\n", length(bad), thresh))
  for(o in bad) cat(sprintf("   Obs %d: CookD=%.4f Year=%d\n", o, cooks[o], year_vec[o]))
}

cat("\nCorrections applied:\n")
if(length(corrections) > 0) { for(i in seq_along(corrections)) cat(sprintf("  %d. %s\n", i, corrections[i]))
} else { cat("  None needed. All Gauss-Markov conditions satisfied. ✓\n") }

cat("\n--- Final Model ---\n"); print(summary(model_final))

# ===========================================================================
# STEP 18: Point & Interval Forecasts
# ===========================================================================
cat("\n========== STEP 18: FORECASTS ==========\n")
fc_vars <- selected; means_fc <- sapply(fc_vars, function(v) mean(df[[v]]))
cat("Mean values of selected factors:\n")
for(v in fc_vars) cat(sprintf("  %-18s = %.4f\n", v, means_fc[v]))

sc10 <- as.data.frame(t(means_fc * 1.10))
sc20 <- as.data.frame(t(means_fc * 1.20))

for(pct in c(10, 20)) {
  sc <- if(pct==10) sc10 else sc20
  pc <- predict(model_best, newdata=sc, interval="confidence", level=0.95)
  pp <- predict(model_best, newdata=sc, interval="prediction", level=0.95)
  cat(sprintf("\n%d%% Growth Scenario:\n", pct))
  cat(sprintf("  Point forecast:             %.4f%%\n", pc[1,"fit"]))
  cat(sprintf("  95%% Confidence interval:   [%.4f%%, %.4f%%]\n", pc[1,"lwr"], pc[1,"upr"]))
  cat(sprintf("  95%% Prediction interval:   [%.4f%%, %.4f%%]\n", pp[1,"lwr"], pp[1,"upr"]))
}

cat("\n--- Per-Factor Sensitivity (+20%) ---\n")
for(v in fc_vars) {
  nd <- as.data.frame(t(means_fc)); nd[[v]] <- means_fc[v] * 1.20
  p <- predict(model_best, newdata=nd, interval="confidence", level=0.95)
  cat(sprintf("  %s +20%%: GDP_growth = %.4f%%\n", v, p[1,"fit"]))
}

# ===========================================================================
# MODEL PERFORMANCE METRICS
# ===========================================================================
cat("\n========== PERFORMANCE METRICS ==========\n")
pred_b <- predict(model_best, df); actual <- df$GDP_growth; resid_bb <- actual - pred_b
mape_b <- mean(abs(resid_bb / pmax(abs(actual), 0.01))) * 100
mae_b  <- mean(abs(resid_bb)); rmse_b <- sqrt(mean(resid_bb^2))

pred_f <- predict(model_full, df); resid_ff <- actual - pred_f
mape_f <- mean(abs(resid_ff / pmax(abs(actual), 0.01))) * 100
rmse_f <- sqrt(mean(resid_ff^2))

cat(sprintf("\n  Metric        Full(8var)      Best(%dvar)\n", length(fc_vars)))
cat(sprintf("  ─────────     ──────────      ──────────\n"))
cat(sprintf("  Adj R²        %.4f           %.4f\n", summary_full$adj.r.squared, summary_best$adj.r.squared))
cat(sprintf("  R²            %.4f           %.4f\n", summary_full$r.squared, summary_best$r.squared))
cat(sprintf("  RMSE          %.4f          %.4f\n", rmse_f, rmse_b))
cat(sprintf("  MAE           —              %.4f\n", mae_b))
cat(sprintf("  MAPE          %.2f%%          %.2f%%\n", mape_f, mape_b))
cat(sprintf("  AIC           %.2f          %.2f\n", AIC(model_full), AIC(model_best)))
cat(sprintf("  F-test p      %.2e       %.2e\n",
            pf(summary_full$fstatistic[1], summary_full$fstatistic[2], summary_full$fstatistic[3], lower=FALSE),
            fpval))

# ===========================================================================
# COMPUTE ALTERNATIVE METRICS FOR REPORT
# ===========================================================================
cat("\n========== ADDITIONAL ACCURACY METRICS ==========\n")

# RMSLE, Theil's U, etc.
n_len <- length(actual)
mse_b <- mean(resid_bb^2)

# Mean Absolute Percentage Deviation (symmetric MAPE — handles near-zero better)
smape_b <- mean(200 * abs(resid_bb) / (abs(actual) + abs(pred_b)))
cat(sprintf("SMAPE (symmetric):  %.2f%%\n", smape_b))

# R² decomposition
tss <- sum((actual - mean(actual))^2)
rss <- sum(resid_bb^2)
cat(sprintf("TSS=%.2f  RSS=%.2f  Explained=%.2f\n", tss, rss, tss-rss))

# Theil's U (U1 — bounded [0,1], lower is better)
theil_u1 <- sqrt(mse_b) / (sqrt(mean(actual^2)) + sqrt(mean(pred_b^2)))
cat(sprintf("Theil's U1:         %.4f  (0=perfect)\n", theil_u1))

# Mean error (bias)
me_b <- mean(resid_bb)
cat(sprintf("Mean Error (bias):  %.4f  (≈0 = unbiased)\n", me_b))

# ===========================================================================
# SAVE RESULTS
# ===========================================================================
cat("\n========== SAVING RESULTS ==========\n")
sink("analysis_results.txt")
cat("==============================================================\n")
cat("  MLR ANALYSIS — MOROCCO GDP GROWTH (V3: Agriculture model)\n")
cat("  Country: Morocco | Source: World Bank WDI | N=58 (1967-2024)\n")
cat("==============================================================\n\n")
cat("--- BEST MODEL (Forward Selection) ---\n"); cat(deparse(formula(model_best)), "\n\n")
cat("--- COEFFICIENTS ---\n"); print(round(summary_best$coefficients, 4))
cat(sprintf("\nR²=%.4f  Adj.R²=%.4f  RMSE=%.4f  AIC=%.1f\n",
            summary_best$r.squared, summary_best$adj.r.squared, rmse_b, AIC(model_best)))
cat(sprintf("DW=%.4f (p=%.4f)  BP=%.4f (p=%.4f)  JB=%.4f (p=%.4f)\n",
            dw$statistic, dw$p.value, bp$statistic, bp$p.value, jb$statistic, jb$p.value))
cat(sprintf("\nFORECAST 10%%: %.4f [%.4f, %.4f]\n",
            predict(model_best, newdata=sc10)[1],
            predict(model_best, newdata=sc10, interval="prediction", level=0.95)[2],
            predict(model_best, newdata=sc10, interval="prediction", level=0.95)[3]))
cat(sprintf("FORECAST 20%%: %.4f [%.4f, %.4f]\n",
            predict(model_best, newdata=sc20)[1],
            predict(model_best, newdata=sc20, interval="prediction", level=0.95)[2],
            predict(model_best, newdata=sc20, interval="prediction", level=0.95)[3]))
sink()

cat("\n==============================================================\n")
cat("  ALL 18 STEPS COMPLETE. Model performance: Adj.R²=%.4f.\n", summary_best$adj.r.squared)
cat("==============================================================\n")
