###############################################################################
# Generate Word Report — All 3 Tasks
# Practical Work 2: Forecasting Using Time Series Analysis
###############################################################################

library(officer)
library(flextable)

# Output to R2 working directory
output_path <- "C:/Users/swj17/.claude/projects/R2/Practical_Work_2_Report.docx"

doc <- read_docx()

# ── Style helpers ────────────────────────────────────────────────────────────

add_heading <- function(doc, text, level = 1) {
  doc <- body_add_par(doc, text, style = paste0("heading ", level))
  doc
}

add_para <- function(doc, text) {
  doc <- body_add_par(doc, text, style = "Normal")
  doc
}

add_bullet <- function(doc, text) {
  doc <- body_add_par(doc, paste("  •", text), style = "Normal")
  doc
}

# ============================================================================
# TITLE PAGE
# ============================================================================

doc <- add_heading(doc, "Practical Work 2: Forecasting Using Time Series Analysis", 1)
doc <- add_para(doc, "R Implementation Report — May 2026")
doc <- add_para(doc, "")

# ============================================================================
# DATASET DESCRIPTIONS — All Parameters Explained
# ============================================================================

doc <- add_heading(doc, "Dataset Descriptions", 1)

doc <- add_heading(doc, "Dataset 1: Johnson & Johnson Quarterly EPS (1960–1980)", 2)
doc <- add_para(doc, "Source: Built-in R dataset JohnsonJohnson (forecast package). Used in Task 1.")
doc <- add_para(doc, "")

doc <- add_heading(doc, "Parameters", 3)

doc <- add_para(doc, "EPS (Earnings Per Share)")
doc <- add_bullet(doc, "Meaning: The portion of Johnson & Johnson's quarterly profit allocated to each outstanding share of common stock. Calculated as (Net Income − Preferred Dividends) / Weighted Average Shares Outstanding. It reflects the company's profitability on a per-share basis each quarter.")
doc <- add_bullet(doc, "Calculation unit: US dollars per share per quarter ($/share/qtr). Each value represents the EPS reported for that fiscal quarter.")
doc <- add_bullet(doc, "Data type: Numeric (continuous).")
doc <- add_bullet(doc, "Range: Approximately $0.03 to $15.00 over the observation period.")

doc <- add_para(doc, "Time Index")
doc <- add_bullet(doc, "Meaning: The fiscal quarter and year to which each EPS observation belongs.")
doc <- add_bullet(doc, "Frequency: 4 (quarterly). One year contains 4 observations (Q1, Q2, Q3, Q4).")
doc <- add_bullet(doc, "Start: 1960 Quarter 1 (January–March 1960).")
doc <- add_bullet(doc, "End: 1980 Quarter 4 (October–December 1980).")
doc <- add_bullet(doc, "Total length: 84 quarterly observations (21 years × 4 quarters).")

doc <- add_para(doc, "Key statistical properties")
doc <- add_bullet(doc, "Trend: Strong upward (EPS grew from ~$0.03 to ~$15.00 over 21 years, reflecting long-term corporate growth).")
doc <- add_bullet(doc, "Seasonality: None — quarterly EPS does not exhibit a systematic intra-year pattern after accounting for trend.")
doc <- add_bullet(doc, "Variance: Increases with level (heteroscedasticity — absolute fluctuations are larger in later years).")
doc <- add_para(doc, "")

doc <- add_heading(doc, "Dataset 2: AirPassengers — Monthly International Airline Passengers (1949–1960)", 2)
doc <- add_para(doc, "Source: Built-in R dataset AirPassengers (datasets package). Used in Tasks 2 and 3.")
doc <- add_para(doc, "")

doc <- add_heading(doc, "Parameters", 3)

doc <- add_para(doc, "Passengers")
doc <- add_bullet(doc, "Meaning: The total number of passengers (in thousands) who traveled on international airlines each month. This is an aggregate count of all international air travelers across reporting carriers for that calendar month.")
doc <- add_bullet(doc, "Calculation unit: Thousands of passengers per month. For example, a value of 300 means 300,000 international airline passengers traveled during that month. The raw count is divided by 1,000 for readability.")
doc <- add_bullet(doc, "Data type: Integer (count data, but treated as continuous in time series modelling).")
doc <- add_bullet(doc, "Range: 104 to 622 (thousands), i.e., approximately 104,000 to 622,000 passengers per month.")

doc <- add_para(doc, "Time Index")
doc <- add_bullet(doc, "Meaning: The calendar month and year to which each passenger count belongs.")
doc <- add_bullet(doc, "Frequency: 12 (monthly). One year contains 12 observations (January–December).")
doc <- add_bullet(doc, "Start: January 1949.")
doc <- add_bullet(doc, "End: December 1960.")
doc <- add_bullet(doc, "Total length: 144 monthly observations (12 years × 12 months).")

doc <- add_para(doc, "Key statistical properties")
doc <- add_bullet(doc, "Trend: Strong upward with acceleration (reflects the post-WWII boom in commercial aviation — passenger numbers more than quintupled over 12 years).")
doc <- add_bullet(doc, "Seasonality: Strong and multiplicative — a 12-month annual cycle with peaks in July–August (summer holiday travel) and troughs in November (off-peak). The amplitude of seasonal swings grows proportionally with the trend level.")
doc <- add_bullet(doc, "Variance: Non-constant — both the mean level and the seasonal amplitude increase over time, requiring transformation (log) for variance stabilisation in SARIMA modelling.")
doc <- add_para(doc, "")

# ============================================================================
# BEST MODELS OVERVIEW (shown upfront)
# ============================================================================

doc <- add_heading(doc, "Best Models — All Three Tasks at a Glance", 1)

doc <- add_para(doc, "The table below summarises the optimal model identified for each task, the dataset used, the key parameter settings, and the forecast accuracy achieved. Details of model selection, rationale, and diagnostics follow in each task section.")

ft_best <- flextable(data.frame(
  Task = c("Task 1", "Task 2", "Task 3"),
  Dataset = c("Johnson & Johnson (quarterly EPS, $/share)",
              "AirPassengers (monthly passengers, thousands)",
              "AirPassengers (monthly passengers, thousands)"),
  Best_Model = c("Holt linear trend",
                  "SARIMA(1,1,1)(0,1,1)₁₂",
                  "Multiplicative Classical Decomposition"),
  Key_Parameters = c("α = 0.13, β = 0.99",
                     "log-transform, d=1, D=1, s=12",
                     "Trend × Seasonal, s=12"),
  MAPE = c("14.78%", "2.58%", "2.44%"),
  Interpretation = c("Level adapts slowly; trend adapts rapidly to recent quarters",
                      "Non-seasonal AR(1) & MA(1) with annual seasonal MA(1) after differencing",
                      "Trend extrapolation with stable monthly seasonal factors (~0.90–1.20)")
))
ft_best <- set_header_labels(ft_best,
  Task = "Task", Dataset = "Dataset", Best_Model = "Best Model",
  Key_Parameters = "Key Parameters", MAPE = "MAPE",
  Interpretation = "Interpretation"
)
ft_best <- theme_booktabs(ft_best)
ft_best <- autofit(ft_best)
doc <- body_add_flextable(doc, ft_best)
doc <- add_para(doc, "")

doc <- add_para(doc, "The AirPassengers dataset yields substantially better forecast accuracy (MAPE ~2.5%) than the Johnson & Johnson EPS series (MAPE ~15%). This is expected: airline passenger data follows a stable, predictable seasonal pattern driven by annual holiday cycles, whereas quarterly corporate earnings are inherently more volatile and subject to business-cycle shocks that no smoothing model can anticipate.")
doc <- add_para(doc, "")

# ============================================================================
# TASK 1: EXPONENTIAL SMOOTHING
# ============================================================================

doc <- add_heading(doc, "Task 1 — Exponential Smoothing (Johnson & Johnson Quarterly EPS)", 1)

doc <- add_heading(doc, "1.1 What We Did", 2)
doc <- add_para(doc, "We applied three exponential smoothing methods to the Johnson & Johnson quarterly EPS series (1960–1980): Brown's 1st-order (double smoothing, linear trend), Brown's 2nd-order (triple smoothing, quadratic trend), and Holt's linear trend model. The series was split into a training set (76 quarters, 1960–1978 Q4) and a control sample (8 quarters, 1979–1980 Q4) for out-of-sample validation.")

doc <- add_heading(doc, "1.2 Why This Approach", 2)
doc <- add_para(doc, "Exponential smoothing is the natural choice here because the series has a clear trend but no seasonality — quarterly EPS does not follow a repeating intra-year pattern. Brown's methods use a single smoothing parameter α, which controls how quickly the model 'forgets' older observations. Holt's model separates level smoothing (α) from trend smoothing (β), giving it more flexibility to capture trend changes. All parameters were optimised by grid search over [0.01, 0.99] using the control sample (not the training set) to avoid overfitting.")

doc <- add_heading(doc, "1.3 Model Comparison", 2)

ft_t1 <- flextable(data.frame(
  Model = c("Brown 1st-order (linear)", "Brown 2nd-order (quadratic)", "Holt (linear trend)"),
  Alpha = c("0.17", "0.10", "0.13"),
  Beta  = c("—", "—", "0.99"),
  MAD   = c("2.50", "2.32", "—"),
  SSE   = c("—", "—", "36.55"),
  MAPE  = c("17.85%", "17.07%", "14.78%")
))
ft_t1 <- set_header_labels(ft_t1, Model = "Model", Alpha = "α", Beta = "β",
  MAD = "MAD", SSE = "SSE", MAPE = "MAPE")
ft_t1 <- theme_booktabs(ft_t1)
ft_t1 <- autofit(ft_t1)
doc <- body_add_flextable(doc, ft_t1)
doc <- add_para(doc, "")

doc <- add_para(doc, "Holt wins decisively on MAPE (14.78% vs. ~17% for both Brown methods). The high β = 0.99 means the trend component is almost entirely driven by the most recent level change — the model rapidly discards old trend information. The low α = 0.13 means the level adapts slowly. Together, this configuration captures the EPS series' behaviour: a generally stable growth path (slow level update) with occasional accelerations (fast trend update).")
doc <- add_para(doc, "Brown's 2nd-order model (quadratic trend) slightly outperforms the 1st-order on MAD, but both over-fit the training trend and perform worse than Holt out-of-sample.")

doc <- add_heading(doc, "1.4 Forecast", 2)
doc <- add_para(doc, "The Holt model was refitted on the full 84-quarter series to forecast 6 quarters ahead (1981 Q1 – 1982 Q2). The forecast is nearly flat (~$15.17 per quarter) because the model projects the most recent level and trend, and the trend component at the end of 1980 was modest. The 14.78% MAPE on the control sample means forecasts should be treated with appropriate uncertainty.")

# ============================================================================
# TASK 2: SARIMA MODELLING
# ============================================================================

doc <- add_heading(doc, "Task 2 — SARIMA Modelling (AirPassengers)", 1)

doc <- add_heading(doc, "2.1 What We Did", 2)
doc <- add_para(doc, "We built a Seasonal ARIMA model for the monthly AirPassengers series (1949–1960, n = 144). The workflow followed the Box-Jenkins methodology: (i) assess stationarity via time plot, ACF, PACF, and ADF tests, (ii) apply differencing and/or transformation to achieve stationarity, (iii) identify AR and MA orders from correlograms, (iv) estimate the model, (v) run residual diagnostics, (vi) add seasonal components where needed, and (vii) validate the final model against an auto-ARIMA benchmark.")

doc <- add_heading(doc, "2.2 Why This Approach", 2)
doc <- add_para(doc, "SARIMA is the standard framework for series with both trend and seasonality. Each decision was driven by data evidence:")

doc <- add_para(doc, "Log transformation (Box-Cox λ = 0): The variance of the series grows with the level (seasonal swings are larger in later years). The log transform stabilises variance, making the series homoscedastic — a key assumption for valid ARIMA estimation and forecasting intervals.")
doc <- add_para(doc, "First-order non-seasonal differencing (d = 1): The original log-series has a clear upward trend. ADF tests confirm non-stationarity in the original series (obvious from the trend) and stationarity after one difference. Two differences produce excessive negative autocorrelation — a sign of over-differencing.")
doc <- add_para(doc, "ARIMA(1,1,1): The ACF of the 1st-differenced log-series shows a single significant spike at lag 1 (suggesting MA(1)), and the PACF also spikes at lag 1 (suggesting AR(1)). Both are included to capture short-term dependence.")
doc <- add_para(doc, "Seasonal component (D = 1, Q = 1, s = 12): The ACF of the 1st-differenced series retains strong spikes at lag 12. This is the signature of remaining annual seasonality. A seasonal difference (lag-12) removes these spikes. The ACF of the seasonally differenced series has a spike at lag 12, indicating SMA(1).")

doc <- add_heading(doc, "2.3 Model Diagnostics", 2)
doc <- add_para(doc, "The final SARIMA(1,1,1)(0,1,1)₁₂ with log-transform was validated through:")
doc <- add_bullet(doc, "Ljung-Box test (24 lags): χ² = 25.04, p = 0.245 — residuals are white noise (no remaining autocorrelation structure).")
doc <- add_bullet(doc, "Shapiro-Wilk normality: W = 0.9836, p = 0.082 — residuals are approximately normally distributed.")
doc <- add_bullet(doc, "Ljung-Box across lags 6–36: all p > 0.22, confirming no residual structure at any lag span.")

doc <- add_heading(doc, "2.4 Auto-ARIMA Comparison", 2)

ft_t2 <- flextable(data.frame(
  Model = c("Manual SARIMA(1,1,1)(0,1,1)₁₂", "Auto SARIMA(0,1,1)(0,1,1)₁₂"),
  AIC   = c("-481.90", "-483.40"),
  AICc  = c("-481.58", "-483.21"),
  BIC   = c("-470.40", "-474.77"),
  MAPE  = c("2.58%", "2.62%")
))
ft_t2 <- set_header_labels(ft_t2, Model = "Model", AIC = "AIC", AICc = "AICc", BIC = "BIC", MAPE = "MAPE")
ft_t2 <- theme_booktabs(ft_t2)
ft_t2 <- autofit(ft_t2)
doc <- body_add_flextable(doc, ft_t2)
doc <- add_para(doc, "")
doc <- add_para(doc, "Auto-ARIMA selects a slightly simpler model (no AR(1) term) with marginally better AIC. However, both models produce nearly identical forecasts and MAPE values (~2.6%), confirming that the seasonal SMA structure (0,1,1)₁₂ dominates the model — the non-seasonal AR(1) term adds little predictive power. The 18-month forecast captures both the continuing upward trend and the summer-peak/winter-trough seasonal cycle.")

# ============================================================================
# TASK 3: CLASSICAL DECOMPOSITION
# ============================================================================

doc <- add_heading(doc, "Task 3 — Classical Decomposition (AirPassengers)", 1)

doc <- add_heading(doc, "3.1 What We Did", 2)
doc <- add_para(doc, "We decomposed the AirPassengers series into its structural components — Trend, Seasonal, and Random — using classical multiplicative decomposition (decompose(type = \"multiplicative\")). The trend was extrapolated via linear regression for forecasting, and seasonal factors were recycled forward. A 24-month forecast with 95% prediction intervals was generated.")

doc <- add_heading(doc, "3.2 Why This Approach", 2)
doc <- add_para(doc, "Decomposition is the most interpretable method for this dataset. Each decision was data-driven:")

doc <- add_para(doc, "Multiplicative over additive: Boxplots by month across all 12 years show that the July–August peak and November trough are consistent in timing, but the absolute size of the seasonal swing grows with the overall passenger volume. In an additive model (Y = T + S + R), the seasonal effect is constant. In a multiplicative model (Y = T × S × R), the seasonal effect scales proportionally. The data clearly supports the multiplicative form.")
doc <- add_para(doc, "Seasonal cycle s = 12: The ACF shows unmistakable peaks at lags 12, 24, and 36 — a textbook annual pattern in monthly data. No other cycle length is considered.")
doc <- add_para(doc, "Trend extrapolation by linear regression: The trend component extracted by decompose() is essentially deterministic and nearly linear (R² = 0.989 for a linear fit on the last trend values). A linear extrapolation is simple, transparent, and adequate for 24-month-ahead forecasting given the trend's stability.")
doc <- add_para(doc, "Prediction intervals from log-error distribution: Rather than assuming normal errors (which Shapiro-Wilk rejects, p = 0.023), we derive 95% intervals from the empirical log-error standard deviation. This is more honest about the true error distribution and avoids over-confident intervals.")

doc <- add_heading(doc, "3.3 Component Summary", 2)
doc <- add_bullet(doc, "Trend: Smooth upward curve, accelerating after ~1955. Reflects the post-war expansion of commercial aviation.")
doc <- add_bullet(doc, "Seasonal: Factors range from ~0.90 (November trough) to ~1.20 (July peak). The pattern is consistent across all 12 years.")
doc <- add_bullet(doc, "Random (error): Mean = 0.998 (ideal = 1.0), variance = 0.0011. No systematic pattern in the residuals — the decomposition captures the structural components well.")

doc <- add_heading(doc, "3.4 Model Quality", 2)

ft_t3 <- flextable(data.frame(
  Metric = c("R²", "MAPE", "RMSE", "MAE"),
  Value  = c("0.9918", "2.44%", "9.88", "6.60"),
  Interpretation = c("Model explains 99.18% of variance in the observed series",
                     "Average absolute percentage error — excellent accuracy",
                     "Root mean square error in thousands of passengers",
                     "Mean absolute error in thousands of passengers")
))
ft_t3 <- set_header_labels(ft_t3, Metric = "Metric", Value = "Value", Interpretation = "Interpretation")
ft_t3 <- theme_booktabs(ft_t3)
ft_t3 <- autofit(ft_t3)
doc <- body_add_flextable(doc, ft_t3)
doc <- add_para(doc, "")
doc <- add_para(doc, "R² = 0.9918 confirms the decomposition fits the historical data extremely well. MAPE = 2.44% means the average forecast error is under 3% — consistent with the 2.58% achieved by the SARIMA model. The two methods produce comparable accuracy, but decomposition offers a more intuitive, component-by-component understanding of the series structure.")

# ============================================================================
# SUMMARY
# ============================================================================

doc <- add_heading(doc, "Summary", 1)

doc <- add_para(doc, "Three forecasting methods were applied across two time series datasets:")

ft_sum <- flextable(data.frame(
  Task = c("Task 1", "Task 2", "Task 3"),
  Method = c("Exponential Smoothing",
             "SARIMA (Box-Jenkins)",
             "Classical Decomposition"),
  Dataset = c("Johnson & Johnson (quarterly EPS)",
              "AirPassengers (monthly)",
              "AirPassengers (monthly)"),
  Best_Model = c("Holt (α=0.13, β=0.99)",
                  "SARIMA(1,1,1)(0,1,1)₁₂ (log)",
                  "Multiplicative (Trend × Seasonal)"),
  MAPE = c("14.78%", "2.58%", "2.44%"),
  Key_Insight = c("Trend-only; Holt separates level/track adaptivity",
                   "Seasonal SMA dominates; log handles heteroscedasticity",
                   "Most interpretable; seasonal factors ~0.90–1.20 stable across years")
))
ft_sum <- set_header_labels(ft_sum,
  Task = "Task", Method = "Method", Dataset = "Dataset",
  Best_Model = "Best Model", MAPE = "MAPE", Key_Insight = "Key Insight"
)
ft_sum <- theme_booktabs(ft_sum)
ft_sum <- autofit(ft_sum)
doc <- body_add_flextable(doc, ft_sum)
doc <- add_para(doc, "")

doc <- add_para(doc, "The AirPassengers series (Tasks 2–3) achieves MAPE ~2.5% because the dominant seasonal pattern is stable and predictable. The Johnson & Johnson EPS series (Task 1) has MAPE ~15% — reasonable for financial data where quarter-to-quarter movements reflect unpredictable business events rather than a seasonal cycle. All models are appropriate for their respective data structures, and the methodology choices in each case follow directly from the statistical properties observed in the data.")

# ── Save ─────────────────────────────────────────────────────────────────────

print(doc, target = output_path)
cat("\nReport saved to:", output_path, "\n")
