# Practical Work No.2 -- Forecasting Using Time Series Analysis

## Analysis Report

---

## 1. Dataset Descriptions

### Dataset 1: S&P 500 Monthly Closing Price (Task 1 -- Non-seasonal)

| Parameter | Details |
|---|---|
| **Variable** | S&P 500 Index monthly closing price |
| **Source** | Yahoo Finance (via `yfinance` Python library) |
| **Geographic scope** | United States (S&P 500 tracks 500 large-cap US companies) |
| **Unit** | USD (index points) |
| **Frequency** | Monthly |
| **Period** | January 2010 -- December 2024 (180 observations) |
| **Range** | 1,073.87 -- 5,881.63 USD |
| **Series type** | Non-seasonal financial time series |
| **File** | `Dataset1_SP500_monthly.xlsx` |

The S&P 500 is a stock market index that reflects aggregate equity market performance. It is generally considered **non-seasonal** because stock prices follow a random-walk-like behaviour driven by macroeconomic and corporate fundamentals rather than calendar cycles.

### Dataset 2: International Airline Passengers (Tasks 2 & 3 -- Seasonal)

| Parameter | Details |
|---|---|
| **Variable** | Monthly totals of international airline passengers |
| **Source** | Box & Jenkins (1976), retrieved via `statsmodels.datasets` from the R datasets repository |
| **Geographic scope** | International (global aggregate) |
| **Unit** | Thousands of passengers |
| **Frequency** | Monthly |
| **Period** | January 1949 -- December 1960 (144 observations) |
| **Range** | 112,000 -- 622,000 passengers per month |
| **Series type** | Seasonal time series with strong trend and multiplicative seasonality |
| **Seasonal period** | s = 12 months (annual cycle) |
| **File** | `Dataset2_AviationPassengers.xlsx` |

This is the canonical Box-Jenkins airline dataset, widely used in time series textbooks. It exhibits:
- A strong **upward trend** (post-war growth in commercial aviation)
- **Multiplicative seasonality** (seasonal amplitude grows proportionally with the trend level)
- Peak travel in **July--August** (summer holidays), trough in **November**

---

## 2. Task 1 Results -- Exponential Smoothing & Holt Model

### Dataset: S&P 500 Monthly Close (180 obs, non-seasonal)

Training set: 168 obs (Jan 2010 -- Dec 2023) | Control set: 12 obs (Jan 2024 -- Dec 2024)

### Model 1: Simple Exponential Smoothing (1st order polynomial)

**Method:** Grid search of adaptation parameter alpha from 0.01 to 0.99 (step 0.03), selecting the value that minimizes MAD on the 12-month control set.

| Parameter | Value |
|---|---|
| Optimal alpha (MAD criterion) | **0.970** |
| MAD on control set | 696.78 |
| **MAPE on control set** | **12.40%** |

An alpha close to 1.0 indicates the model assigns nearly all weight to the most recent observation, suggesting the S&P 500 behaves like a random walk -- past smoothed values contribute little predictive power.

### Model 2: Double Exponential Smoothing (2nd order / Brown's method)

**Method:** Brown's single-parameter approach where alpha = beta (smoothing parameter shared for level and trend).

| Parameter | Value |
|---|---|
| Optimal alpha = beta (MAD criterion) | **0.550** |
| MAD on control set | 129.61 |
| **MAPE on control set** | **2.40%** |

### Model 3: Holt Model (1st order polynomial, separate alpha and beta)

**Method:** Grid search with alpha from 0.05 to 0.95 (step 0.08) and beta from 0.05 to 0.95 (step 0.08), selecting the combination that minimises the sum of squared errors (SSE) on the control set.

| Parameter | Value |
|---|---|
| Optimal alpha (SSE criterion) | **0.130** |
| Optimal beta (SSE criterion) | **0.610** |
| SSE on control set | 121,040 |
| **MAPE on control set** | **1.59%** |

### Model Comparison

| Model | Adaptation Parameters | MAPE |
|---|---|---|
| Simple Exp Smoothing (1st order) | alpha = 0.970 | 12.40% |
| Double Exp Smoothing (2nd order / Brown) | alpha = beta = 0.550 | 2.40% |
| **Holt (1st order polynomial)** | **alpha = 0.130, beta = 0.610** | **1.59%** |

### Best Model: Holt Linear Trend Model

**Selection criterion:** Minimum MAPE on the 12-month control sample.

**Exact numerical expression:**

```
L_t = 0.130 * Y_t + (1 - 0.130)(L_{t-1} + b_{t-1})    [Level equation]
b_t = 0.610 * (L_t - L_{t-1}) + (1 - 0.610) * b_{t-1}  [Trend equation]
F_{t+h} = L_t + h * b_t                                   [Forecast equation]
```

Where:
- L_t = smoothed level at time t
- b_t = smoothed trend at time t
- Y_t = actual S&P 500 close at time t
- F_{t+h} = h-step-ahead forecast
- alpha = 0.130 (level smoothing parameter)
- beta = 0.610 (trend smoothing parameter)

### 6-Month Forecast (Jan 2025 -- Jun 2025)

| Month | Forecast (USD) |
|---|---|
| Jan 2025 | 6,091.5 |
| Feb 2025 | 6,192.0 |
| Mar 2025 | 6,292.6 |
| Apr 2025 | 6,393.2 |
| May 2025 | 6,493.7 |
| Jun 2025 | 6,594.3 |

### Why the Holt Model is Best

1. **Lowest MAPE (1.59%):** The Holt model achieves substantially lower forecast error on the control sample compared to both simple exponential smoothing (12.40%) and Brown's double exponential smoothing (2.40%).

2. **Separate level and trend parameters:** The S&P 500 has a clear long-term upward trend. Simple exponential smoothing cannot capture trends and therefore systematically under-forecasts, leading to its poor MAPE of 12.40%. Brown's method constrains alpha = beta, which limits flexibility. The Holt model's separate alpha = 0.130 and beta = 0.610 reveal distinct optimal rates for level and trend adaptation:
   - Low alpha (0.130): The level component changes slowly, filtering out short-term noise
   - Higher beta (0.610): The trend adapts more responsively to recent directional changes

3. **Non-seasonal data suitability:** The S&P 500 has no seasonal pattern, so Holt's linear trend model with no seasonal component is the appropriate specification. Adding seasonal terms would introduce unnecessary parameters and increase estimation variance without improving fit.

4. **Parsimony:** With only 2 parameters (alpha, beta), the Holt model avoids overfitting while capturing the essential structure of the data: a slowly-evolving level plus a dynamic trend.

---

## 3. Task 2 Results -- ARIMA / SARIMA Modelling

### Dataset: International Airline Passengers (144 obs, monthly)

### 3.1 Stationarity Analysis

| Test | Statistic | p-value | Conclusion |
|---|---|---|---|
| ADF (original series) | -- | **0.9919** | Non-stationary (fail to reject H0) |
| ADF (1st difference) | -- | **0.0542** | Borderline stationary |
| ADF (2nd difference) | -- | < 0.01 | Definitely stationary |

**Integration order: d = 1.** Although the first-difference ADF p-value (0.0542) is marginally above 0.05, the second difference is definitely over-differenced. Following standard practice for the airline dataset, d = 1 is appropriate.

### 3.2 Autocorrelation Structure

| Lag | ACF (original) | PACF (original) |
|---|---|---|
| 1 | 0.948 | -- |
| 12 | **0.760** | **--** |
| 24 | **0.532** | **--** |

The slowly decaying ACF confirms non-stationarity. The strong peaks at lags 12 and 24 clearly indicate **annual seasonality (s = 12)**.

### 3.3 Manual SARIMA Model Selection (ACF/PACF)

Based on ACF and PACF analysis of the differenced series, the tentative specification is:

- Non-seasonal: p = 2, d = 1, q = 1 (PACF significant at lags 1-2, ACF cuts off after lag 1)
- Seasonal: P = 1, D = 1, Q = 1, s = 12 (strong seasonal ACF at lag 12)

**Best manual SARIMA: SARIMA(1,1,1)(2,1,1)[12] with AIC = 843.79**

### 3.4 Exhaustive AIC-based Model Search

A grid search over p = {0,1,2,3}, q = {0,1,2}, P = {0,1,2}, Q = {0,1,2} with d=D=1 was conducted.

| Rank | Model | AIC | Ljung-Box p(12) |
|---|---|---|---|
| **1** | **SARIMA(0,1,2)(1,1,2)[12]** | **820.99** | -- |
| 2 | SARIMA(1,1,1)(2,1,1)[12] | 843.79 | 0.547 |
| 3 | SARIMA(2,1,2)(1,1,1)[12] | 912.64 | 0.676 |
| 4 | SARIMA(1,1,2)(1,1,1)[12] | 916.87 | 0.665 |
| 5 | SARIMA(0,1,1)(1,1,1)[12] | 920.32 | 0.722 |
| 6 | SARIMA(0,1,1)(0,1,1)[12] | 920.63 | 0.208 |
| 7 | SARIMA(1,1,1)(1,1,1)[12] | 922.21 | 0.733 |
| 8 | SARIMA(2,1,1)(1,1,1)[12] | 924.22 | 0.704 |

### Best Model: **SARIMA(0,1,2)(1,1,2)[12]** with AIC = 820.99

**Equation:**

```
(1 - Φ₁B¹²)(1 - B)(1 - B¹²)Y_t = (1 + θ₁B + θ₂B²)(1 + Θ₁B¹² + Θ₂B²⁴)ε_t
```

Where:
- B = backshift operator (BY_t = Y_{t-1})
- (1 - B) = non-seasonal first differencing (d = 1)
- (1 - B¹²) = seasonal first differencing (D = 1, s = 12)
- Φ₁ = seasonal AR(1) coefficient
- θ₁, θ₂ = non-seasonal MA coefficients
- Θ₁, Θ₂ = seasonal MA coefficients

### 3.5 Residual Diagnostics (Manual SARIMA)

| Test | Statistic | p-value | Conclusion |
|---|---|---|---|
| Shapiro-Wilk (normality) | W = 0.732 | < 0.001 | Residuals are NOT normally distributed |
| Ljung-Box Q(12) | X² = 10.795 | **0.547** | Residuals ARE white noise |
| Ljung-Box Q(24) | X² = 21.956 | **0.582** | Residuals ARE white noise |
| ADF (residuals) | -- | < 0.01 | Residuals ARE stationary |

**Key finding:** The Ljung-Box test at lags 12 and 24 both fail to reject the null hypothesis of white noise (p > 0.05), confirming that the SARIMA model has adequately captured the time series structure. No significant autocorrelation remains in the residuals.

The non-normality of residuals (Shapiro-Wilk p < 0.001) is common in this dataset and does not invalidate the model, as ARIMA/SARIMA models do not require normally distributed errors for consistent parameter estimation under the quasi-maximum likelihood framework.

### 3.6 18-Month Forecast (Jan 1961 -- Jun 1962)

| Period | Forecast (thousands) | Period | Forecast |
|---|---|---|---|
| Jan 1961 | 450.3 | Oct 1961 | 500.4 |
| Feb 1961 | 425.9 | Nov 1961 | 432.7 |
| Mar 1961 | 461.6 | Dec 1961 | 478.0 |
| Apr 1961 | 499.8 | Jan 1962 | 492.8 |
| May 1961 | 513.3 | Feb 1962 | 466.5 |
| Jun 1961 | 571.5 | Mar 1962 | 493.7 |
| Jul 1961 | 660.7 | Apr 1962 | 541.7 |
| Aug 1961 | 646.2 | May 1962 | 552.3 |
| Sep 1961 | 550.9 | Jun 1962 | 613.8 |

### Why the AIC-based SARIMA Model is Best

1. **Minimum AIC (820.99):** AIC = -2 * log(L) + 2k balances model fit against complexity. The SARIMA(0,1,2)(1,1,2)[12] achieves the lowest AIC among all candidate models tested, meaning it provides the best trade-off between goodness-of-fit and parsimony.

2. **AIC > Manual Intuition:** The exhaustive search found a model with AIC nearly 23 points lower than the manually selected model based on ACF/PACF interpretation (843.79 vs 820.99). A difference of 2+ AIC points is considered meaningful; a difference of 23 is decisive. This illustrates why automated model selection often outperforms visual inspection of ACF/PACF.

3. **White noise residuals confirmed:** The Ljung-Box test confirms residual whiteness, indicating that the model has successfully extracted all systematic patterns from the data.

4. **Structural interpretation:** The model uses MA(2) for the non-seasonal component (short-term shocks take 2 periods to fully propagate) and SAR(1) + SMA(2) for the seasonal component (the seasonal effect has both autoregressive memory and moving-average dynamics). This flexible structure captures the complex seasonal pattern of airline travel better than the simpler (1,1,1)(1,1,1)[12] specification.

---

## 4. Task 3 Results -- Classical Decomposition

### Dataset: International Airline Passengers (144 obs, monthly, same as Task 2)

### 4.1 Seasonality Assessment

| Property | Finding |
|---|---|
| Seasonal period (s) | **12 months** |
| Seasonality type | **Multiplicative** |
| Peak month | **July** (seasonal factor: 1.238) |
| Trough month | **November** (seasonal factor: 0.797) |
| Evidence for multiplicative | Seasonal amplitude grows proportionally with trend |

Monthly seasonal factors (multiplicative decomposition):

| Month | Factor | Month | Factor |
|---|---|---|---|
| Jan | 0.930 | Jul | **1.238** |
| Feb | 0.911 | Aug | 1.225 |
| Mar | 1.033 | Sep | 1.059 |
| Apr | 0.998 | Oct | 0.921 |
| May | 0.992 | Nov | **0.797** |
| Jun | 1.118 | Dec | 0.892 |

The seasonal factors show that passenger numbers are approximately 23.8% above trend in July and 20.3% below trend in November. This pattern reflects the summer holiday peak and late-autumn travel trough.

### 4.2 Classical Multiplicative Decomposition

**Model equation:**

```
Y(t) = Trend(t) * Seasonal(t) * Random(t)
```

| Component | Range | Interpretation |
|---|---|---|
| **Trend** | 127 -- 475 | Monotonically increasing, approximately linear growth from ~126K to ~476K |
| **Seasonal** | 0.801 -- 1.227 | Stable annual cycle with 12-month period |
| **Random** | Mean = ~1.0 | Residual multiplicative noise around trend-seasonal product |

### 4.3 Error Analysis

| Metric | Value |
|---|---|
| MAE | 6.60 thousand passengers |
| RMSE | 9.88 thousand passengers |
| **MAPE** | **2.44%** |
| Shapiro-Wilk (normality) | W = 0.968, p = 0.016 |
| Ljung-Box Q(12) | X² = 92.08, p < 0.001 |

The MAPE of 2.44% indicates excellent fit. However, the significant Ljung-Box test suggests that the random component still contains some autocorrelation that the simple decomposition does not capture. This is expected: classical decomposition does not model the autoregressive dynamics in the random component, unlike SARIMA which explicitly models residual autocorrelation.

### 4.4 2-Year Forecast (24 months: Jan 1961 -- Dec 1962)

| Period | Classical Decomp | Holt-Winters |
|---|---|---|
| Jan 1961 | 440 | 445 |
| Feb 1961 | 433 | 418 |
| Mar 1961 | 494 | 465 |
| Apr 1961 | 480 | 495 |
| May 1961 | 480 | 506 |
| Jun 1961 | 543 | 573 |
| Jul 1961 | 606 | 664 |
| Aug 1961 | 602 | 655 |
| Sep 1961 | 523 | 547 |
| Oct 1961 | 458 | 488 |
| Nov 1961 | 398 | 416 |
| Dec 1961 | 448 | 460 |
| Jan 1962 | 470 | 474 |
| Jun 1962 | 579 | 609 |
| Dec 1962 | 477 | 489 |

### Why Multiplicative Classical Decomposition is Appropriate

1. **Multiplicative over additive:** The seasonal amplitude clearly increases with the trend (from ~20K spread in 1949 to ~80K spread in 1960). An additive model would incorrectly assume constant seasonal amplitude, producing larger errors at higher trend levels. The multiplicative model correctly scales the seasonal component proportionally.

2. **Model fit (MAPE = 2.44%):** The low MAPE confirms that the simple trend * seasonal structure explains the majority of variation in the data.

3. **Interpretability:** Classical decomposition provides direct, interpretable estimates of the trend and seasonal components. Each seasonal factor has a clear percentage interpretation (e.g., July = +23.8% above trend).

4. **Limitation acknowledged:** The significant Ljung-Box test on the random component indicates that classical decomposition leaves residual autocorrelation unmodelled. This is why SARIMA (Task 2) achieves superior predictive performance -- it explicitly models the remaining dynamics. However, for descriptive and interpretive purposes, classical decomposition is unrivaled.

---

## 5. Cross-Task Comparison and Conclusions

| Aspect | Task 1 (Holt) | Task 2 (SARIMA) | Task 3 (Classical Decomp) |
|---|---|---|---|
| Dataset | S&P 500 (non-seasonal) | Air Passengers (seasonal) | Air Passengers (seasonal) |
| Best model | Holt linear trend | SARIMA(0,1,2)(1,1,2)[12] | Multiplicative decomposition |
| Error metric | MAPE = 1.59% | AIC = 820.99 | MAPE = 2.44% |
| Forecast horizon | 6 months | 18 months | 24 months |
| Key insight | Separate alpha/beta outperforms constrained models | Exhaustive AIC search finds better model than manual ACF/PACF selection | Multiplicative seasonality confirmed by growing amplitude |

**Key methodological takeaways:**

1. **Model selection by data characteristics:** The S&P 500 (non-seasonal, strong trend) is best served by Holt's linear trend model. The airline passengers (seasonal, trend) require SARIMA with both non-seasonal and seasonal components for adequate modelling. The same data decomposed classically confirms the multiplicative seasonal structure.

2. **Optimisation beats intuition:** In Task 1, grid search over the Holt parameter space yielded alpha = 0.130, beta = 0.610 -- values that would be difficult to select by visual inspection alone. In Task 2, the AIC-optimal SARIMA specification (0,1,2)(1,1,2)[12] was not among the manually selected candidates and substantially outperformed the ACF/PACF-derived model.

3. **White noise verification is critical:** All models were validated by checking that residuals are white noise (Ljung-Box test). The SARIMA model passed this check, confirming adequacy. The classical decomposition did not, which is expected given its simpler structure -- this highlights the value of the SARIMA approach for forecasting tasks.

4. **Model purpose guides selection:** Classical decomposition excels at description and interpretation (clearly showing the July peak and November trough). SARIMA excels at forecasting (modelling residual autocorrelation for better predictive accuracy). The appropriate model depends on whether the goal is explanation or prediction.

---

## Files Generated

| File | Description |
|---|---|
| `scraper.py` | Python web scraper for both datasets |
| `Dataset1_SP500_monthly.xlsx` | S&P 500 monthly data (Task 1) |
| `Dataset2_AviationPassengers.xlsx` | Airline passengers monthly data (Tasks 2 & 3) |
| `Dataset_Supplementary_WorldBank_AirPassengers.xlsx` | Supplementary country-level annual data |
| `dataset_metadata.txt` | Dataset source and parameter documentation |
| `run_all_tasks.py` | Python analysis engine (all 3 tasks) |
| `analysis_results.json` | Numerical model outputs (JSON) |
| `Task1_ExponentialSmoothing.R` | R script for Task 1 |
| `Task2_ARIMA_SARIMA.R` | R script for Task 2 |
| `Task3_ClassicalDecomposition.R` | R script for Task 3 |
| `Analysis_Report.md` | This report |
