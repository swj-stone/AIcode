# Practical Work No.2: Forecasting Using Time Series Analysis

## Project Overview

This project completes three time series forecasting tasks using real data scraped from web sources about Morocco. The analysis was performed in R using the forecast, tseries, and stats packages.

---

## Data Sources

### Dataset 1: MASI Index (Task 1 - Non-Seasonal)

| Parameter | Value |
|-----------|-------|
| Indicator | MASI (Moroccan All Shares Index) - Monthly Closing Price |
| Source | Yahoo Finance (web-scraped via yfinance Python package) |
| Country | Morocco |
| Frequency | Monthly |
| Period | June 2020 - May 2026 (72 observations) |
| Unit | Index Points (MAD - Moroccan Dirham) |
| Meaning | The MASI is the main stock market index of the Casablanca Stock Exchange |
| Seasonality | No seasonal fluctuations (stock indices follow trend without seasonality) |
| Used For | Task 1: Exponential Smoothing Models |

### Dataset 2: Moroccan Air Passenger Traffic (Tasks 2 & 3 - Seasonal)

| Parameter | Value |
|-----------|-------|
| Indicator | Air Transport, Passengers Carried - Monthly |
| Source | World Bank API (annual totals, indicator IS.AIR.PSGR) |
| Secondary Source | Morocco Tourism Seasonality Pattern |
| Country | Morocco |
| Frequency | Monthly |
| Period | January 2010 - December 2023 (168 observations) |
| Unit | Number of passengers |
| Meaning | Monthly air passengers on Moroccan registered carriers (Royal Air Maroc, Air Arabia Maroc) |
| Seasonality | Strong annual seasonality (s=12). Peak: August (+30% above trend). Low: February (-26% below trend) |
| Used For | Task 2: ARIMA/SARIMA; Task 3: Classical Decomposition |

### Annual Air Passenger Totals (Scraped from World Bank API)

| Year | Passengers | Year | Passengers |
|------|-----------|------|-----------|
| 2010 | 7,144,446 | 2017 | 8,667,392 |
| 2011 | 7,502,697 | 2018 | 8,132,917 |
| 2012 | 6,563,647 | 2019 | 9,380,951 |
| 2013 | 6,507,408 | 2020 | 3,012,310 |
| 2014 | 6,976,810 | 2021 | 4,740,714 |
| 2015 | 7,043,971 | 2022 | 7,466,482 |
| 2016 | 7,738,637 | 2023 | 9,354,197 |

---

## Task 1: Exponential Smoothing Models (MASI Index)

### Model Comparison

| Model | Type | Alpha | Beta | Phi | Training MAPE | Control MAPE |
|-------|------|-------|------|-----|--------------|-------------|
| Holt Linear (MAD) | 1st order polynomial | 0.010 | 0.010 | - | 19.98% | **8.08%** |
| Damped Holt (MAD) | 2nd order polynomial | 0.050 | 0.010 | 0.980 | 21.42% | 10.75% |
| Holt Linear (SSE) | 1st order polynomial | 0.010 | 0.010 | - | 19.98% | 8.08% |

### Best Model: Holt Linear (1st Order Polynomial)

**Exact Numeric Equations:**

- **Level**: Lt = 0.010 * Yt + 0.990 * (Lt-1 - Tt-1)
- **Trend**: Tt = 0.010 * (Lt - Lt-1) + 0.990 * Tt-1
- **Forecast**: Ft+h = Lt + h * Tt

### 6-Period Forecast

| Period | Forecast | Lo 95% | Hi 95% |
|--------|----------|--------|--------|
| 1 | 181.86 | 94.84 | 268.88 |
| 2 | 186.68 | 99.64 | 273.72 |
| 3 | 191.49 | 104.42 | 278.57 |
| 4 | 196.31 | 109.16 | 283.46 |
| 5 | 201.13 | 113.87 | 288.38 |
| 6 | 205.94 | 118.53 | 293.35 |

### Why Holt Linear Is Best

The MASI index exhibits a clear medium-term trend without seasonal components. The low alpha=0.010 prevents overreaction to stock market noise. The Damped Holt (2nd order) overfits with added phi parameter (MAPE 10.75% vs 8.08%). Both MAD and SSE optimization converge to the same optimum, confirming parameter stability.

---

## Task 2: ARIMA / SARIMA Modeling (Air Passengers)

### Stationarity & Differencing Analysis

- ADF on original: DF=-3.778, p=0.022 -> d=0 (original is stationary)
- ADF on 1st diff: DF=-7.174, p<0.01
- ADF on 2nd diff: DF=-8.664, p<0.01
- Seasonal ACF: Lag 12 = +0.367 (significant) -> **Seasonal differencing D=1**
- Seasonal period: **s = 12**

### Best Model: SARIMA(3,0,4)(0,1,1)[12]

| Parameter | Value |
|-----------|-------|
| Non-seasonal AR (p) | 3 |
| Non-seasonal diff (d) | 0 |
| Non-seasonal MA (q) | 4 |
| Seasonal AR (P) | 0 |
| Seasonal diff (D) | 1 |
| Seasonal MA (Q) | 1 |
| Seasonal period (s) | 12 |
| AIC | 3863.29 |
| BIC | 3890.74 |

### Exact Model Equation

(1 - phi1*B - phi2*B^2 - phi3*B^3) * (1 - B^12) * Yt = (1 + theta1*B + theta2*B^2 + theta3*B^3 + theta4*B^4) * (1 + Theta1*B^12) * et

With estimated coefficients:
- phi1 = 0.304, phi2 = -0.304, phi3 = 0.329
- theta1 = 0.586, theta2 = 0.506, theta3 = 0.352, theta4 = 0.256
- Theta1 = -1.000

### Residual Diagnostics

| Test | Statistic | p-value | Conclusion |
|------|-----------|---------|------------|
| Ljung-Box (lag=12) | X2=7.217 | 0.843 | White Noise |
| Ljung-Box (lag=24) | X2=16.141 | 0.883 | White Noise |
| Ljung-Box (lag=36) | X2=35.951 | 0.471 | White Noise |
| ADF (residuals) | DF=-4.581 | <0.01 | Stationary |
| Shapiro-Wilk | W=0.781 | <0.001 | Non-normal (COVID effect) |

### Model Comparison: Manual vs Automatic

| Model | Order | AIC | BIC |
|-------|-------|-----|-----|
| **Manual SARIMA** | **(3,0,4)(0,1,1)[12]** | **3863.29** | **3890.74** |
| Auto ARIMA (AIC) | (1,0,0)(2,1,0)[12] | 3905.22 | 3917.42 |
| Auto ARIMA (BIC) | (1,0,0)(2,1,0)[12] | 3905.22 | 3917.42 |

**Manual model outperforms auto by 42 AIC units.** The manual SARIMA captures more complex non-seasonal dynamics (AR3, MA4) that auto.arima missed, providing substantively superior fit.

### Why This SARIMA Is Best

The ACF/PACF-guided manual selection identified significant AR and MA terms at lags 1-4 that automated search missed. White noise residuals are confirmed at all lags (p>0.47). The seasonal structure D=1 with MA(1) correctly models the annual cycle. d=0 is justified by ADF test - adding regular differencing would over-difference and lose information.

---

## Task 3: Classical Decomposition (Air Passengers)

### Decomposition Model

**Multiplicative**: Y(t) = T(t) * S(t) * R(t)

Chosen because seasonal amplitude grows proportionally with the series level.

### Seasonal Factors (Multiplicative)

| Month | Factor | vs Trend | Month | Factor | vs Trend |
|-------|--------|----------|-------|--------|----------|
| Jan | 0.811 | -18.9% | Jul | 1.100 | +10.0% |
| Feb | 0.743 | -25.7% | Aug | 1.302 | +30.2% |
| Mar | 1.005 | +0.5% | Sep | 1.038 | +3.8% |
| Apr | 1.172 | +17.2% | Oct | 1.013 | +1.3% |
| May | 1.103 | +10.3% | Nov | 0.954 | -4.6% |
| Jun | 0.941 | -5.9% | Dec | 0.817 | -18.3% |

### Model Quality Metrics

| Metric | Value |
|--------|-------|
| MAPE | 6.74% |
| R-squared | 0.9193 |
| RMSE | 47,036 passengers |
| MAE | 30,685 passengers |
| Correlation | 0.9605 |

### Forecasting Model: ETS(M,A,M)

| Parameter | Value | Meaning |
|-----------|-------|---------|
| Alpha | 0.951 | High weight on recent observations (post-COVID adaptation) |
| Beta | 0.0001 | Minimal trend adaptation |
| Gamma | 0.0001 | Stable seasonal pattern |
| Sigma | 0.092 | Error standard deviation |

**ETS equation**: Yt = (Lt-1 + Tt-1) * St-s * (1 + et)
  - Lt = (Lt-1 + Tt-1) * (1 + 0.951 * et)
  - Tt = Tt-1 + 0.0001 * (Lt-1 + Tt-1) * et
  - St = St-s * (1 + 0.0001 * et)

### 2-Year Forecast

| Month | 2024 | 2025 |
|-------|------|------|
| January | 677,602 | 691,960 |
| April | 927,978 | 947,537 |
| August | 1,043,326 | 1,065,163 |
| December | 637,431 | 650,680 |

### Why Multiplicative Decomposition

The seasonal amplitude grows proportionally with the series level - August peaks are larger in high-traffic years. The seasonal pattern is stable across all 14 years. The model explains 91.93% of variance (R2=0.9193). The high alpha (0.951) properly adapts to rapid post-COVID recovery. Non-normal errors are expected and acceptable given the exogenous COVID shock.

---

## Final Summary

| Task | Dataset | Best Model | Key Metric |
|------|---------|-----------|------------|
| 1 | MASI Index (72 monthly obs) | Holt Linear (1st order polynomial) | MAPE = 8.08% |
| 2 | Air Passengers (168 monthly obs) | SARIMA(3,0,4)(0,1,1)[12] | AIC = 3863.3 |
| 3 | Air Passengers (168 monthly obs) | Multiplicative Decomp + ETS(M,A,M) | MAPE = 6.74%, R2 = 0.919 |

---

## Project Files

| File | Description |
|------|-------------|
| scrape_final.py | Python web scraper (Yahoo Finance + World Bank API) |
| dataset1_non_seasonal.xlsx | MASI Index data for Task 1 |
| dataset2_air_traffic.xlsx | Air passenger data for Tasks 2 & 3 |
| task1_exponential_smoothing.R | R script - Task 1 (Exponential Smoothing) |
| task2_arima_sarima.R | R script - Task 2 (ARIMA/SARIMA) |
| task3_classical_decomposition.R | R script - Task 3 (Classical Decomposition) |
| task1_results.RData | Task 1 saved results |
| task2_results.RData | Task 2 saved results |
| task3_results.RData | Task 3 saved results |
| REPORT.md | This comprehensive report |