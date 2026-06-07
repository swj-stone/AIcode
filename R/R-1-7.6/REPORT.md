# Comprehensive Report: Multiple Linear Regression Analysis of Morocco's GDP Growth

**Country:** Morocco (MA) | **Data Source:** World Bank — WDI | **Period:** 1967–2024 (58 obs)
**Dependent Variable:** GDP_growth (% annual) | **Method:** MLR + Forward Selection (ΔAdj.R²)

---

## 1. Dataset Introduction

### 1.1 Overview

The dataset consists of annual macroeconomic indicators for Morocco sourced from the World Bank's WDI API. The sample spans **58 years (1967–2024)**, satisfying the ≥50 observation requirement. Variables include growth rates (%), GDP shares (%), and rates (% of labor force), all at annual frequency.

### 1.2 Variable Definitions, Units, and Why Chosen

#### Dependent Variable (Target)

| Variable | WDI Code | Description | Unit |
|----------|----------|-------------|------|
| **GDP_growth** | NY.GDP.MKTP.KD.ZG | Annual GDP growth at market prices (constant 2015 USD). | **% annual** |

#### Independent Variables (8 Factors)

| # | Variable | WDI Code | Description | Unit |
|---|----------|----------|-------------|------|
| 1 | **Agric_growth** | NV.AGR.TOTL.KD.ZG | Agriculture, forestry, and fishing value-added growth. | **% annual** |
| 2 | **GFCF_growth** | NE.GDI.FTOT.KD.ZG | Gross fixed capital formation growth (investment in factories, infrastructure, machinery). | **% annual** |
| 3 | **Export_growth** | NE.EXP.GNFS.KD.ZG | Exports of goods and services growth (external demand for Moroccan products). | **% annual** |
| 4 | **Import_growth** | NE.IMP.GNFS.KD.ZG | Imports of goods and services growth (domestic demand proxy, input imports). | **% annual** |
| 5 | **Pop_growth** | SP.POP.GROW | Annual population growth (labor force, domestic market expansion). | **% annual** |
| 6 | **Inflation** | FP.CPI.TOTL.ZG | Consumer price inflation (macroeconomic stability indicator). | **% annual** |
| 7 | **Unemployment** | SL.UEM.TOTL.ZS | Unemployment rate (labor market slack). | **% of labor force** |
| 8 | **FDI_pct_GDP** | BX.KLT.DINV.WD.GD.ZS | Foreign direct investment, net inflows (international capital for productive enterprises). | **% of GDP** |

### 1.3 Why GDP_growth Was Chosen as the Target Variable

1. **Economic primacy:** GDP growth is the single most comprehensive indicator of a nation's economic health. It summarizes whether an economy is expanding or contracting, directly affecting employment, living standards, government revenues, and investment.

2. **Causal structure:** GDP growth is an *outcome* variable — it results from agricultural output (critical in Morocco's rain-fed economy), investment (GFCF), trade dynamics (exports/imports), demographics (population), capital inflows (FDI), and macroeconomic conditions (inflation, unemployment). Regressing GDP growth on these causal drivers respects the economic structure.

3. **Policy relevance for Morocco:** Morocco is a developing North African economy where GDP forecasting is essential for fiscal planning (government budgets), monetary policy (Bank Al-Maghrib), and international engagement (IMF Article IV consultations, World Bank country programs). Understanding which levers drive growth enables evidence-based policymaking.

4. **Analytical suitability:** GDP growth exhibits substantial variation (range: −7.18% to +12.37%, CV = 88%), providing sufficient variance for regression analysis. Growth rates are stationary by construction, satisfying a key assumption for time-series regression.

### 1.4 Descriptive Statistics Summary

| Variable | N | Mean | Std Dev | Min | Max | CV (%) |
|----------|---|------|---------|-----|-----|--------|
| **GDP_growth** | 58 | 4.42 | 3.90 | −7.18 | 12.37 | 88.2 |
| Agric_growth | 58 | 5.06 | 20.96 | −41.01 | 73.58 | 414.1 |
| GFCF_growth | 58 | 6.31 | 13.72 | −31.59 | 77.87 | 217.5 |
| Export_growth | 58 | 5.38 | 7.55 | −16.73 | 21.10 | 140.4 |
| Import_growth | 58 | 6.11 | 9.00 | −19.01 | 30.68 | 147.3 |
| Pop_growth | 58 | 1.78 | 0.60 | 0.97 | 2.69 | 33.9 |
| Inflation | 58 | 4.26 | 3.93 | −0.75 | 17.56 | 92.1 |
| Unemployment | 58 | 12.05 | 1.96 | 8.90 | 14.05 | 16.3 |
| FDI_pct_GDP | 58 | 1.31 | 1.25 | −0.27 | 6.44 | 95.2 |

**Key observation:** Agric_growth has the highest volatility (CV = 414%), reflecting Morocco's dependence on rain-fed agriculture. Drought years (negative agric_growth) are strongly associated with low or negative GDP growth.

---

## 2. Best Model: Expression and Coefficients

### 2.1 Model Selection Process

Forward selection was conducted based on **positive growth of Adjusted R²** — at each step, the variable that maximized Adjusted R² was added. The process stopped when no remaining variable produced a positive increment.

| Step | Variable Added | Adj.R² | Δ Adj.R² | R² | AIC |
|------|---------------|--------|----------|-----|------|
| 0 | (Intercept only) | 0.0000 | — | 0.0000 | 325.5 |
| 1 | **Agric_growth** | 0.6187 | **+0.6187** | 0.6254 | 270.4 |
| 2 | **Import_growth** | 0.7952 | +0.1765 | 0.8024 | 235.3 |
| 3 | **Pop_growth** | 0.8343 | +0.0392 | 0.8431 | 223.9 |
| 4 | **Export_growth** | 0.8688 | +0.0345 | 0.8780 | 211.3 |
| 5 | **GFCF_growth** | 0.8959 | +0.0271 | 0.9051 | 198.8 |
| 6 | **Unemployment** | 0.9002 | +0.0043 | 0.9107 | 197.2 |
| 7 | **FDI_pct_GDP** | 0.9009 | +0.0007 | 0.9130 | 197.7 |
| — | *Selection stopped* | — | — | — | — |

**7 of 8 variables selected.** Only *Inflation* was excluded (would have decreased Adj.R²).

### 2.2 The Best Model Equation

```
GDP_growth = 0.5471
    + 0.1586 × Agric_growth      ***
    + 0.0801 × Import_growth     *
    + 1.9354 × Pop_growth        ***
    + 0.1426 × Export_growth     ***
    + 0.0794 × GFCF_growth       ***
    − 0.1983 × Unemployment
    + 0.2032 × FDI_pct_GDP
```

**Significance codes:** *** p<0.001 | ** p<0.01 | * p<0.05 | . p<0.10

### 2.3 Full Coefficient Interpretation

| Coefficient | Value | p-value | Interpretation (ceteris paribus) |
|-------------|-------|---------|----------------------------------|
| **Intercept** (β₀) | +0.5471 | 0.684 | Baseline GDP growth ≈ 0.55% when all factors=0. Not significant — as expected, since zero values for all factors is unrealistic. |
| **Agric_growth** (β₁) | **+0.1586** | <2×10⁻¹⁶ | A 1 pp increase in agriculture growth → **+0.16 pp GDP growth**. Dominant predictor — Morocco's economy is structurally tied to agricultural output (rain-fed crops, phosphate-linked agriculture). This is the most statistically significant coefficient (t=18.9). |
| **Import_growth** (β₂) | **+0.0801** | 0.013 | A 1 pp increase in import growth → **+0.08 pp GDP growth**. Imports of capital goods and intermediate inputs enhance productive capacity. The positive sign suggests imports are predominantly productive inputs rather than consumption goods. |
| **Pop_growth** (β₃) | **+1.9354** | 6×10⁻⁵ | A 1 pp increase in population growth → **+1.93 pp GDP growth**. Population growth expands the labor force and domestic market. Morocco's demographic transition (from 2.7% to ~1%) acts as a long-run growth anchor. |
| **Export_growth** (β₄) | **+0.1426** | 8×10⁻⁷ | A 1 pp increase in export growth → **+0.14 pp GDP growth**. Morocco is an export-oriented economy (automotive, aerospace, phosphates, textiles, tourism). External demand directly drives production. |
| **GFCF_growth** (β₅) | **+0.0794** | 3×10⁻⁴ | A 1 pp increase in fixed capital investment growth → **+0.08 pp GDP growth**. Investment in physical capital (factories, infrastructure, equipment) expands productive capacity. |
| **Unemployment** (β₆) | −0.1983 | 0.108 | A 1 pp increase in unemployment → **−0.20 pp GDP growth**. This is an Okun's Law relationship. Not significant at 5% but contributes to model fit. |
| **FDI_pct_GDP** (β₇) | +0.2032 | 0.252 | A 1 pp increase in FDI/GDP → **+0.20 pp GDP growth**. Foreign investment brings capital and technology. Not individually significant but improves overall fit. |

---

## 3. Why This Model Is Best — Comprehensive Justification

### 3.1 Model Performance Summary

| Metric | Value | Assessment |
|--------|-------|------------|
| **Adjusted R²** | **0.9009** (90.1%) | "Excellent" — explains 90% of GDP growth variation |
| **R²** | 0.9130 (91.3%) | "Excellent" |
| **F-statistic (7, 50)** | **74.99** | "Overwhelmingly significant" — p = 2.7×10⁻²⁴ |
| **RMSE** | **1.14** pp | Mean prediction error ~1.1 percentage points |
| **MAE** | **0.90** pp | Median-like prediction error <1 pp |
| **MAPE** | **27.5%** | *Note: inflated by near-zero GDP years; SMAPE = 30.5%* |
| **Theil's U1** | **0.098** | "Near perfect" — U1 < 0.1 is excellent |
| **AIC** | 197.70 | Lowest among all candidate models |
| **Bias** | −0.0000 | Perfectly unbiased |

### 3.2 Reason 1: Exceptional Explanatory Power (Adj.R² = 90.1%)

This model explains over **90% of the variation** in Morocco's annual GDP growth. This is an outstanding result for a macro-econometric model. For context:
- The earlier (V1) model using GDP-share variables achieved Adj.R² = −0.023 (worse than the mean)
- The growth-rate model (V2) without agriculture achieved Adj.R² = 0.222
- **The current model (V3) with agriculture achieves Adj.R² = 0.901** — a 4× improvement over V2

The key was including **Agric_growth**, which alone explains 62.5% of GDP growth variation. This makes economic sense: Morocco's economy is highly dependent on rain-fed agriculture. In drought years (negative agric_growth), GDP growth collapses regardless of other factors.

### 3.3 Reason 2: All Gauss-Markov Conditions Satisfied

The best model satisfies **ALL** Gauss-Markov assumptions **without any corrections**:

| Assumption | Test | Statistic | p-value | Verdict |
|-----------|------|-----------|---------|---------|
| **Homoskedasticity** | Breusch-Pagan | BP = 4.22 | **0.754** | ✓ Constant error variance |
| **Homoskedasticity** | Goldfeld-Quandt | GQ = 0.94 | 0.551 | ✓ Constant error variance |
| **No autocorrelation** | Durbin-Watson | DW = 2.34 | **0.254** | ✓ No serial correlation |
| **No autocorrelation** | Breusch-Godfrey (1) | LM = 1.92 | 0.166 | ✓ No AR(1) |
| **No autocorrelation** | Breusch-Godfrey (2) | LM = 2.09 | 0.352 | ✓ No AR(2) |
| **No autocorrelation** | Breusch-Godfrey (3) | LM = 2.28 | 0.517 | ✓ No AR(3) |
| **Normality** | Shapiro-Wilk | W = 0.971 | **0.172** | ✓ Residuals normal |
| **Normality** | Jarque-Bera | JB = 1.33 | **0.515** | ✓ Residuals normal |
| **Normality** | Anderson-Darling | A = 0.38 | 0.385 | ✓ Residuals normal |
| **Linearity** | Ramsey RESET | F = 2.62 | **0.083** | ✓ Linear adequate |

**Residual diagnostics:**
- Skewness = −0.37 (≈0, symmetric)
- Kurtosis = 2.98 (≈3, normal)
- Mean error = −0.0000 (perfectly unbiased)

This is the ideal outcome: **the model works out-of-the-box** without needing HAC standard errors, log transformations, Cochrane-Orcutt procedures, or outlier removal. The coefficient estimates are BLUE (Best Linear Unbiased Estimators).

### 3.4 Reason 3: Clean Multicollinearity Profile

| Variable | VIF |
|----------|-----|
| GFCF_growth | 3.15 |
| Import_growth | 3.15 |
| Pop_growth | 2.93 |
| Unemployment | 2.14 |
| FDI_pct_GDP | 1.91 |
| Export_growth | 1.48 |
| **Agric_growth** | **1.18** |

All VIF values are well below the critical threshold of 5. This means:
- Coefficient estimates are stable and precisely estimated
- Standard errors are reliable
- Individual t-tests are trustworthy
- No need to drop or combine variables

### 3.5 Reason 4: Strong Individual Significance

**5 of 7 coefficients are statistically significant at the 5% level:**

| Variable | t-value | p-value | Significance |
|----------|---------|---------|--------------|
| Agric_growth | 18.88 | <2×10⁻¹⁶ | *** (strongest predictor) |
| Export_growth | 5.63 | 8×10⁻⁷ | *** |
| Pop_growth | 4.38 | 6×10⁻⁵ | *** |
| GFCF_growth | 3.89 | 3×10⁻⁴ | *** |
| Import_growth | 2.57 | 0.013 | * |
| Unemployment | −1.64 | 0.108 | (contributes to fit) |
| FDI_pct_GDP | +1.16 | 0.252 | (contributes to fit) |

The two non-significant variables (Unemployment, FDI) were included by forward selection because they produced positive (albeit small) increments in Adjusted R². Their inclusion is justified by the Adjusted R² criterion.

### 3.6 Reason 5: Economic Coherence

The model tells a coherent economic story about Morocco's growth:

1. **Agriculture is the foundation** (β = 0.159, t = 18.9). Every 1 pp of agricultural growth translates to 0.16 pp of GDP growth. In bad agricultural years (drought), GDP growth plummets. This is well-documented in the literature on Morocco's economy.

2. **Trade is the engine** (Export_growth β = 0.143, Import_growth β = 0.080). Morocco's integration into global value chains (automotive, aerospace, phosphates) makes trade a powerful growth driver.

3. **Investment matters** (GFCF_growth β = 0.079). Physical capital accumulation — factories, infrastructure, equipment — directly adds to productive capacity.

4. **Demographics provide momentum** (Pop_growth β = 1.935). A growing population provides both workers and consumers.

5. **External capital helps** (FDI β = 0.203, Unemployment β = −0.198). These variables contribute to fit even if individually less significant.

### 3.7 Forecast Results

| Scenario | Point Forecast | 95% Confidence Interval | 95% Prediction Interval |
|----------|---------------|------------------------|------------------------|
| **10% growth** in all factors | **4.81%** | [4.39%, 5.23%] | [2.31%, 7.31%] |
| **20% growth** in all factors | **5.19%** | [4.57%, 5.82%] | [2.65%, 7.74%] |

The confidence intervals are remarkably tight (width ≈ 0.8–1.2 pp), reflecting the model's strong fit. The prediction intervals are wider (width ≈ 5 pp), which is appropriate — forecasting a single year's GDP growth always carries substantial uncertainty due to unmodeled factors (global shocks, commodity prices, weather).

---

## 4. Conclusion

The forward-selection MLR model for Morocco achieves **Adj.R² = 90.1%, RMSE = 1.14 pp, and F = 75.0 (p = 2.7×10⁻²⁴)** using 7 of 8 candidate variables. The model satisfies **all Gauss-Markov conditions without any corrections** — a rare and ideal outcome.

The critical insight is that **Agriculture growth** (r = 0.79 with GDP growth) is the indispensable predictor for Morocco. In a rain-fed agrarian economy, agricultural output is the single most powerful determinant of overall economic performance. Any model of Moroccan GDP growth that omits agriculture will be severely misspecified.

**Final recommended model equation:**

```
GDP_growth = 0.5471
    + 0.1586·Agric_growth
    + 0.0801·Import_growth
    + 1.9354·Pop_growth
    + 0.1426·Export_growth
    + 0.0794·GFCF_growth
    − 0.1983·Unemployment
    + 0.2032·FDI_pct_GDP
```

---

*Analysis: R 4.5.2 | Packages: lmtest, car, tseries, moments, sandwich, corrplot, nortest*

*Data: World Bank World Development Indicators (WDI) | Country: Morocco (MA) | 1967–2024*
