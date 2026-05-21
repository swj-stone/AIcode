"""
Run all 3 tasks in Python (statsmodels) to get numerical results for the report.
Uses efficient search grids to complete within a reasonable time.
"""
import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import json
import os
from scipy import stats
from statsmodels.tsa.stattools import adfuller, acf, pacf
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.tsa.holtwinters import SimpleExpSmoothing, Holt, ExponentialSmoothing
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.tsa.statespace.sarimax import SARIMAX
from statsmodels.tsa.seasonal import seasonal_decompose
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUTPUT_DIR = r"C:\Users\swj17\.claude\projects\R"

# Load datasets
ds1 = pd.read_excel(f"{OUTPUT_DIR}\\Dataset1_SP500_monthly.xlsx", index_col=0, parse_dates=True)
ds2 = pd.read_excel(f"{OUTPUT_DIR}\\Dataset2_AviationPassengers.xlsx", index_col=0, parse_dates=True)

sp500 = ds1.iloc[:, 0].dropna()
airpass = ds2.iloc[:, 0].dropna()

print("=" * 70)
print("RUNNING ALL 3 TASKS")
print("=" * 70)

results = {}

# =============================================================
# TASK 1: Exponential Smoothing
# =============================================================
print("\n--- TASK 1: Exponential Smoothing (S&P 500) ---")

train = sp500.iloc[:-12]
control = sp500.iloc[-12:]

# SES (1st order) - coarser grid for speed
alphas = np.arange(0.01, 1.0, 0.03)
best_mad, best_a_ses = float('inf'), 0.5
for a in alphas:
    try:
        m = SimpleExpSmoothing(train).fit(smoothing_level=a, optimized=False)
        fc = m.forecast(12)
        mad_v = np.mean(np.abs(control.values - fc.values))
        if mad_v < best_mad:
            best_mad, best_a_ses = mad_v, a
    except: pass
ses_mape = 100 * np.mean(np.abs((control.values - SimpleExpSmoothing(train).fit(smoothing_level=best_a_ses, optimized=False).forecast(12).values) / control.values))
print(f"SES: alpha={best_a_ses:.3f}, MAD={best_mad:.2f}, MAPE={ses_mape:.2f}%")

# Double ES (Brown)
best_mad2, best_a_br = float('inf'), 0.5
for a in alphas:
    try:
        m = Holt(train).fit(smoothing_level=a, smoothing_trend=a, optimized=False)
        fc = m.forecast(12)
        mad_v = np.mean(np.abs(control.values - fc.values))
        if mad_v < best_mad2:
            best_mad2, best_a_br = mad_v, a
    except: pass
brown_mape = 100 * np.mean(np.abs((control.values - Holt(train).fit(smoothing_level=best_a_br, smoothing_trend=best_a_br, optimized=False).forecast(12).values) / control.values))
print(f"Brown: alpha=beta={best_a_br:.3f}, MAD={best_mad2:.2f}, MAPE={brown_mape:.2f}%")

# Holt
alphas_h = np.arange(0.05, 1.0, 0.08)
betas_h = np.arange(0.05, 1.0, 0.08)
best_sse, best_ah, best_bh = float('inf'), 0.5, 0.1
for a in alphas_h:
    for b in betas_h:
        try:
            m = Holt(train).fit(smoothing_level=a, smoothing_trend=b, optimized=False)
            fc = m.forecast(12)
            sse = np.sum((control.values - fc.values)**2)
            if sse < best_sse:
                best_sse, best_ah, best_bh = sse, a, b
        except: pass
holt_mape = 100 * np.mean(np.abs((control.values - Holt(train).fit(smoothing_level=best_ah, smoothing_trend=best_bh, optimized=False).forecast(12).values) / control.values))
print(f"Holt: alpha={best_ah:.3f}, beta={best_bh:.3f}, SSE={best_sse:.2f}, MAPE={holt_mape:.2f}%")

# Best model
models_t1 = {
    'SES (1st order)': ses_mape,
    'Double ES (2nd order / Brown)': brown_mape,
    'Holt (1st order polynomial)': holt_mape,
}
best_t1 = min(models_t1, key=models_t1.get)
print(f">>> Best: {best_t1} (MAPE = {models_t1[best_t1]:.2f}%)")

# Fit best on full data & forecast 6
if 'SES' in best_t1 and 'Double' not in best_t1:
    final_t1 = SimpleExpSmoothing(sp500).fit(smoothing_level=best_a_ses)
elif 'Brown' in best_t1 or 'Double' in best_t1:
    final_t1 = Holt(sp500).fit(smoothing_level=best_a_br, smoothing_trend=best_a_br)
else:
    final_t1 = Holt(sp500).fit(smoothing_level=best_ah, smoothing_trend=best_bh)

fcst_t1 = final_t1.forecast(6)
print(f"6-month forecast: {[f'{v:.0f}' for v in fcst_t1.values]}")

results['task1'] = {
    'model_name': best_t1,
    'ses_alpha': round(best_a_ses, 3), 'ses_mape': round(ses_mape, 2),
    'brown_alpha': round(best_a_br, 3), 'brown_mape': round(brown_mape, 2),
    'holt_alpha': round(best_ah, 3), 'holt_beta': round(best_bh, 3), 'holt_mape': round(holt_mape, 2),
    'forecast_6m': [round(v, 1) for v in fcst_t1.values],
}

# =============================================================
# TASK 2: ARIMA / SARIMA
# =============================================================
print("\n--- TASK 2: ARIMA / SARIMA (Air Passengers) ---")

# Stationarity
adf_orig = adfuller(airpass)
adf_d1 = adfuller(airpass.diff().dropna())
print(f"ADF orig: p={adf_orig[1]:.4f}, diff1: p={adf_d1[1]:.4f}")
d_order = 1

acf_vals = acf(airpass, nlags=36)
pacf_vals = pacf(airpass, nlags=36)
print(f"ACF: lag1={acf_vals[1]:.3f}, lag12={acf_vals[12]:.3f}, lag24={acf_vals[24]:.3f}")

# Fit key SARIMA models
sarima_models = []
orders_to_try = [
    ((2,1,1),(1,1,1)), ((1,1,1),(1,1,1)), ((0,1,1),(1,1,1)),
    ((2,1,2),(1,1,1)), ((1,1,1),(2,1,1)), ((2,1,1),(0,1,1)),
    ((0,1,1),(0,1,1)), ((1,1,2),(1,1,1)),
]

best_manual_aic, best_manual_model, best_manual_order, best_manual_seas = float('inf'), None, None, None
for (p,d,q), (P,D,Q) in orders_to_try:
    try:
        m = SARIMAX(airpass, order=(p,d,q), seasonal_order=(P,D,Q,12),
                     enforce_stationarity=False, enforce_invertibility=False)
        fitted = m.fit(disp=False, maxiter=100, method='lbfgs')
        aic_val = fitted.aic
        res = fitted.resid.dropna()
        lb = acorr_ljungbox(res, lags=[12], return_df=True)
        lb_p = lb['lb_pvalue'].values[0]
        print(f"  SARIMA({p},{d},{q})({P},{D},{Q})[12] AIC={aic_val:.1f} LB_p={lb_p:.4f}")
        sarima_models.append({
            'order': f"({p},{d},{q})({P},{D},{Q})[12]", 'aic': round(aic_val, 2), 'lb_pvalue': round(lb_p, 4)
        })
        if aic_val < best_manual_aic:
            best_manual_aic = aic_val
            best_manual_model = fitted
            best_manual_order = (p,d,q)
            best_manual_seas = (P,D,Q,12)
    except Exception as e:
        print(f"  SARIMA({p},{d},{q})({P},{D},{Q})[12] FAILED: {str(e)[:50]}")

print(f"\nBest manual: SARIMA{best_manual_order} x {best_manual_seas}, AIC={best_manual_aic:.1f}")

# Exhaustive search with limited grid
best_auto_aic, best_auto_order, best_auto_seas = float('inf'), None, None
for p in [0,1,2,3]:
    for q in [0,1,2]:
        for P in [0,1,2]:
            for Q in [0,1,2]:
                try:
                    m = SARIMAX(airpass, order=(p,1,q), seasonal_order=(P,1,Q,12),
                                 enforce_stationarity=False, enforce_invertibility=False)
                    fitted = m.fit(disp=False, maxiter=80, method='lbfgs')
                    if fitted.aic < best_auto_aic:
                        best_auto_aic = fitted.aic
                        best_auto_order = (p,1,q)
                        best_auto_seas = (P,1,Q,12)
                except: pass
print(f"Best auto: SARIMA{best_auto_order} x {best_auto_seas}, AIC={best_auto_aic:.1f}")

# Residual diagnostics
resid = best_manual_model.resid.dropna()
sw_stat, sw_p = stats.shapiro(resid[:100])
lb_res = acorr_ljungbox(resid, lags=[12, 24], return_df=True)
print(f"Shapiro-Wilk: W={sw_stat:.3f}, p={sw_p:.4f}")
print(f"LB(12): X2={lb_res['lb_stat'].values[0]:.3f}, p={lb_res['lb_pvalue'].values[0]:.4f}")

# Forecast
fcst_obj = best_manual_model.get_forecast(18)
fcst_mean = fcst_obj.predicted_mean
fcst_ci = fcst_obj.conf_int(alpha=0.05)

# Determine final model
if best_auto_aic < best_manual_aic:
    final_order = best_auto_order
    final_seas = best_auto_seas
    final_aic = best_auto_aic
    sarima_source = "Exhaustive AIC search"
else:
    final_order = best_manual_order
    final_seas = best_manual_seas
    final_aic = best_manual_aic
    sarima_source = "Manual ACF/PACF analysis"

results['task2'] = {
    'd_order': d_order,
    'adf_orig_p': round(adf_orig[1], 4), 'adf_d1_p': round(adf_d1[1], 4),
    'acf_lag12': round(acf_vals[12], 3), 'acf_lag24': round(acf_vals[24], 3),
    'best_manual_order': str(best_manual_order), 'best_manual_seas': str(best_manual_seas),
    'best_manual_aic': round(best_manual_aic, 2),
    'best_auto_order': str(best_auto_order), 'best_auto_seas': str(best_auto_seas),
    'best_auto_aic': round(best_auto_aic, 2),
    'final_model': f"SARIMA{final_order} x {final_seas}",
    'final_aic': round(final_aic, 2),
    'sarima_source': sarima_source,
    'shapiro_w': round(sw_stat, 3), 'shapiro_p': round(sw_p, 4),
    'lb12_stat': round(lb_res['lb_stat'].values[0], 3), 'lb12_p': round(lb_res['lb_pvalue'].values[0], 4),
    'lb24_stat': round(lb_res['lb_stat'].values[1], 3), 'lb24_p': round(lb_res['lb_pvalue'].values[1], 4),
    'forecast_18m': [round(v, 1) for v in fcst_mean.values],
    'all_models': sarima_models,
}

# =============================================================
# TASK 3: Classical Decomposition
# =============================================================
print("\n--- TASK 3: Classical Decomposition (Air Passengers) ---")

decomp = seasonal_decompose(airpass, model='multiplicative', period=12)
trend_m = decomp.trend.dropna()
seasonal_m = decomp.seasonal.dropna()
resid_m = decomp.resid.dropna()

# Seasonal factors
monthly_factors = np.array([np.mean(airpass.values[i::12] / np.polyval(np.polyfit(np.arange(len(airpass)), airpass.values, 1), np.arange(len(airpass))[i::12])) for i in range(12)])

print(f"Trend: {trend_m.min():.0f} – {trend_m.max():.0f}")
print(f"Seasonal factor range: {seasonal_m.min():.3f} – {seasonal_m.max():.3f}")
print(f"Peak month: {np.argmax(monthly_factors)+1}, Trough: {np.argmin(monthly_factors)+1}")

# Error analysis - align all components to common index
common_idx = trend_m.index.intersection(seasonal_m.index).intersection(resid_m.index)
fitted_vals = (trend_m.loc[common_idx] * seasonal_m.loc[common_idx]).values
actual_vals = airpass.loc[common_idx].values
errors = resid_m.loc[common_idx].values
mae = np.mean(np.abs(actual_vals - fitted_vals))
rmse_val = np.sqrt(np.mean((actual_vals - fitted_vals)**2))
mape_val = 100 * np.mean(np.abs((actual_vals - fitted_vals) / actual_vals))

sw_e_stat, sw_e_p = stats.shapiro(errors[:100])
lb_e = acorr_ljungbox(errors, lags=[12], return_df=True)

print(f"MAE={mae:.2f}, RMSE={rmse_val:.2f}, MAPE={mape_val:.2f}%")
print(f"Shapiro-Wilk errors: W={sw_e_stat:.3f}, p={sw_e_p:.4f}")

# Classical decomposition forecast
t = np.arange(len(airpass))
trend_reg = np.polyfit(t, airpass.values, 1)
future_t = np.arange(len(airpass), len(airpass) + 24)
trend_fcst = np.polyval(trend_reg, future_t)

# Rotate seasonal factors to match forecast start
start_m = len(airpass) % 12
rotated_seas = np.roll(monthly_factors, -start_m)
classical_fcst = trend_fcst * np.tile(rotated_seas, 2)[:24]

# Also try Holt-Winters multiplicative for better forecast
hw_m = ExponentialSmoothing(airpass, seasonal_periods=12, trend='add', seasonal='mul',
                              initialization_method='estimated').fit()
hw_fcst = hw_m.forecast(24)

print(f"Classical fcst (first 6): {[f'{v:.0f}' for v in classical_fcst[:6]]}")
print(f"HW fcst (first 6): {[f'{v:.0f}' for v in hw_fcst.values[:6]]}")

results['task3'] = {
    'trend_range': f"{trend_m.min():.0f} – {trend_m.max():.0f}",
    'seasonal_min': round(seasonal_m.min(), 3), 'seasonal_max': round(seasonal_m.max(), 3),
    'peak_month': int(np.argmax(monthly_factors)) + 1,
    'trough_month': int(np.argmin(monthly_factors)) + 1,
    'seasonal_factors': [round(v, 3) for v in monthly_factors],
    'mae': round(mae, 2), 'rmse': round(rmse_val, 2), 'mape': round(mape_val, 2),
    'sw_error_w': round(sw_e_stat, 3), 'sw_error_p': round(sw_e_p, 4),
    'lb_error_stat': round(lb_e['lb_stat'].values[0], 3), 'lb_error_p': round(lb_e['lb_pvalue'].values[0], 4),
    'classical_forecast_24m': [round(v, 1) for v in classical_fcst],
    'hw_forecast_24m': [round(v, 1) for v in hw_fcst.values],
}

# ---- Generate ACF/PACF plots ----
os.makedirs(f"{OUTPUT_DIR}\\Task2_plots_py", exist_ok=True)
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
plot_acf(airpass, lags=36, ax=axes[0,0], title='ACF - Airline Passengers')
plot_pacf(airpass, lags=36, ax=axes[0,1], title='PACF - Airline Passengers')
diff1 = airpass.diff().dropna()
plot_acf(diff1, lags=36, ax=axes[1,0], title='ACF - 1st Difference')
plot_pacf(diff1, lags=36, ax=axes[1,1], title='PACF - 1st Difference')
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}\\Task2_plots_py\\acf_pacf.png", dpi=100)
plt.close()

# ---- Forecast plot for SARIMA ----
fig, ax = plt.subplots(figsize=(12, 6))
best_manual_fcst = best_manual_model.get_forecast(18)
ax.plot(np.arange(len(airpass)), airpass.values, 'b-', linewidth=2, label='Original')
ax.plot(np.arange(len(airpass)), best_manual_model.fittedvalues, 'g-', linewidth=1.5, label='Fitted')
fcst_idx = np.arange(len(airpass), len(airpass) + 18)
ax.plot(fcst_idx, fcst_mean.values, 'r-', linewidth=2, label='Forecast')
ci = best_manual_fcst.conf_int(alpha=0.05)
ax.fill_between(fcst_idx, ci.iloc[:,0], ci.iloc[:,1], alpha=0.3, color='gray', label='95% CI')
ax.set_title(f'SARIMA{best_manual_order} x {best_manual_seas} - 18-Month Forecast')
ax.legend()
ax.grid(alpha=0.3)
plt.savefig(f"{OUTPUT_DIR}\\Task2_plots_py\\sarima_forecast.png", dpi=100)
plt.close()

# ---- Decomposition plot ----
fig, axes = plt.subplots(4, 1, figsize=(12, 10))
axes[0].plot(airpass.index, airpass.values, 'b-', linewidth=1.5)
axes[0].set_title('Original Series')
axes[1].plot(trend_m.index, trend_m.values, 'g-', linewidth=1.5)
axes[1].set_title('Trend')
axes[2].plot(seasonal_m.index[:24], seasonal_m.values[:24], 'orange', linewidth=1.5)
axes[2].set_title('Seasonal (first 2 years)')
axes[3].plot(resid_m.index, resid_m.values, 'r-', linewidth=0.8)
axes[3].set_title('Random (Residual)')
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}\\Task2_plots_py\\decomposition.png", dpi=100)
plt.close()

# Save results
with open(f"{OUTPUT_DIR}\\analysis_results.json", 'w') as f:
    json.dump(results, f, indent=2)

print("\n" + "=" * 70)
print("ALL RESULTS SAVED to analysis_results.json")
print("=" * 70)

# Print key results for report
print("\n" + "=" * 70)
print("KEY RESULTS FOR REPORT")
print("=" * 70)
print(f"\nTask 1 Best Model: {best_t1}")
print(f"  SES alpha={best_a_ses:.3f}, MAPE={ses_mape:.2f}%")
print(f"  Brown alpha={best_a_br:.3f}, MAPE={brown_mape:.2f}%")
print(f"  Holt alpha={best_ah:.3f}, beta={best_bh:.3f}, MAPE={holt_mape:.2f}%")
print(f"\nTask 2 Best SARIMA: SARIMA{final_order} x {final_seas}")
print(f"  AIC={final_aic:.2f}")
print(f"  Manual AIC={best_manual_aic:.2f}, Auto AIC={best_auto_aic:.2f}")
print(f"\nTask 3 Decomposition:")
print(f"  MAPE={mape_val:.2f}%, MAE={mae:.2f}, RMSE={rmse_val:.2f}")
print(f"  Peak month: {np.argmax(monthly_factors)+1}, Trough: {np.argmin(monthly_factors)+1}")
