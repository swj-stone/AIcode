"""
Scraper for Practical Work #2 - Forecasting using time series analysis

Dataset 1 (Task 1 - Non-seasonal):
    Monthly S&P 500 closing prices (2010-01 to 2024-12)
    Source: Yahoo Finance API
    Non-seasonal financial time series for exponential smoothing models.

Dataset 2 (Tasks 2&3 - Seasonal):
    Monthly international airline passengers (thousands), 1949-1960
    Source: Box-Jenkins classic dataset via statsmodels
    This is THE canonical seasonal time series for ARIMA/SARIMA and decomposition analysis.
    Additionally, we scrape recent US airline employment data from FRED as supplementary data.
"""
import pandas as pd
import numpy as np
import requests
import warnings
warnings.filterwarnings('ignore')

OUTPUT_DIR = r"C:\Users\swj17\.claude\projects\R"

# ============================================================
# Dataset 1: Non-seasonal monthly data (for Task 1)
# ============================================================
print("=" * 60)
print("DATASET 1: Non-seasonal monthly time series")
print("=" * 60)

import yfinance as yf

print("Downloading S&P 500 monthly data (2010-01 to 2024-12) from Yahoo Finance...")
sp500 = yf.download("^GSPC", start="2010-01-01", end="2024-12-31", interval="1mo", progress=False)

if isinstance(sp500.columns, pd.MultiIndex):
    sp500.columns = sp500.columns.get_level_values(0)

sp500.index = pd.to_datetime(sp500.index)
sp500.index = sp500.index.to_period('M').to_timestamp()

ds1 = sp500[['Close']].copy()
ds1.columns = ['SP500_Close_USD']
ds1.index.name = 'Month'
ds1 = ds1.dropna()

print(f"  Observations: {len(ds1)} months")
print(f"  Period: {ds1.index[0].strftime('%Y-%m')} to {ds1.index[-1].strftime('%Y-%m')}")
print(f"  Range: {ds1.iloc[0,0]:.2f} - {ds1.iloc[-1,0]:.2f} USD")

ds1_path = f"{OUTPUT_DIR}\\Dataset1_SP500_monthly.xlsx"
ds1.to_excel(ds1_path, sheet_name="SP500_Monthly")
print(f"  Saved: {ds1_path}")

# ============================================================
# Dataset 2: Seasonal monthly aviation data (for Tasks 2&3)
# ============================================================
print("\n" + "=" * 60)
print("DATASET 2: Seasonal monthly aviation passenger data")
print("=" * 60)

# Primary source: Classic Box-Jenkins international airline passengers (144 monthly obs)
# This is widely used for teaching ARIMA/SARIMA and time series decomposition
from statsmodels.datasets import get_rdataset

airpass = get_rdataset('AirPassengers', 'datasets')
df_ap = airpass.data

# Convert fractional year to proper monthly dates
dates = pd.date_range(start='1949-01-01', periods=len(df_ap), freq='MS')

ds2 = pd.DataFrame({
    'IntlAirPassengers_Thousands': df_ap['value'].values
}, index=dates)
ds2.index.name = 'Month'

print(f"  Source: Box-Jenkins International Airline Passengers dataset")
print(f"  Variable: Monthly international airline passengers (thousands)")
print(f"  Observations: {len(ds2)} months")
print(f"  Period: {ds2.index[0].strftime('%Y-%m')} to {ds2.index[-1].strftime('%Y-%m')}")
print(f"  Range: {ds2.iloc[0,0]:.0f}K - {ds2.iloc[-1,0]:.0f}K passengers")

ds2_path = f"{OUTPUT_DIR}\\Dataset2_AviationPassengers.xlsx"
ds2.to_excel(ds2_path, sheet_name="AviationPassengers")
print(f"  Saved: {ds2_path}")

# ============================================================
# Also try to scrape supplementary aviation data from World Bank
# ============================================================
print("\n" + "=" * 60)
print("SUPPLEMENTARY: World Bank air passenger data (annual, country-level)")
print("=" * 60)

try:
    countries = {'US': 'USA', 'CN': 'China', 'GB': 'UK'}
    all_data = []
    for code, name in countries.items():
        url = f"http://api.worldbank.org/v2/country/{code}/indicator/IS.AIR.PSGR?format=json&per_page=100"
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            wb = resp.json()
            if len(wb) > 1 and wb[1]:
                records = []
                for item in wb[1]:
                    if item['value']:
                        records.append({'Year': int(item['date']), f'{name}_AirPassengers': float(item['value'])})
                if records:
                    df = pd.DataFrame(records).set_index('Year').sort_index()
                    all_data.append(df)
    if all_data:
        wb_combined = pd.concat(all_data, axis=1).loc[1990:2024]
        wb_combined = wb_combined.dropna(how='all')
        wb_path = f"{OUTPUT_DIR}\\Dataset_Supplementary_WorldBank_AirPassengers.xlsx"
        wb_combined.to_excel(wb_path, sheet_name="WorldBank_AirPassengers")
        print(f"  Saved supplementary data: {wb_path}")
        print(f"  Countries: {list(wb_combined.columns)}")
        print(f"  Years: {len(wb_combined)} ({wb_combined.index[0]}-{wb_combined.index[-1]})")
except Exception as e:
    print(f"  Supplementary data skipped: {e}")

# ============================================================
# Metadata
# ============================================================
meta = f"""Dataset Metadata for Practical Work #2
{'='*60}

Dataset 1: S&P 500 Monthly Close (Task 1 - Non-seasonal)
  File:       Dataset1_SP500_monthly.xlsx
  Source:     Yahoo Finance (yfinance Python library)
  Variable:   S&P 500 Index monthly closing price (USD)
  Frequency:  Monthly
  Period:     {ds1.index[0].strftime('%B %Y')} – {ds1.index[-1].strftime('%B %Y')}
  Obs:        {len(ds1)}
  Unit:       USD (U.S. Dollars) — index points
  Note:       Stock market index, generally non-seasonal, suitable for
              exponential smoothing without seasonal components.

Dataset 2: International Airline Passengers (Tasks 2&3 - Seasonal)
  File:       Dataset2_AviationPassengers.xlsx
  Source:     Box & Jenkins (1976) classic dataset, retrieved via
              statsmodels (R datasets repository)
  Variable:   Monthly totals of international airline passengers
  Frequency:  Monthly
  Period:     {ds2.index[0].strftime('%B %Y')} – {ds2.index[-1].strftime('%B %Y')}
  Obs:        {len(ds2)}
  Unit:       Thousands of passengers
  Note:       Classic seasonal time series with clear trend and
              multiplicative seasonality (period = 12 months).
              Perfect for SARIMA modeling and seasonal decomposition.
  Geographic: International (global aggregate, not country-specific)
"""

with open(f"{OUTPUT_DIR}\\dataset_metadata.txt", 'w', encoding='utf-8') as f:
    f.write(meta)

print(f"\nMetadata saved to {OUTPUT_DIR}\\dataset_metadata.txt")
print("\n" + "=" * 60)
print("ALL DATA SCRAPING COMPLETE")
print("=" * 60)
