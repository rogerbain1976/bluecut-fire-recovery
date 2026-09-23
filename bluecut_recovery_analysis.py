# Blue Cut Fire Recovery: Statistical Analysis of GEE Sample Points
# Author: Roger Bain
# Input: BlueCut_Recovery_Points.csv exported by BlueCut_Recovery_v1 (Earth Engine)
# Run in Google Colab: upload the CSV, then run this file top to bottom.

import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from scipy.spatial import cKDTree

# ---------------------------------------------------------------
# 1. Load and clean
# ---------------------------------------------------------------
d = pd.read_csv('BlueCut_Recovery_Points.csv').drop(columns=['system:index', '.geo'])
print('Points exported:', len(d))

# Remove developed land (NLCD 21-24): I-15, rail lines, structures
d = d[~d.NLCD2016.isin([21, 22, 23, 24])].copy()
print('Points after removing developed land:', len(d))

# -9999 = no clear pixel in that composite. Only Month6 and MTBS have any.
print('Missing values:\n', (d == -9999).sum()[lambda s: s > 0])

# ---------------------------------------------------------------
# 2. Recovery trajectory (summer composites)
# ---------------------------------------------------------------
periods = ['Prefire', 'Postfire', 'Year1', 'Year5', 'Year9', 'Year10']
print('\nMean NDVI by period:')
for p in periods:
    print(f'  {p:9s} {d[p + "_NDVI"].mean():.3f}')

print('\nRecovery relative to pre-fire NDVI:')
for y in ['Year1', 'Year5', 'Year9', 'Year10']:
    r = d[y + '_NDVI'] / d.Prefire_NDVI * 100
    print(f'  {y:7s} median {r.median():6.1f}%   at/above pre-fire {(r >= 100).mean()*100:5.1f}%'
          f'   within 90% {(r >= 90).mean()*100:5.1f}%')

# ---------------------------------------------------------------
# 3. Burn severity: dNBR (USGS classes) and check against MTBS
# ---------------------------------------------------------------
bins = [-9, 0.1, 0.27, 0.44, 0.66, 9]
labels = ['Unburned', 'Low', 'Mod-low', 'Mod-high', 'High']
d['sev'] = pd.cut(d.dNBR, bins, labels=labels)
print('\ndNBR severity (% of points):\n', (d.sev.value_counts(normalize=True).sort_index() * 100).round(1))

m = d[d.MTBS_Severity.between(1, 4)]
print('\nMy dNBR by MTBS severity class:\n', m.groupby('MTBS_Severity').dNBR.agg(['count', 'mean']).round(3))
print('Spearman correlation, dNBR vs MTBS:', round(m[['dNBR', 'MTBS_Severity']].corr('spearman').iloc[0, 1], 3))

# ---------------------------------------------------------------
# 4. What explains recovery? (change in NDVI, pre-fire to Year 10)
# Aspect is circular, so it enters as northness (cos) and eastness (sin).
# ---------------------------------------------------------------
d['north'] = np.cos(np.radians(d.Aspect_deg))
d['east'] = np.sin(np.radians(d.Aspect_deg))
d['chg'] = d.Year10_NDVI - d.Prefire_NDVI

# Local km coordinates and 1 km spatial blocks for clustered standard errors
d['x'] = (d.lon - d.lon.mean()) * 92.4
d['y'] = (d.lat - d.lat.mean()) * 111.0
d['block'] = np.floor(d.x).astype(int).astype(str) + '_' + np.floor(d.y).astype(int).astype(str)

print('\nVariance explained (R2), one predictor at a time:')
for f in ['Slope_deg', 'Elevation_m', 'north + east', 'dNBR', 'Prefire_NDVI']:
    print(f'  {f:14s} {smf.ols("chg ~ " + f, d).fit().rsquared:.4f}')

terrain = smf.ols('chg ~ Slope_deg + Elevation_m + north + east', d).fit()
print('\nTerrain only R2:', round(terrain.rsquared, 4))

full = smf.ols('chg ~ Slope_deg + Elevation_m + north + east + dNBR + Prefire_NDVI', d).fit(
    cov_type='cluster', cov_kwds={'groups': d.block})
print('Full model R2:', round(full.rsquared, 4))
print(full.summary().tables[1])

# ---------------------------------------------------------------
# 5. Spatial autocorrelation of residuals (Moran's I, 8 nearest neighbors)
# ---------------------------------------------------------------
xy = d[['x', 'y']].values
tree = cKDTree(xy)
_, idx = tree.query(xy, k=9)
idx = idx[:, 1:]
z = full.resid.values - full.resid.values.mean()
I = np.sum(z[:, None] * z[idx]) / (8 * np.sum(z ** 2))
print("\nMoran's I of residuals:", round(I, 3))

# ---------------------------------------------------------------
# 6. Class summaries (for StoryMap charts)
# ---------------------------------------------------------------
d['r9'] = d.Year9_NDVI / d.Prefire_NDVI * 100
d['r10'] = d.Year10_NDVI / d.Prefire_NDVI * 100
d['slope_class'] = pd.cut(d.Slope_deg, [0, 5, 15, 30, 45, 90],
                          labels=['Flat', 'Gentle', 'Moderate', 'Steep', 'Very Steep'], include_lowest=True)
d['aspect_class'] = pd.cut((d.Aspect_deg + 45) % 360, [0, 90, 180, 270, 360],
                           labels=['North', 'East', 'South', 'West'], include_lowest=True)
for g in ['slope_class', 'aspect_class', 'sev']:
    print('\n', d.groupby(g, observed=True).agg(n=('r10', 'size'),
          median_Year9=('r9', 'median'), median_Year10=('r10', 'median')).round(1))
