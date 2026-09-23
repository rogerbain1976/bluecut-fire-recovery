# Blue Cut Fire Vegetation Recovery, 2016 to 2026

Ten-year analysis of vegetation recovery after the August 2016 Blue Cut Fire in Cajon Pass, San Bernardino County, California, using Landsat 8/9 imagery in Google Earth Engine and statistical analysis in Python.

**Story:** [Ten Years After the Blue Cut Fire (ArcGIS StoryMap)](ADD-STORYMAP-LINK)

## Key findings

- The fire cut mean summer greenness (NDVI) roughly in half, from 0.319 to 0.159 across all 9,229 burned points.
- Relative to matched unburned land, the burned area recovered to 65% after one year, 89% after five, 93% after nine, and 98% after ten (95% CI: 94 to 102%). Ten years on, the burn is statistically indistinguishable from comparable unburned land.
- Unburned land also greened substantially in 2026. Without a control, the burned area would appear to have recovered to 117% of its pre-fire greenness.
- Burn severity (dNBR) was the strongest predictor of recovery. Slope, elevation, and aspect together explained less than 1% of point-level variation.
- dNBR agreed closely with MTBS burn severity (Spearman correlation 0.69).

## Repository contents

| File | Description |
|---|---|
| `bluecut_recovery_gee.js` | Earth Engine script, version 1. Scaled Landsat composites, NDVI/NBR/dNBR, 10,000 random points in the burn, exports sample CSV and GeoTIFFs. |
| `bluecut_recovery_gee_v2.js` | Earth Engine script, version 2. Adds an unburned control ring (500 m to 3 km outside the perimeter), MTBS fire history, and StoryMap asset exports. |
| `bluecut_recovery_analysis.py` | Python analysis of burned points: recovery trajectory, dNBR validation against MTBS, regression of recovery drivers with spatially clustered standard errors, Moran's I. |
| `bluecut_control_analysis.py` | Python analysis comparing burned points with matched unburned controls, with 1 km spatial block bootstrap confidence intervals. |
| `data/BlueCut_Recovery_Points.csv` | Sample data exported by version 1 (10,000 burned points). |
| `data/BlueCut_Recovery_Points_v2.csv` | Sample data exported by version 2 (10,000 burned and 8,000 control points). |

## Methods summary

1. **Imagery.** Landsat 8 and 9 Collection 2 Level-2 surface reflectance, scaled to true reflectance (DN x 0.0000275 - 0.2), with clouds, shadows, cirrus, snow, and saturated pixels masked.
2. **Composites.** Median summer composites for 2016 (pre-fire and post-fire), 2017, 2021, 2025, and 2026, plus a winter 2017 composite used for display only.
3. **Indices.** NDVI, NBR, and dNBR (pre-fire minus post-fire NBR), classified with USGS thresholds.
4. **Sampling.** 10,000 random points inside the MTBS perimeter (seed 2016) and 8,000 in the control ring (seed 2026). Points on developed land (NLCD 2016 classes 21 to 24) were removed, leaving 9,229 burned points.
5. **Controls.** Control points on developed land, water, barren land, or land MTBS mapped as burned since 1996 were excluded. Remaining controls were trimmed to the burned area's 5th to 95th percentile of elevation, slope, and pre-fire NDVI, then weighted to match across strata, leaving 987.
6. **Statistics.** Recovery relative to pre-fire NDVI; burned vs. control relative recovery; OLS models with standard errors clustered by 1 km spatial blocks; Moran's I on residuals; 1,000-rep spatial block bootstrap for confidence intervals.

## How to reproduce

**Earth Engine.** Register a noncommercial Earth Engine account, open the [Code Editor](https://code.earthengine.google.com), paste `bluecut_recovery_gee_v2.js`, and click Run. Run `BlueCut_Recovery_Points_v2` from the Tasks tab. Run `bluecut_recovery_gee.js` the same way for the version 1 export. Files are saved to a `BlueCut_GEE` folder in Google Drive.

**Python.** Requires pandas, numpy, scipy, and statsmodels (all preinstalled in Google Colab). From the repository root, run:

```
python bluecut_recovery_analysis.py
python bluecut_control_analysis.py
```

In Google Colab, clone the repository first:

```
!git clone https://github.com/ADD-GITHUB-USERNAME/bluecut-fire-recovery.git
%cd bluecut-fire-recovery
!python bluecut_recovery_analysis.py
!python bluecut_control_analysis.py
```

Results print to the console. Minor differences in the fourth decimal place can occur across pandas versions.

## Limitations

NDVI measures greenness, not species composition. Annual grasses can match the greenness of mature chaparral, so recovered greenness does not show that the original shrub community returned. Landsat pixels (30 m) blend mixed vegetation, and summer composites capture a single point in the growing season. NLCD 2016 classes inside the burn appear to reflect post-fire conditions and were used only to remove developed land.

## Background

This project rebuilds my 2025 graduate capstone (MS in GIS Technology, University of Arizona). The rebuild corrects the original workflow, which used unscaled reflectance values and fit statistical models to class averages rather than individual sample points.

## Data sources

- Landsat 8 and 9 Collection 2 Level-2 (USGS), via Google Earth Engine
- Monitoring Trends in Burn Severity perimeters and severity mosaics (USFS and USGS)
- USGS 3D Elevation Program, 10 m
- National Land Cover Database 2016 (MRLC)

## Author

Roger Bain, GIS analyst, Victorville, California. [LinkedIn](ADD-LINKEDIN-LINK)
