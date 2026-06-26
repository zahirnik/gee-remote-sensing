# 03 - Drought, Soil Moisture and Groundwater

This folder collects Google Earth Engine scripts I use for operational
drought, soil moisture and groundwater monitoring. The workflows combine
multiple satellite precipitation archives (CHIRPS, GSMaP, PERSIANN-CDR,
TRMM), MODIS evapotranspiration, MERRA-2 root-zone soil moisture and a
non-parametric Mann-Kendall trend test with Sen-slope estimation, so that
the same study area can be characterised through complementary lenses
(meteorological, agricultural and hydrological drought) and screened for
statistically significant long-term change.

## Scripts

| Script | Datasets | Technique | Region |
| --- | --- | --- | --- |
| `standardized_precipitation_index.js` | CHIRPS daily | Rolling-window cumulative precipitation, day-of-year conditioned mean/std-dev, standardization for SPI-1/3/6/12/24/48 | User AOI (default Iran synoptic) |
| `precipitation_dataset_comparison.js` | CHIRPS, GSMaP v6 (reanalysis + operational), PERSIANN-CDR, TRMM 3B42 | Harmonized daily compositing, per-polygon `reduceRegions`, side-by-side time-series chart | Iran administrative polygons |
| `mann_kendall_trend_analysis.js` | MODIS MOD13A1 EVI (swappable) | Pairwise sign statistic S, Sen-slope median, tie-corrected variance, Z-score and p-value | MENA |
| `modis_evapotranspiration_drought.js` | MODIS MOD16A2 ET/PET | Monthly ET and PET composites, Evapotranspiration Stress Index ESI = ET / PET | Iran baseflow basins |
| `merra_soil_moisture_aquifers.js` | MERRA-2 land (`NASA/GSFC/MERRA/lnd/2`) | Monthly mean of root-zone moisture `RZMC` over aquifer polygons | MENA aquifers |

## How to run

1. Open <https://code.earthengine.google.com> with a Google account that
   has Earth Engine access.
2. Paste the script you want to run into a new Code Editor file.
3. Replace the `AOI` constant at the top with your own `FeatureCollection`
   asset path (or a drawn polygon).
4. Adjust the all-caps configuration block (`START_DATE`, `END_DATE`,
   `TIME_WINDOW_MONTHS`, `EXPORT_FOLDER`, etc.) to your study.
5. Click **Run**. Charts appear in the Console panel; CSV and GeoTIFF
   exports are queued in the **Tasks** tab and have to be started from
   there manually.

## Source provenance

Polished from the working scripts in `gee-scripts/ESII_2022/`,
`gee-scripts/ESI_Personal/`, `gee-scripts/ALI_SAM/`, `gee-scripts/TATAR/`
and `gee-scripts-writer/GroundWater/`. The SPI script in `ESII_2022/`
was the strongest source and forms the backbone of
`standardized_precipitation_index.js`.
