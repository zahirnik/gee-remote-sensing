# 04 - Fire, Heat & Urban Climate

Thermal remote sensing pipelines built in Google Earth Engine. The category
covers four complementary uses of thermal-infrared and land-cover data: (1)
active-fire hotspot detection from MODIS FIRMS, (2) high-resolution Land
Surface Temperature (LST) retrieval from Landsat's thermal band with a
proper NDVI-based emissivity correction, (3) Surface Urban Heat Island
(SUHI) time series stratified by IGBP land cover within FAO GAUL admin
boundaries, (4) day-vs-night LST sampling across Terra and Aqua MODIS, and
(5) supervised Landsat classification of urban land cover with an LST
feature.

## Scripts

| Script | Datasets | Technique | Region |
| --- | --- | --- | --- |
| [`firms_active_fire_detection.js`](firms_active_fire_detection.js) | FIRMS (MODIS Terra/Aqua thermal anomalies) | Daily AOI clipping + seasonal max-brightness / max-confidence compositing + confidence thresholding | Configurable AOI (default: Australia 2019-2020 fire season) |
| [`landsat_surface_temperature_timeseries.js`](landsat_surface_temperature_timeseries.js) | Landsat 8 SR (B10 thermal), Landsat 8 SR (B4/B5 for NDVI) | Brightness-temperature scaling, scene-wide NDVI min/max -> fractional vegetation cover -> Sobrino-style emissivity -> single-channel LST | Configurable polygon (example: NW Iran) |
| [`surface_urban_heat_island_monthly.js`](surface_urban_heat_island_monthly.js) | MODIS MOD11A1 daily LST, MODIS MCD12Q1 yearly IGBP land cover, FAO GAUL level-2 | Daily LST collapsed to monthly median, paired with that year's land cover, SUHI = mean(urban LST) - mean(rural LST) over the city admin boundary | Any FAO GAUL ADM2 city (default: Tehran, 2001-2020) |
| [`modis_lst_day_night_seasonal.js`](modis_lst_day_night_seasonal.js) | MODIS MOD11A1 (Terra) + MYD11A1 (Aqua), day and night LST | Per-station mean LST extracted from four platform-overpass combinations, with view-time and QC bands, exported as CSV for diurnal/seasonal analysis | Configurable ROI + station points (example: Iranian plateau) |
| [`urban_land_cover_classification_landsat.js`](urban_land_cover_classification_landsat.js) | Landsat 5/7/8 TOA, SRTM (DEM/slope/aspect), Landsat SMW LST (Ermida et al. 2020), 21 spectral indices | Stratified per-class polygon sampling + Random Forest (50 trees) + OOB error + error-matrix validation + variable-importance chart | Isfahan (Iran), 1985 / 1993 / 2000 / 2008 / 2013 / 2019 epochs |

## How to run

1. Open the [GEE Code Editor](https://code.earthengine.google.com).
2. Paste the contents of a script into a new file.
3. Adjust the all-caps user constants at the top:
   - `CITY_NAME`, `STUDY_AREA`, `ROI` to localize the analysis.
   - `START_DATE` / `END_DATE` or `START_YEAR` / `END_YEAR` for the time window.
   - Asset paths (`STUDY_AREA_ASSET`, `TRAINING_ASSET`, `STATION_POINTS_ASSET`)
     to point at your own ingested polygons / points.
4. Click **Run**. Map layers and console charts appear immediately; CSV /
   GeoTIFF exports queue under the **Tasks** tab and have to be confirmed
   to start.

Some scripts depend on community modules that need to be authorized once
per Google account:

- `users/fitoprincipe/geetools:batch` (FIRMS batch export)
- `users/sofiaermida/landsat_smw_lst:modules/Landsat_LST.js` (Landsat SMW LST,
  required by `urban_land_cover_classification_landsat.js`)

## Source provenance

Derived from earlier internal GEE scripts; cleaned, renamed, and re-commented
for portfolio use. The polished versions drop debug prints, remove
commented-out variants, rename variables to descriptive English
(`dailyLst`, `urbanMask`, `monthlyUhiTimeSeries`, ...), and add all-caps
constants and section / per-line comments explaining the non-obvious physics
(scale factors, emissivity model, SUHI definition).
