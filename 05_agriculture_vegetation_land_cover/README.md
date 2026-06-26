# 05 - Agriculture, Vegetation & Land Cover

Earth-Engine workflows for vegetation, land-cover and agricultural monitoring.
The scripts combine supervised classification (Random Forest with confusion-matrix
validation) and multi-temporal vegetation / hydrology indices, drawing on Landsat,
Sentinel-2, MODIS, AVHRR, USDA NASS CDL and ERA5-Land. They are written to
demonstrate operational pipelines: cloud / quality masking, harmonised
multi-sensor compositing, feature-stack construction, stratified sampling,
train / test splits with independent error matrices, and Drive-export-ready
outputs (rasters and CSV tables).

## Scripts

| Script | Datasets | Technique | Region |
| --- | --- | --- | --- |
| `random_forest_landcover_classification.js` | Landsat 5 TOA, MODIS MCD12Q1 (IGBP labels) | Supervised Random Forest with resubstitution + held-out validation matrices | San Francisco Bay (configurable) |
| `multi_temporal_ndvi_landsat.js` | Landsat 5 / 7 / 8 Collection 1 T1 (harmonised) | Sensor-aware NDVI, monthly median composites, annual max-NDVI stack, long-term time series | Lake Urmia basin (configurable) |
| `evapotranspiration_pet_avhrr.js` | MOD16A2 PET, TerraClimate AET, AVHRR CDR NDVI | Monthly PET / AET compositing and long-record NDVI baseline | MENA grasslands (configurable) |
| `crop_type_classification_usda_cdl.js` | USDA NASS CDL, Sentinel-2 SR (harmonised), TIGER counties | Random Forest crop-type classification, stratified sampling, county-level area summary | Illinois cropland (configurable) |
| `era5_snow_cover_timeseries.js` | ERA5-Land HOURLY (SWE, snow cover fraction) | Monthly snow water equivalent and snow-cover-fraction composites with climatology and variability maps | Alborz Mountains (configurable) |

## How to run

1. Open [code.earthengine.google.com](https://code.earthengine.google.com).
2. Create a new script and paste the contents of the desired `.js` file.
3. Adjust the all-caps configuration block at the top of the script
   (`STUDY_AREA`, `TARGET_YEAR`, `NUM_TREES_RF`, sample sizes, date windows).
4. Click **Run**. Results appear on the map canvas, charts are printed to the
   Console, and any `Export.*` tasks are queued in the **Tasks** tab.
5. For the classification scripts, the validation confusion matrix prints to
   the Console; for the time-series scripts, charts of the AOI mean are
   produced and stacks are exported to Drive (`GEE_exports/`).

## Source provenance

Derived from earlier internal GEE scripts; cleaned, renamed, and re-commented
for portfolio use. The polished versions normalise variable names, add
structured headers and section comments, expose all tuning parameters as
all-caps constants at the top of each file, remove debug `print` statements
and dead code paths, and preserve the original algorithmic choices (Random
Forest classifier hyper-parameters, compositing rules, masking thresholds and
validation strategy).
