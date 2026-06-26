# 01 — Surface Water & Lake Dynamics

This category bundles Google Earth Engine (GEE) scripts that quantify how open
water bodies and their drainage basins behave through time: lake area from
multi-decadal Landsat archives, basin-scale water-balance forcings
(precipitation, evapotranspiration, soil moisture), salt-crust extent over
shrinking endorheic lakes, and basin-aggregated hydroclimate covariates
(precipitation, temperature, snow). Together they cover the full toolkit a
hydrologist or remote-sensing analyst needs to monitor lake dynamics, attribute
change to climate vs. land-use drivers, and feed downstream models.

## Scripts

| Script | Datasets | Technique | Region |
| --- | --- | --- | --- |
| `lake_area_jrc_global_surface_water.js` | JRC GSW Yearly History (v1.3) | Per-pixel area summation of permanent-water class, stacked across years and reduced per lake polygon | Any global lake FeatureCollection |
| `multi_lake_water_balance_modis_chirps.js` | CHIRPS daily precip, MOD16A2 ET, SMAP root-zone soil moisture, ALOS AW3D30 DSM | Hydrological-year aggregation + DSM-based lake-surface masking + zonal statistics over basins | User-supplied lake / basin polygons (HydroLAKES-style) |
| `otsu_water_threshold_landsat.js` | Landsat 5/7/8 C01 SR (harmonised) | Otsu between-class-variance thresholding on NIR-band histogram for water/non-water classification | Demo: Pyramid Lake (NV, USA); retargetable |
| `salt_lake_optical_indices.js` | Landsat 5/7/8 C01 SR, JRC GSW Monthly History | Visible-band Salinity Index `SI = (Green + Red)/2`, binned salt-crust area time series | Urmia Lake bed (NW Iran); applicable to Aral Sea, Caspian, Great Salt Lake |
| `river_basin_hydroclimate_aggregation.js` | CHIRPS precip, ERA5-Land 2 m temperature, MOD10A1 NDSI snow, FLDAS SnowCover, HydroSHEDS level-3 basins | Monthly compositing + permanent-water masking + per-basin zonal statistics in long-format CSV | Demo: upper Columbia / Snake River headwaters; retargetable |

## How to run

Each script is a self-contained `.js` file ready to paste into the
[GEE Code Editor](https://code.earthengine.google.com/). The first ~15 lines
of every file are an ALL_CAPS constants block (e.g. `STUDY_AREA`,
`START_DATE`, `END_DATE`, `WATER_THRESHOLD`, `LAKE_COLLECTION`). Edit those
values to retarget the script, then click **Run**:

- Map outputs render in the GEE Map panel.
- Chart outputs render in the Console panel.
- Tabular outputs appear in the **Tasks** tab and must be started manually
  to export the CSV to Google Drive.

Asset paths are written as `users/<YOUR_USERNAME>/<YOUR_ASSET>` placeholders;
replace them with your own ingested lake / basin FeatureCollections (for
example HydroLAKES polygons, HydroBASINS level-3/4 catchments, or a manually
digitised salt-lake bed).

## Source provenance

Derived from earlier internal GEE scripts; cleaned, renamed, and re-commented
for portfolio use.
