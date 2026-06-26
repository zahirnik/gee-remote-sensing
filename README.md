# Earth Engine Remote Sensing Portfolio

A curated collection of **24 polished Google Earth Engine scripts**, organised
into five application areas. Every script ships with a header docblock
(datasets, region, output, how-to-run), all-caps tunable constants at the top,
numbered section comments, and inline explanations of non-obvious choices —
so anyone can drop a script into the [GEE Code Editor](https://code.earthengine.google.com),
adjust a few constants, and re-run the analysis.

Total: **4,047 lines** of Earth-Engine JavaScript across five categories.

---

## Categories

### [01 — Surface water & lake dynamics](01_surface_water_dynamics/) — 5 scripts, 831 lines
Multi-decadal lake-area accounting, OTSU water classification, multi-source
water-balance, salt-lake optical indices, and river-basin hydroclimate
aggregation. Datasets: JRC Global Surface Water, Landsat 5/7/8, MODIS,
CHIRPS, SMAP, FLDAS, HydroSHEDS.

### [02 — Flood & inundation mapping](02_flood_inundation_mapping/) — 4 scripts, 661 lines
Sentinel-1 SAR change-detection flood mask, GHSL population × MODIS land
cover exposure overlay, Sentinel-2 NDWI/MNDWI optical confirmation, and a
unified damage-assessment report exporter. Default AOI: Cyclone Idai / Beira,
Mozambique 2019.

### [03 — Drought, soil moisture & groundwater](03_drought_groundwater/) — 5 scripts, 933 lines
Standardised Precipitation Index (1/3/6/12/24/48-month windows from CHIRPS),
multi-product precipitation comparison (CHIRPS vs GSMaP vs PERSIANN-CDR vs
TRMM), Mann-Kendall + Sen-slope trend analysis, MODIS MOD16 ET-based ESI,
and MERRA-2 root-zone soil moisture over aquifer polygons.

### [04 — Fire, heat & urban climate](04_fire_heat_urban_climate/) — 5 scripts, 877 lines
MODIS FIRMS active-fire detection, Landsat 8 single-channel LST with
Sobrino-style emissivity, 20-year monthly Surface Urban Heat Island
(MOD11A1 × MCD12Q1 × FAO GAUL), MODIS day-vs-night LST sampling at station
points, and Random-Forest urban land-cover classification with SMW LST,
spectral indices, and DEM as features.

### [05 — Vegetation, land cover & agriculture](05_vegetation_land_cover/) — 5 scripts, 745 lines
Supervised Random Forest land-cover classification on Landsat / Sentinel-2,
multi-temporal NDVI compositing from harmonised Landsat 5/7/8, AVHRR/MODIS
PET and ET, crop-type classification using USDA NASS CDL as labels, and
ERA5-Land snow-cover monthly time series.

---

## How to use a script

1. Click into a category folder above.
2. Open any `.js` file — the **header docblock** tells you what datasets it
   uses, what region it targets, and what it outputs.
3. Copy the file's contents into the [GEE Code Editor](https://code.earthengine.google.com).
4. Edit the **ALL-CAPS constants** at the top (study area, dates, thresholds,
   asset paths). The placeholder `users/<YOUR_GEE_USERNAME>/...` points at
   private assets you must replace with your own ingestions.
5. Click **Run**. Map layers + console charts appear immediately; CSV /
   GeoTIFF exports queue under the **Tasks** tab and need a manual confirm.

A few scripts depend on community modules that GEE will ask you to authorise
once per Google account — each script's header notes which ones.

---

## Why this layout

- **Categorised by application**, not by technique. Time-series compositing
  and zonal stats appear everywhere — they're not what differentiates these
  scripts; the *target phenomenon* (lake area, flood extent, drought, fire,
  vegetation) is.
- **Curated, not exhaustive.** These 24 scripts are the strongest of a much
  larger personal corpus. Quality over volume.
- **Human-readable code.** Variable names like `lakeAreaTimeSeries`,
  `monthlyUhiTimeSeries`, and `randomForestClassifier` instead of `lats`,
  `ut`, `clf` — long names, deliberate.
- **One self-contained file per analysis.** Each script can be run on its
  own; nothing depends on anything else in the repo. The cost of a tiny bit
  of duplication is worth the gain of obvious portability.

---

## Source provenance

All scripts are cleaned and re-organised versions of earlier internal GEE
work. Each polished script preserves the underlying scientific algorithm
unchanged — only naming, structure, and commentary have been overhauled for
clarity and reusability.

---

## License

Not yet specified — add a `LICENSE` (e.g. MIT) before broad re-use.
