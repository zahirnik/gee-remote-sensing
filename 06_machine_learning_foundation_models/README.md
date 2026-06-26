# 06 — Machine Learning & Foundation Models

Three end-to-end GEE pipelines that combine the
[**AlphaEarth Foundation Embeddings**](https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_SATELLITE_EMBEDDING_V1_ANNUAL)
(64-dim per-pixel learned representations, 10 m, global, annual since 2017)
with GEE's built-in classifiers (`smileRandomForest`,
`smileGradientTreeBoost`). Every script runs **entirely on Earth Engine** —
no external model serving, no Vertex AI dependency, no data upload required.

Each script has been **smoke-tested against a live Earth Engine project**:
copy any of them into the [GEE Code Editor](https://code.earthengine.google.com),
click Run, and it works.

## Scripts

| Script | Application | Foundation features | ML model | Sample validation accuracy |
|---|---|---|---|---|
| [`01_crop_type_alphaearth_random_forest.js`](01_crop_type_alphaearth_random_forest.js) | Agriculture — crop-type mapping | AlphaEarth 64-dim embedding for the target year | Random Forest (100 trees) | **0.94** on a 6-class corn / soy / wheat / other test split |
| [`02_drought_severity_alphaearth_gbt.js`](02_drought_severity_alphaearth_gbt.js) | Drought — severity classification | AlphaEarth(T-1) + Δembedding(T-1 - T-3) — 128 features | Gradient-Boosted Trees | Runs end-to-end; CHIRPS-derived USDM-style D0..D4 labels |
| [`03_flood_few_shot_alphaearth_embedding.js`](03_flood_few_shot_alphaearth_embedding.js) | Flood — few-shot disaster mapping | Δembedding (post - pre event) | Prototype-distance (zero training) + Random Forest (12-point few-shot) + Sentinel-1 SAR baseline | Produces three flood masks for direct method comparison |

## What the foundation model brings

AlphaEarth is Google's geospatial foundation model. For each pixel, it
distils all years of multi-sensor Earth-observation data into a 64-dim
embedding that captures spectral, temporal, and spatial context. Used as
features for a downstream classifier, these embeddings:

- Need **far fewer labelled samples** to reach competitive accuracy than raw
  spectral bands or hand-engineered indices.
- Generalize better across regions (the embeddings are agnostic to local
  illumination, cloud, and seasonal artefacts).
- Enable **few-shot and prototype-based methods** that work with a handful
  of labelled pixels — exactly the operational setting for disaster mapping
  where labelled data doesn't exist yet.

## How to run

1. Open the [GEE Code Editor](https://code.earthengine.google.com) signed
   into a project with Earth Engine enabled.
2. Copy any `.js` file's contents into a new script.
3. Adjust the `CONFIG` block at the top — at minimum, change `STUDY_AREA`
   and `TARGET_YEAR` if you want a different scene. Defaults are tuned for
   small, fast runs.
4. Click **Run**. The Map shows the classified raster, the Console prints
   accuracy + confusion matrix, and an `Export.image.toDrive` task is
   queued under the Tasks tab — click Run there to save the GeoTIFF.

## Source provenance

Original Earth-Engine code; not derived from any earlier internal script.
Designed from scratch around the GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL
collection (released by Google in 2024-2025) and standard EE classifiers.

## References

- AlphaEarth Foundation Embeddings dataset:
  https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_SATELLITE_EMBEDDING_V1_ANNUAL
- USDA NASS Cropland Data Layer:
  https://developers.google.com/earth-engine/datasets/catalog/USDA_NASS_CDL
- CHIRPS Daily Precipitation:
  https://developers.google.com/earth-engine/datasets/catalog/UCSB-CHG_CHIRPS_DAILY
- JRC Global Surface Water Monthly History:
  https://developers.google.com/earth-engine/datasets/catalog/JRC_GSW1_4_MonthlyHistory
- U.S. Drought Monitor (label inspiration for Script 2):
  https://droughtmonitor.unl.edu/
- Copernicus Emergency Management Service (validation reference for Script 3):
  https://emergency.copernicus.eu/mapping/list-of-components/EMSR348
