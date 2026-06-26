# 02. Flood & Inundation Mapping

SAR-based rapid flood mapping with multi-sensor exposure assessment in Google
Earth Engine. The scripts implement a Sentinel-1 GRD change-detection pipeline
(speckle filtering, ratio thresholding, JRC permanent-water removal, slope
masking, connectivity cleanup) for the flood mask, then quantify human and
land-use exposure by intersecting that mask with JRC GHSL population density
and MODIS IGBP land cover. A Sentinel-2 NDWI / MNDWI pipeline is also provided
for optical confirmation, and a reporting script bundles all damage statistics
into a single CSV / GeoJSON ready for downstream GIS use.

## Scripts

| Script | Datasets | Technique | Output |
| --- | --- | --- | --- |
| `sar_flood_extent_sentinel1.js` | Sentinel-1 GRD (VH), JRC GSW, HydroSHEDS DEM | Before/after ratio thresholding with speckle, connectivity and slope masking | Binary flood mask (GeoTIFF + SHP), flooded area in hectares |
| `flood_exposure_population_landcover.js` | Sentinel-1 GRD, JRC GHSL Population (250 m), MODIS MCD12Q1 (500 m) | Reprojected flood mask intersected with population density, IGBP cropland and urban classes | Exposed-population raster, affected cropland & urban rasters, console summary |
| `optical_flood_extraction_sentinel2.js` | Sentinel-2 MSI L1C, JRC GSW | NDWI + MNDWI pre/post change detection (dual-index agreement) with permanent-water removal | Optical flood mask GeoTIFF, new-water area in hectares |
| `flood_damage_assessment_report.js` | Sentinel-1 GRD, JRC GHSL Population, MODIS MCD12Q1, JRC GSW, HydroSHEDS DEM | End-to-end pipeline packaging area, population and land-use statistics into a single Feature | CSV + GeoJSON damage report exported to Drive |

## How to run

1. Open the [GEE Code Editor](https://code.earthengine.google.com/).
2. Paste any single script into a new file.
3. Edit the constants at the top of the file - `AOI`, `BEFORE_DATE_RANGE` /
   `AFTER_DATE_RANGE` (or `PRE_EVENT_RANGE` / `POST_EVENT_RANGE` for the
   optical script), and the SAR/optical thresholds if you need to tune for a
   different scene.
4. Click **Run**. Map layers update immediately and any `Export.*` tasks land
   in the **Tasks** panel; click **Run** there to write outputs to Drive.

The default AOIs and dates are pre-set to the **Cyclone Idai / Beira,
Mozambique (March-April 2019)** event for the SAR and reporting scripts, and
to the **Lake Urmia basin, NW Iran (spring 2019)** for the optical script, so
each script can be run end-to-end as a demo.

## Source provenance

Derived from earlier internal GEE scripts; cleaned, renamed, and re-commented
for portfolio use. Algorithms preserve the original UN-SPIDER SAR
change-detection logic; only naming, structure and commentary have been
polished.
