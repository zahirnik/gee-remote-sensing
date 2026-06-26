/*
 * Optical Flood Extraction with Sentinel-2 (NDWI / MNDWI Change Detection)
 *
 * Datasets:    Sentinel-2 MSI Level-1C (B3 Green, B8 NIR, B11 SWIR-1),
 *              JRC Global Surface Water (seasonality)
 * Region:      Lake Urmia basin (NW Iran) - spring 2019 floods
 * Output:      Binary flood / inundation mask from optical change detection,
 *              NDWI and MNDWI pre / post composites
 *
 * Method:
 *   1. Build cloud-screened pre- and post-event Sentinel-2 composites.
 *   2. Compute NDWI = (Green - NIR) / (Green + NIR) and
 *      MNDWI = (Green - SWIR1) / (Green + SWIR1) for both composites.
 *   3. Flag pixels that are non-water before and water after as new water.
 *   4. Remove permanent water bodies (JRC seasonality >= 10 months).
 *
 * Usage:
 *   1. Paste into the GEE Code Editor.
 *   2. Set AOI, PRE_EVENT_RANGE, POST_EVENT_RANGE.
 *   3. Click Run; the optical flood mask is added to the map; export queued.
 */

// ------------------------------------------------------------------
// 0. User parameters
// ------------------------------------------------------------------
var AOI = ee.Geometry.Polygon([[
  [45.22954558261228, 37.13801934100997],
  [45.89187258505374, 36.6423766377166],
  [46.45526378386106, 37.20481167237601],
  [45.85758093337326, 37.50955338899794]
]]);

var PRE_EVENT_RANGE  = ['2018-04-01', '2018-05-01'];
var POST_EVENT_RANGE = ['2019-04-01', '2019-05-01'];

var CLOUD_PERCENT_MAX        = 40;     // discard scenes cloudier than this
var NDWI_WATER_THRESHOLD     = 0.0;    // NDWI > 0 commonly used to flag water
var MNDWI_WATER_THRESHOLD    = 0.0;    // MNDWI > 0 typical for mixed land/water
var PERMANENT_WATER_MONTHS   = 10;

// ------------------------------------------------------------------
// 1. Helper: cloud-screened Sentinel-2 composite for a date range
// ------------------------------------------------------------------
function buildSentinel2Composite(dateRange) {
  return ee.ImageCollection('COPERNICUS/S2')
    .filterBounds(AOI)
    .filterDate(dateRange[0], dateRange[1])
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_PERCENT_MAX))
    .median()
    .clip(AOI);
}

// ------------------------------------------------------------------
// 2. Helper: NDWI and MNDWI bands
// ------------------------------------------------------------------
function addWaterIndices(image) {
  var ndwi  = image.normalizedDifference(['B3', 'B8']).rename('NDWI');
  var mndwi = image.normalizedDifference(['B3', 'B11']).rename('MNDWI');
  return image.addBands([ndwi, mndwi]);
}

// ------------------------------------------------------------------
// 3. Build pre / post composites with water indices
// ------------------------------------------------------------------
var sentinel2Pre  = addWaterIndices(buildSentinel2Composite(PRE_EVENT_RANGE));
var sentinel2Post = addWaterIndices(buildSentinel2Composite(POST_EVENT_RANGE));

var ndwiPre  = sentinel2Pre.select('NDWI');
var ndwiPost = sentinel2Post.select('NDWI');

var mndwiPre  = sentinel2Pre.select('MNDWI');
var mndwiPost = sentinel2Post.select('MNDWI');

// ------------------------------------------------------------------
// 4. Change-detection flood mask
// ------------------------------------------------------------------
// A pixel is "newly inundated" when it was dry (water index below threshold)
// before the event and wet (above threshold) after. We require agreement of
// both NDWI and MNDWI to reduce false positives from vegetation moisture.
var waterPre  = ndwiPre.gt(NDWI_WATER_THRESHOLD).or(mndwiPre.gt(MNDWI_WATER_THRESHOLD));
var waterPost = ndwiPost.gt(NDWI_WATER_THRESHOLD).and(mndwiPost.gt(MNDWI_WATER_THRESHOLD));

var floodMaskOptical = waterPost.and(waterPre.not()).rename('flood');

// Remove permanent water bodies so we report only event-driven inundation
var jrcSeasonality = ee.Image('JRC/GSW1_0/GlobalSurfaceWater').select('seasonality');
var permanentWaterMask = jrcSeasonality.gte(PERMANENT_WATER_MONTHS);
floodMaskOptical = floodMaskOptical.where(permanentWaterMask, 0).updateMask(floodMaskOptical);

// ------------------------------------------------------------------
// 5. Compute new-water area (hectares)
// ------------------------------------------------------------------
var floodPixelArea = floodMaskOptical.multiply(ee.Image.pixelArea());
var floodAreaStats = floodPixelArea.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: AOI,
  scale: 10,
  bestEffort: true
});
var floodAreaHectares = floodAreaStats.getNumber('flood').divide(10000).round();
print('Optical-derived new water area (hectares):', floodAreaHectares);

// ------------------------------------------------------------------
// 6. Visualisation
// ------------------------------------------------------------------
Map.centerObject(AOI, 8);
var ndwiPalette = {min: -1, max: 1, palette: ['a3a3a3', 'e9ff34', '53db46', '6699CC']};

Map.addLayer(ndwiPre,           ndwiPalette,                'NDWI Pre-event',      false);
Map.addLayer(ndwiPost,          ndwiPalette,                'NDWI Post-event',     false);
Map.addLayer(mndwiPre,          ndwiPalette,                'MNDWI Pre-event',     false);
Map.addLayer(mndwiPost,         ndwiPalette,                'MNDWI Post-event',    false);
Map.addLayer(floodMaskOptical,  {palette: '0000FF'},        'Optical Flood Mask',  true);

// ------------------------------------------------------------------
// 7. Export
// ------------------------------------------------------------------
Export.image.toDrive({
  image: floodMaskOptical,
  description: 'Optical_Flood_Mask_Sentinel2',
  fileNamePrefix: 'optical_flood_mask',
  scale: 10,
  region: AOI,
  maxPixels: 1e10
});
