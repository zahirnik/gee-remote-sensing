/*
 * Flood Exposure Assessment: Population & Land Cover
 *
 * Datasets:    Sentinel-1 GRD (VH, IW) for flood mask,
 *              JRC GHSL Population Density 2015 (250 m),
 *              MODIS MCD12Q1 IGBP Land Cover (500 m),
 *              JRC Global Surface Water (seasonality)
 * Region:      Beira, Mozambique - Cyclone Idai event (March 2019)
 * Output:      Exposed-population raster, affected cropland mask,
 *              affected urban mask, console summary of exposure totals
 *
 * Method:
 *   1. Build a binary flood mask from Sentinel-1 before/after differencing.
 *   2. Reproject the flood mask onto the GHSL grid (250 m) and intersect
 *      with population density to count exposed people.
 *   3. Reproject the flood mask onto the MODIS grid (500 m) and intersect
 *      with IGBP cropland (classes 12 & 14) and urban (class 13) classes.
 *
 * Usage:
 *   1. Paste into the GEE Code Editor.
 *   2. Set AOI, BEFORE_DATE_RANGE, AFTER_DATE_RANGE.
 *   3. Click Run; summary numbers print to the console; raster exports queue.
 */

// ------------------------------------------------------------------
// 0. User parameters
// ------------------------------------------------------------------
var AOI = ee.Geometry.Polygon([[
  [35.53377589953368, -19.6674648789114],
  [34.50106105578368, -18.952058786515526],
  [33.63314113390868, -19.87423907259203],
  [34.74825343859618, -20.61123742951084]
]]);

var BEFORE_DATE_RANGE = ['2019-03-10', '2019-03-31'];
var AFTER_DATE_RANGE  = ['2019-04-11', '2019-04-20'];

var SAR_POLARIZATION         = 'VH';
var SAR_PASS_DIRECTION       = 'DESCENDING';
var SAR_DIFFERENCE_THRESHOLD = 1.25;
var SPECKLE_SMOOTHING_RADIUS = 50;
var MIN_CONNECTED_PIXELS     = 8;
var SLOPE_MASK_DEGREES       = 5;
var PERMANENT_WATER_MONTHS   = 10;

// MODIS IGBP class codes
var MODIS_CLASS_CROPLAND        = 12; // pure cropland
var MODIS_CLASS_CROPLAND_MOSAIC = 14; // cropland / natural vegetation mosaic
var MODIS_CLASS_URBAN           = 13; // urban / built-up

// ------------------------------------------------------------------
// 1. Build Sentinel-1 flood mask (compact pipeline)
// ------------------------------------------------------------------
var sentinel1Collection = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', SAR_POLARIZATION))
  .filter(ee.Filter.eq('orbitProperties_pass', SAR_PASS_DIRECTION))
  .filter(ee.Filter.eq('resolution_meters', 10))
  .filterBounds(AOI)
  .select(SAR_POLARIZATION);

var sarBeforeFiltered = sentinel1Collection.filterDate(BEFORE_DATE_RANGE[0], BEFORE_DATE_RANGE[1])
  .mosaic().clip(AOI)
  .focal_mean(SPECKLE_SMOOTHING_RADIUS, 'circle', 'meters');

var sarAfterFiltered = sentinel1Collection.filterDate(AFTER_DATE_RANGE[0], AFTER_DATE_RANGE[1])
  .mosaic().clip(AOI)
  .focal_mean(SPECKLE_SMOOTHING_RADIUS, 'circle', 'meters');

var floodCandidate = sarAfterFiltered.divide(sarBeforeFiltered).gt(SAR_DIFFERENCE_THRESHOLD);

// Refine: remove permanent water, isolated pixels, and steep slopes
var jrcSeasonality = ee.Image('JRC/GSW1_0/GlobalSurfaceWater').select('seasonality');
var permanentWaterMask = jrcSeasonality.gte(PERMANENT_WATER_MONTHS)
  .updateMask(jrcSeasonality.gte(PERMANENT_WATER_MONTHS));

var floodMask = floodCandidate.where(permanentWaterMask, 0);
floodMask = floodMask.updateMask(floodMask);
floodMask = floodMask.updateMask(floodMask.connectedPixelCount().gte(MIN_CONNECTED_PIXELS));

var terrainSlope = ee.Algorithms.Terrain(ee.Image('WWF/HydroSHEDS/03VFDEM')).select('slope');
floodMask = floodMask.updateMask(terrainSlope.lt(SLOPE_MASK_DEGREES));

// ------------------------------------------------------------------
// 2. Population exposure (GHSL 250 m)
// ------------------------------------------------------------------
var populationCount = ee.Image('JRC/GHSL/P2016/POP_GPW_GLOBE_V1/2015').clip(AOI);
var ghslProjection = populationCount.projection();

// Reproject flood mask onto the GHSL grid so masking is pixel-aligned
var floodMaskOnGhslGrid = floodMask.reproject({crs: ghslProjection});

var populationExposed = populationCount
  .updateMask(floodMaskOnGhslGrid)
  .updateMask(populationCount);

var populationExposedStats = populationExposed.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: AOI,
  scale: 250,
  maxPixels: 1e9
});

var populationExposedTotal = populationExposedStats.getNumber('population_count').round();
print('Exposed population (people):', populationExposedTotal);

// ------------------------------------------------------------------
// 3. Land-cover exposure (MODIS MCD12Q1, 500 m)
// ------------------------------------------------------------------
var modisLandCover = ee.ImageCollection('MODIS/006/MCD12Q1')
  .filterDate('2014-01-01', AFTER_DATE_RANGE[1])
  .sort('system:index', false)
  .select('LC_Type1')
  .first()
  .clip(AOI);

var modisProjection = modisLandCover.projection();
var floodMaskOnModisGrid = floodMask.reproject({crs: modisProjection});

// --- Cropland (classes 12 & 14) ---
var croplandMask = modisLandCover.eq(MODIS_CLASS_CROPLAND)
  .or(modisLandCover.eq(MODIS_CLASS_CROPLAND_MOSAIC));
var croplandLayer = modisLandCover.updateMask(croplandMask);
var croplandAffected = floodMaskOnModisGrid.updateMask(croplandMask);

var croplandAffectedPixelArea = croplandAffected.multiply(ee.Image.pixelArea());
var croplandAffectedStats = croplandAffectedPixelArea.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: AOI,
  scale: 500,
  maxPixels: 1e9
});
var croplandAffectedHectares = croplandAffectedStats.getNumber(SAR_POLARIZATION).divide(10000).round();
print('Affected cropland (hectares):', croplandAffectedHectares);

// --- Urban (class 13) ---
var urbanMask = modisLandCover.eq(MODIS_CLASS_URBAN);
var urbanLayer = modisLandCover.updateMask(urbanMask);
var urbanAffected = urbanLayer.mask(floodMaskOnModisGrid).updateMask(urbanLayer);

var urbanAffectedPixelArea = urbanAffected.multiply(ee.Image.pixelArea());
var urbanAffectedStats = urbanAffectedPixelArea.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: AOI,
  scale: 500,
  bestEffort: true
});
var urbanAffectedHectares = urbanAffectedStats.getNumber('LC_Type1').divide(10000).round();
print('Affected urban area (hectares):', urbanAffectedHectares);

// ------------------------------------------------------------------
// 4. Visualisation
// ------------------------------------------------------------------
Map.centerObject(AOI, 8);

var populationCountVis = {min: 0, max: 200, palette: ['060606', '337663', '337663', 'ffffff']};
var populationExposedVis = {min: 0, max: 200, palette: ['yellow', 'orange', 'red']};
var landCoverVis = {
  min: 1, max: 17,
  palette: [
    '05450a', '086a10', '54a708', '78d203', '009900', 'c6b044', 'dcd159',
    'dade48', 'fbff13', 'b6ff05', '27ff87', 'c24f44', 'a5a5a5', 'ff6d4c',
    '69fff8', 'f9ffa4', '1c0dff'
  ]
};

Map.addLayer(floodMask,         {palette: '0000FF'},               'Flood Mask');
Map.addLayer(populationCount,   populationCountVis,                'Population Density',  false);
Map.addLayer(populationExposed, populationExposedVis,              'Exposed Population');
Map.addLayer(modisLandCover,    landCoverVis,                      'MODIS Land Cover',    false);
Map.addLayer(croplandLayer,     {min: 0, max: 14, palette: ['30b21c']}, 'Cropland',       false);
Map.addLayer(croplandAffected,  {min: 0, max: 14, palette: ['30b21c']}, 'Affected Cropland');
Map.addLayer(urbanLayer,        {min: 0, max: 13, palette: ['grey']},   'Urban',          false);
Map.addLayer(urbanAffected,     {min: 0, max: 13, palette: ['grey']},   'Affected Urban');

// ------------------------------------------------------------------
// 5. Exports
// ------------------------------------------------------------------
Export.image.toDrive({
  image: populationExposed,
  description: 'Exposed_Population_Raster',
  fileNamePrefix: 'population_exposed',
  scale: 250,
  region: AOI,
  maxPixels: 1e10
});

Export.image.toDrive({
  image: croplandAffected,
  description: 'Affected_Cropland_Raster',
  fileNamePrefix: 'cropland_affected',
  scale: 500,
  region: AOI,
  maxPixels: 1e10
});

Export.image.toDrive({
  image: urbanAffected,
  description: 'Affected_Urban_Raster',
  fileNamePrefix: 'urban_affected',
  scale: 500,
  region: AOI,
  maxPixels: 1e10
});
