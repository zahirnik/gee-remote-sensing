/*
 * Flood Damage Assessment - Tabular Reporting
 *
 * Datasets:    Sentinel-1 GRD (VH) flood mask,
 *              JRC GHSL Population Density 2015 (250 m),
 *              MODIS MCD12Q1 IGBP Land Cover (500 m),
 *              JRC Global Surface Water,
 *              WWF HydroSHEDS DEM
 * Region:      Beira, Mozambique - Cyclone Idai event (March 2019)
 * Output:      Summary FeatureCollection exported as CSV / GeoJSON containing
 *              flooded area, exposed population, affected cropland and urban
 *              area for the AOI.
 *
 * Method:
 *   1. Build a SAR-derived flood mask (same pipeline as the SAR script).
 *   2. Intersect with GHSL population to count exposed people.
 *   3. Intersect with MODIS cropland (classes 12 & 14) and urban (class 13)
 *      to compute affected agricultural and built-up area.
 *   4. Pack all numbers into a single Feature and export to Drive.
 *
 * Usage:
 *   1. Paste into the GEE Code Editor.
 *   2. Set AOI, BEFORE_DATE_RANGE, AFTER_DATE_RANGE and the EVENT_NAME label.
 *   3. Click Run. Two export tasks (CSV + GeoJSON) appear in the Tasks panel.
 */

// ------------------------------------------------------------------
// 0. User parameters
// ------------------------------------------------------------------
var EVENT_NAME = 'CycloneIdai_Beira_2019';

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

var MODIS_CLASS_CROPLAND        = 12;
var MODIS_CLASS_CROPLAND_MOSAIC = 14;
var MODIS_CLASS_URBAN           = 13;

// ------------------------------------------------------------------
// 1. SAR flood mask
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

var jrcSeasonality = ee.Image('JRC/GSW1_0/GlobalSurfaceWater').select('seasonality');
var permanentWaterMask = jrcSeasonality.gte(PERMANENT_WATER_MONTHS)
  .updateMask(jrcSeasonality.gte(PERMANENT_WATER_MONTHS));

var floodMask = floodCandidate.where(permanentWaterMask, 0);
floodMask = floodMask.updateMask(floodMask);
floodMask = floodMask.updateMask(floodMask.connectedPixelCount().gte(MIN_CONNECTED_PIXELS));

var terrainSlope = ee.Algorithms.Terrain(ee.Image('WWF/HydroSHEDS/03VFDEM')).select('slope');
floodMask = floodMask.updateMask(terrainSlope.lt(SLOPE_MASK_DEGREES));

// ------------------------------------------------------------------
// 2. Total flooded area (hectares)
// ------------------------------------------------------------------
var floodPixelArea = floodMask.select(SAR_POLARIZATION).multiply(ee.Image.pixelArea());
var floodAreaHectares = ee.Number(floodPixelArea.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: AOI,
  scale: 10,
  bestEffort: true
}).get(SAR_POLARIZATION)).divide(10000).round();

// ------------------------------------------------------------------
// 3. Exposed population (GHSL 250 m)
// ------------------------------------------------------------------
var populationCount = ee.Image('JRC/GHSL/P2016/POP_GPW_GLOBE_V1/2015').clip(AOI);
var floodMaskOnGhslGrid = floodMask.reproject({crs: populationCount.projection()});

var populationExposed = populationCount
  .updateMask(floodMaskOnGhslGrid)
  .updateMask(populationCount);

var populationExposedTotal = ee.Number(populationExposed.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: AOI,
  scale: 250,
  maxPixels: 1e9
}).get('population_count')).round();

// ------------------------------------------------------------------
// 4. Affected cropland and urban area (MODIS 500 m)
// ------------------------------------------------------------------
var modisLandCover = ee.ImageCollection('MODIS/006/MCD12Q1')
  .filterDate('2014-01-01', AFTER_DATE_RANGE[1])
  .sort('system:index', false)
  .select('LC_Type1')
  .first()
  .clip(AOI);

var floodMaskOnModisGrid = floodMask.reproject({crs: modisLandCover.projection()});

var croplandMask = modisLandCover.eq(MODIS_CLASS_CROPLAND)
  .or(modisLandCover.eq(MODIS_CLASS_CROPLAND_MOSAIC));
var croplandAffected = floodMaskOnModisGrid.updateMask(croplandMask);
var croplandAffectedHectares = ee.Number(
  croplandAffected.multiply(ee.Image.pixelArea())
    .reduceRegion({reducer: ee.Reducer.sum(), geometry: AOI, scale: 500, maxPixels: 1e9})
    .get(SAR_POLARIZATION)
).divide(10000).round();

var urbanMask = modisLandCover.eq(MODIS_CLASS_URBAN);
var urbanAffected = modisLandCover.updateMask(urbanMask).mask(floodMaskOnModisGrid)
  .updateMask(modisLandCover.updateMask(urbanMask));
var urbanAffectedHectares = ee.Number(
  urbanAffected.multiply(ee.Image.pixelArea())
    .reduceRegion({reducer: ee.Reducer.sum(), geometry: AOI, scale: 500, bestEffort: true})
    .get('LC_Type1')
).divide(10000).round();

// ------------------------------------------------------------------
// 5. Assemble report Feature
// ------------------------------------------------------------------
var reportFeature = ee.Feature(AOI.centroid(10), {
  event_name:                  EVENT_NAME,
  before_start:                BEFORE_DATE_RANGE[0],
  before_end:                  BEFORE_DATE_RANGE[1],
  after_start:                 AFTER_DATE_RANGE[0],
  after_end:                   AFTER_DATE_RANGE[1],
  sar_polarization:            SAR_POLARIZATION,
  sar_pass_direction:          SAR_PASS_DIRECTION,
  sar_difference_threshold:    SAR_DIFFERENCE_THRESHOLD,
  flooded_area_hectares:       floodAreaHectares,
  exposed_population_people:   populationExposedTotal,
  affected_cropland_hectares:  croplandAffectedHectares,
  affected_urban_hectares:     urbanAffectedHectares
});

var reportCollection = ee.FeatureCollection([reportFeature]);
print('Damage assessment report:', reportFeature);

// ------------------------------------------------------------------
// 6. Exports (CSV + GeoJSON)
// ------------------------------------------------------------------
Export.table.toDrive({
  collection: reportCollection,
  description: 'Flood_Damage_Report_CSV',
  fileFormat: 'CSV',
  fileNamePrefix: 'flood_damage_report_' + EVENT_NAME
});

Export.table.toDrive({
  collection: reportCollection,
  description: 'Flood_Damage_Report_GeoJSON',
  fileFormat: 'GeoJSON',
  fileNamePrefix: 'flood_damage_report_' + EVENT_NAME
});
