/*
 * SAR Flood Extent Mapping with Sentinel-1 GRD (Change Detection)
 *
 * Datasets:    Sentinel-1 GRD (VH, IW, DESCENDING),
 *              JRC Global Surface Water (seasonality),
 *              WWF HydroSHEDS DEM (slope mask)
 * Region:      Beira, Mozambique - Cyclone Idai flood event (March 2019)
 * Output:      Binary flood mask (raster GeoTIFF + vector SHP),
 *              total flooded area (hectares) printed to console
 *
 * Method:
 *   Before / after Sentinel-1 GRD pairs are mosaicked, speckle-filtered, and
 *   ratio-differenced. A threshold on the ratio image yields the candidate
 *   flood mask, which is then refined with permanent-water, connectivity and
 *   slope masks to remove false positives.
 *
 * Usage:
 *   1. Paste into the GEE Code Editor.
 *   2. Adjust AOI, BEFORE_DATE_RANGE, AFTER_DATE_RANGE and the SAR thresholds.
 *   3. Click Run. The binary flood mask is added to the map and two Export
 *      tasks (raster + vector) are queued in the Tasks panel.
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

// Sentinel-1 acquisition filter
var SAR_POLARIZATION   = 'VH';          // VH is generally preferred for flood mapping over land
var SAR_PASS_DIRECTION = 'DESCENDING';  // restrict to a single pass to avoid geometry-induced artefacts

// Algorithm parameters (tuned by trial-and-error on multiple events)
var SAR_DIFFERENCE_THRESHOLD = 1.25;    // ratio threshold (after / before); pixels above are candidate flood
var SPECKLE_SMOOTHING_RADIUS = 50;      // metres, focal-mean kernel radius for speckle suppression
var MIN_CONNECTED_PIXELS     = 8;       // connectivity filter to drop isolated speckle hits
var SLOPE_MASK_DEGREES       = 5;       // mask terrain steeper than this; water seldom collects on steep slopes
var PERMANENT_WATER_MONTHS   = 10;      // JRC seasonality threshold; >=10 months/year means "permanent water"

// ------------------------------------------------------------------
// 1. Load Sentinel-1 GRD before / after collections
// ------------------------------------------------------------------
var sentinel1Collection = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', SAR_POLARIZATION))
  .filter(ee.Filter.eq('orbitProperties_pass', SAR_PASS_DIRECTION))
  .filter(ee.Filter.eq('resolution_meters', 10))
  .filterBounds(AOI)
  .select(SAR_POLARIZATION);

var sarBeforeCollection = sentinel1Collection.filterDate(BEFORE_DATE_RANGE[0], BEFORE_DATE_RANGE[1]);
var sarAfterCollection  = sentinel1Collection.filterDate(AFTER_DATE_RANGE[0],  AFTER_DATE_RANGE[1]);

// ------------------------------------------------------------------
// 2. Pre-process Sentinel-1 (mosaic + speckle filter)
// ------------------------------------------------------------------
// GRD imagery already includes thermal-noise removal, radiometric calibration
// and terrain correction, so a speckle filter is the only remaining step.
var sarBeforeMosaic = sarBeforeCollection.mosaic().clip(AOI);
var sarAfterMosaic  = sarAfterCollection.mosaic().clip(AOI);

var sarBeforeFiltered = sarBeforeMosaic.focal_mean(SPECKLE_SMOOTHING_RADIUS, 'circle', 'meters');
var sarAfterFiltered  = sarAfterMosaic.focal_mean(SPECKLE_SMOOTHING_RADIUS, 'circle', 'meters');

// ------------------------------------------------------------------
// 3. Change detection: ratio after / before
// ------------------------------------------------------------------
// A flood event darkens land (specular reflection of the new water surface),
// so the after / before ratio of linear backscatter rises above 1 over flooded
// pixels. Working in the linear domain avoids dB-difference sign confusion.
var sarRatioDifference = sarAfterFiltered.divide(sarBeforeFiltered);
var floodCandidate = sarRatioDifference.gt(SAR_DIFFERENCE_THRESHOLD);

// ------------------------------------------------------------------
// 4. Refine flood mask with permanent water, connectivity and slope
// ------------------------------------------------------------------
// (a) Remove pixels that are permanent water in the JRC climatology
var jrcSeasonality = ee.Image('JRC/GSW1_0/GlobalSurfaceWater').select('seasonality');
var permanentWaterMask = jrcSeasonality.gte(PERMANENT_WATER_MONTHS)
  .updateMask(jrcSeasonality.gte(PERMANENT_WATER_MONTHS));

var floodMask = floodCandidate.where(permanentWaterMask, 0);
floodMask = floodMask.updateMask(floodMask);

// (b) Drop isolated noisy pixels via 8-neighbour connectivity
var pixelConnectivity = floodMask.connectedPixelCount();
floodMask = floodMask.updateMask(pixelConnectivity.gte(MIN_CONNECTED_PIXELS));

// (c) Drop steep slopes (water does not accumulate on slopes > 5 degrees)
var hydroshedsDem = ee.Image('WWF/HydroSHEDS/03VFDEM');
var terrainSlope = ee.Algorithms.Terrain(hydroshedsDem).select('slope');
floodMask = floodMask.updateMask(terrainSlope.lt(SLOPE_MASK_DEGREES));

// ------------------------------------------------------------------
// 5. Compute flooded area (hectares)
// ------------------------------------------------------------------
var floodPixelArea = floodMask.select(SAR_POLARIZATION).multiply(ee.Image.pixelArea());

var floodAreaStats = floodPixelArea.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: AOI,
  scale: 10,
  bestEffort: true
});

var floodAreaHectares = floodAreaStats.getNumber(SAR_POLARIZATION).divide(10000).round();
print('Flooded area (hectares):', floodAreaHectares);

// ------------------------------------------------------------------
// 6. Visualisation
// ------------------------------------------------------------------
Map.centerObject(AOI, 8);
Map.addLayer(sarBeforeFiltered, {min: -25, max: 0}, 'Sentinel-1 Before', false);
Map.addLayer(sarAfterFiltered,  {min: -25, max: 0}, 'Sentinel-1 After',  true);
Map.addLayer(sarRatioDifference, {min: 0, max: 2}, 'Ratio (after / before)', false);
Map.addLayer(floodMask, {palette: '0000FF'}, 'Flood Mask');

// ------------------------------------------------------------------
// 7. Exports (raster + vector)
// ------------------------------------------------------------------
Export.image.toDrive({
  image: floodMask,
  description: 'SAR_Flood_Extent_Raster',
  fileNamePrefix: 'sar_flood_extent',
  region: AOI,
  scale: 10,
  maxPixels: 1e10
});

var floodVectors = floodMask.reduceToVectors({
  scale: 10,
  geometryType: 'polygon',
  geometry: AOI,
  eightConnected: false,
  bestEffort: true,
  tileScale: 2
});

Export.table.toDrive({
  collection: floodVectors,
  description: 'SAR_Flood_Extent_Vector',
  fileFormat: 'SHP',
  fileNamePrefix: 'sar_flood_extent_polygons'
});
