/*
 * Multi-year lake-area time series from JRC Global Surface Water (Yearly History)
 *
 * Datasets:    JRC/GSW1_3/YearlyHistory (Pekel et al. 2016)
 * Region:      Any global lake FeatureCollection (HydroLAKES-style polygons)
 * Output:      One CSV row per lake feature with annual permanent-water area (km^2)
 *              for every year in the analysis window, exported to Google Drive.
 *
 * Usage:
 *   1. Paste this script into the GEE Code Editor.
 *   2. Adjust the LAKE_COLLECTION asset path and DATE_RANGE constants at the top.
 *   3. Click Run; the task panel will show an Export-to-Drive task to start.
 */

// ------------------------------------------------------------------
// 0. User-editable constants
// ------------------------------------------------------------------
// Asset path to a FeatureCollection of lake polygons. Each feature should
// carry a unique identifier so rows can later be joined back to attributes.
// Replace with your own ingested lake-polygon asset.
var LAKE_COLLECTION    = 'users/<YOUR_USERNAME>/<YOUR_LAKE_POLYGONS>';
var LAKE_ID_PROPERTY   = 'GLWD_ID';     // unique-id field copied through to output
var START_YEAR         = 1995;
var END_YEAR           = 2020;          // inclusive
// JRC GSW YearlyHistory band 'waterClass' encoding:
//   0 no data, 1 not water, 2 seasonal water, 3 permanent water.
// >= 3 keeps only permanent water; use >= 2 to include seasonal water as well.
var WATER_CLASS_THRESHOLD = 3;
var REDUCTION_SCALE_M  = 30;            // native Landsat / JRC GSW resolution
var OUTPUT_FILENAME    = 'lake_area_jrc_gsw_yearly';

// ------------------------------------------------------------------
// 1. Load lake polygons and the JRC GSW yearly history collection
// ------------------------------------------------------------------
var lakeCollection = ee.FeatureCollection(LAKE_COLLECTION);

var gswYearly = ee.ImageCollection('JRC/GSW1_3/YearlyHistory')
  .filterDate(START_YEAR + '-01-01', (END_YEAR + 1) + '-01-01');

// ------------------------------------------------------------------
// 2. Convert each annual classification image to a per-pixel area (km^2)
// ------------------------------------------------------------------
// Multiplying a 0/1 water mask by ee.Image.pixelArea() preserves true area
// per pixel under the image's native projection (handles latitude distortion).
var annualWaterAreaCollection = gswYearly.map(function (classifiedImage) {
  var waterMask = classifiedImage.gte(WATER_CLASS_THRESHOLD);
  var areaKm2 = waterMask.multiply(ee.Image.pixelArea()).divide(1e6);
  return areaKm2.copyProperties(classifiedImage, ['system:time_start', 'system:time_end']);
});

// ------------------------------------------------------------------
// 3. Stack the per-year area images into a single multi-band image
// ------------------------------------------------------------------
// Stacking into bands lets us run a single reduceRegions() per lake instead
// of one reduce per year, which is much faster on big lake collections.
var yearList = ee.List.sequence(START_YEAR, END_YEAR);
var imageList = annualWaterAreaCollection.toList(annualWaterAreaCollection.size());

var stackedAnnualArea = ee.Image(
  yearList.iterate(function (year, accumulator) {
    var index = ee.Number(year).subtract(START_YEAR);
    var bandName = ee.Number(year).format('%d');
    var yearImage = ee.Image(imageList.get(index)).rename(bandName);
    return ee.Image(accumulator).addBands(yearImage);
  }, ee.Image().select())   // start with an empty image so the first addBands works
);

// ------------------------------------------------------------------
// 4. For each lake, sum the water-area pixels for every year
// ------------------------------------------------------------------
var extractLakeAreaTimeSeries = function (lakeFeature) {
  var perYearSums = stackedAnnualArea.reduceRegions({
    collection: ee.FeatureCollection([lakeFeature]),
    reducer: ee.Reducer.sum(),
    scale: REDUCTION_SCALE_M
  });

  // Roll up the (single-row) FeatureCollection into a flat dictionary
  // keyed by year, then attach to a new output feature.
  var outputFeature = ee.Feature(null);
  outputFeature = ee.Feature(
    yearList.iterate(function (year, accumulator) {
      var bandName = ee.Number(year).format('%d');
      var value = perYearSums.aggregate_sum(bandName);
      return ee.Feature(accumulator).set(bandName, value);
    }, outputFeature)
  );

  return outputFeature.copyProperties(lakeFeature, [LAKE_ID_PROPERTY]);
};

var lakeAreaTimeSeries = lakeCollection.map(extractLakeAreaTimeSeries);

// ------------------------------------------------------------------
// 5. Export the resulting table to Google Drive as CSV
// ------------------------------------------------------------------
Export.table.toDrive({
  collection: lakeAreaTimeSeries,
  description: OUTPUT_FILENAME,
  fileFormat: 'CSV'
});
