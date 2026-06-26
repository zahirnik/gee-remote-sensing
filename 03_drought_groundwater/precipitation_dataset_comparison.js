/*
 * Multi-source precipitation dataset comparison (CHIRPS / GSMaP / PERSIANN-CDR / TRMM)
 *
 * Datasets:    UCSB-CHG/CHIRPS/DAILY            (CHIRPS, 1981-present, 0.05 deg)
 *              JAXA/GPM_L3/GSMaP/v6/reanalysis  (GSMaP reanalysis, 2000-2014)
 *              JAXA/GPM_L3/GSMaP/v6/operational (GSMaP operational, 2014-present)
 *              NOAA/PERSIANN-CDR                (PERSIANN-CDR, 1983-present)
 *              TRMM/3B42                        (TRMM 3B42, 1998-2019)
 * Region:      Iran synoptic stations / administrative polygons
 * Output:      One-row-per-polygon-per-image CSV table per product +
 *              AOI-mean time-series chart overlaying the four products
 *
 * Usage:
 *   1. Paste into the GEE Code Editor.
 *   2. Adjust AOI, START_DATE, END_DATE constants at the top.
 *   3. Click Run; charts appear in the console; CSV export tasks are queued.
 */

// ---------------------------------------------------------------------------
// USER CONFIGURATION
// ---------------------------------------------------------------------------

var AOI            = ee.FeatureCollection('users/<YOUR_GEE_USERNAME>/Iran_Snoptics');
var POLYGON_ID     = 'ST_CODE';   // attribute used as polygon identifier in the exported CSV

var START_DATE     = '2000-06-01';
var END_DATE       = '2014-03-01';

var EXPORT_FOLDER  = 'GEE_Precip_Comparison';
var EXPORT_PREFIX  = 'precip_';

// Bounding box that brackets all CHIRPS / GSMaP / PERSIANN / TRMM filterBounds
// calls -- only used to filter the input collections, not to clip output.
var BBOX = ee.Geometry.Rectangle([41.0, 21.6, 71.8, 41.5], null, false);

// ---------------------------------------------------------------------------
// LOAD AND HARMONIZE THE FOUR PRECIPITATION ARCHIVES
// ---------------------------------------------------------------------------
// Each product reports precipitation with its own band name, native time step
// and unit (mm/day vs mm/hour). We aggregate each one to a common daily mean
// image collection and rename the variable to 'precipitation' (mm) so the
// downstream code is identical for every product.

function dailyComposite (collection, startDate, endDate, reducerName) {
  var origin     = ee.Date(startDate);
  var numDays    = ee.Date(endDate).difference(origin, 'days');
  var daySeq     = ee.List.sequence(0, numDays.subtract(1));

  return ee.ImageCollection(daySeq.map(function (offset) {
    var dayStart = origin.advance(ee.Number(offset),       'days');
    var dayEnd   = origin.advance(ee.Number(offset).add(1), 'days');

    var window = collection.filterDate(dayStart, dayEnd);
    var image  = reducerName === 'sum' ? window.sum() : window.mean();

    return image
      .rename('precipitation')
      .set('system:time_start', dayStart.millis());
  }));
}

// --- CHIRPS (already daily mm/day) ---
var chirpsDaily = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .filterDate(START_DATE, END_DATE)
  .filterBounds(BBOX)
  .select(['precipitation'], ['precipitation']);

// --- GSMaP reanalysis (hourly mm/h, aggregate to daily mean) ---
var gsmapRaw = ee.ImageCollection('JAXA/GPM_L3/GSMaP/v6/reanalysis')
  .filterDate(START_DATE, END_DATE)
  .filterBounds(BBOX)
  .select(['hourlyPrecipRateGC'], ['precipitation']);
var gsmapDaily = dailyComposite(gsmapRaw, START_DATE, END_DATE, 'mean');

// --- PERSIANN-CDR (already daily mm/day) ---
var persiannDaily = ee.ImageCollection('NOAA/PERSIANN-CDR')
  .filterDate(START_DATE, END_DATE)
  .filterBounds(BBOX)
  .select(['precipitation'], ['precipitation']);

// --- TRMM 3B42 (3-hourly mm/h, aggregate to daily sum) ---
var trmmRaw = ee.ImageCollection('TRMM/3B42')
  .filterDate(START_DATE, END_DATE)
  .filterBounds(BBOX)
  .select(['precipitation'], ['precipitation']);
var trmmDaily = dailyComposite(trmmRaw, START_DATE, END_DATE, 'sum');

// ---------------------------------------------------------------------------
// REDUCE TO POLYGON-LEVEL TRIPLETS (polygon x image x mean precipitation)
// ---------------------------------------------------------------------------

function buildTriplets (collection, scaleMeters) {
  return collection.map(function (image) {
    return image.reduceRegions({
      collection: AOI.select([POLYGON_ID]),
      reducer:    ee.Reducer.mean(),
      scale:      scaleMeters
    }).map(function (feature) {
      return feature
        .set('imageId',   image.id())
        .set('imageDate', image.date());
    });
  }).flatten();
}

var chirpsTriplets   = buildTriplets(chirpsDaily,   5000);
var gsmapTriplets    = buildTriplets(gsmapDaily,   10000);
var persiannTriplets = buildTriplets(persiannDaily, 25000);
var trmmTriplets     = buildTriplets(trmmDaily,    25000);

// ---------------------------------------------------------------------------
// AOI-MEAN TIME SERIES COMPARISON CHART
// ---------------------------------------------------------------------------
// To put the four products on one axis we average each daily image across
// the entire AOI, label it with the product name and merge the four into a
// single FeatureCollection. ui.Chart.feature.groups then plots one line per
// product.

function aoiDailyTimeSeries (collection, scaleMeters, productName) {
  return collection.map(function (image) {
    var meanDict = image.reduceRegion({
      reducer:   ee.Reducer.mean(),
      geometry:  AOI.geometry(),
      scale:     scaleMeters,
      maxPixels: 1e10
    });
    return ee.Feature(null, {
      product: productName,
      date:    image.date().millis(),
      precip:  meanDict.get('precipitation')
    });
  });
}

var combinedSeries = aoiDailyTimeSeries(chirpsDaily,    5000, 'CHIRPS')
  .merge(aoiDailyTimeSeries(gsmapDaily,   10000, 'GSMaP'))
  .merge(aoiDailyTimeSeries(persiannDaily, 25000, 'PERSIANN-CDR'))
  .merge(aoiDailyTimeSeries(trmmDaily,    25000, 'TRMM-3B42'));

var comparisonChart = ui.Chart.feature.groups(combinedSeries, 'date', 'precip', 'product')
  .setChartType('LineChart')
  .setOptions({
    title:  'AOI-mean daily precipitation -- CHIRPS vs GSMaP vs PERSIANN-CDR vs TRMM',
    hAxis:  { title: 'Date',            format: 'YYYY-MM' },
    vAxis:  { title: 'Precipitation [mm/day]' },
    lineWidth: 1
  });
print(comparisonChart);

// ---------------------------------------------------------------------------
// EXPORTS: one CSV per product
// ---------------------------------------------------------------------------

Export.table.toDrive({
  collection:  chirpsTriplets,
  description: EXPORT_PREFIX + 'CHIRPS',
  folder:      EXPORT_FOLDER,
  fileFormat:  'CSV'
});

Export.table.toDrive({
  collection:  gsmapTriplets,
  description: EXPORT_PREFIX + 'GSMaP',
  folder:      EXPORT_FOLDER,
  fileFormat:  'CSV'
});

Export.table.toDrive({
  collection:  persiannTriplets,
  description: EXPORT_PREFIX + 'PERSIANN_CDR',
  folder:      EXPORT_FOLDER,
  fileFormat:  'CSV'
});

Export.table.toDrive({
  collection:  trmmTriplets,
  description: EXPORT_PREFIX + 'TRMM_3B42',
  folder:      EXPORT_FOLDER,
  fileFormat:  'CSV'
});
