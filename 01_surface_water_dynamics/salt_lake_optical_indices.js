/*
 * Salinity index time series over a salt / endorheic lake (Urmia Lake by default)
 *
 * Datasets:    LANDSAT/LT05/C01/T1_SR, LANDSAT/LE07/C01/T1_SR, LANDSAT/LC08/C01/T1_SR
 *              JRC/GSW1_2/MonthlyHistory (for the water-body mask)
 * Region:      Urmia Lake bed and salt-flats (NW Iran). Easily retargetable to
 *              the Aral Sea, Great Salt Lake, Caspian shoreline, etc.
 * Output:      Four time-series charts of salt-crust area (km^2) per SI bin
 *              (0.2-0.4, 0.4-0.6, 0.6-0.8, 0.8-1.0) over the chosen window.
 *
 * Usage:
 *   1. Paste this script into the GEE Code Editor.
 *   2. Replace STUDY_AREA_SALT_LAKE with the polygon of your salt-lake bed.
 *   3. Adjust DATE_RANGE and SI_BINS as needed.
 *   4. Click Run; four LineChart prints appear in the Console.
 *
 * Method:
 *   The Salinity Index used here is the simple visible-band formulation
 *       SI = (Green + Red) / 2
 *   evaluated in surface-reflectance units (0-1). High SI means a bright,
 *   highly reflective surface in the visible -- typical of dry salt crusts.
 *   We bin SI into four 0.2-wide intervals and report total area per bin.
 */

// ------------------------------------------------------------------
// 0. User-editable constants
// ------------------------------------------------------------------
var STUDY_AREA_SALT_LAKE = ee.FeatureCollection('users/<YOUR_USERNAME>/<YOUR_SALT_LAKE_POLYGON>');

var START_DATE = '2000-01-01';
var END_DATE   = '2020-12-31';
var COMPOSITE_INTERVAL_MONTHS = 1;
var COMPOSITE_COUNT = 252;            // = (END_DATE - START_DATE) in months

// Bin edges in surface-reflectance units. Each entry produces one chart.
var SALINITY_INDEX_BINS = [
  { label: '0.2-0.4', min: 0.2, max: 0.4 },
  { label: '0.4-0.6', min: 0.4, max: 0.6 },
  { label: '0.6-0.8', min: 0.6, max: 0.8 },
  { label: '0.8-1.0', min: 0.8, max: 1.0 }
];

// ------------------------------------------------------------------
// 1. Cloud / shadow / snow masking helpers (Landsat C01 SR)
// ------------------------------------------------------------------
var maskLandsat457 = function (image) {
  var qa = image.select('pixel_qa');
  var cloud = qa.bitwiseAnd(1 << 5)
    .and(qa.bitwiseAnd(1 << 7))
    .or(qa.bitwiseAnd(1 << 3))
    .and(qa.bitwiseAnd(1 << 4));
  var edgeMask = image.mask().reduce(ee.Reducer.min());
  return image.updateMask(cloud.not()).updateMask(edgeMask);
};

var maskLandsat8 = function (image) {
  var qa = image.select('pixel_qa');
  var mask = qa.bitwiseAnd(1 << 3).eq(0)
    .and(qa.bitwiseAnd(1 << 5).eq(0))
    .and(qa.bitwiseAnd(1 << 4).eq(0));
  return image.updateMask(mask);
};

// ------------------------------------------------------------------
// 2. Build a harmonised Landsat 5/7/8 surface-reflectance collection
// ------------------------------------------------------------------
var landsat5 = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')
  .map(maskLandsat457)
  .select(['B7', 'B6', 'B5', 'B4', 'B3', 'B2', 'B1'])
  .map(function (image) { return image.clip(STUDY_AREA_SALT_LAKE); });

var landsat7 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
  .map(maskLandsat457)
  .select(['B7', 'B6', 'B5', 'B4', 'B3', 'B2', 'B1'])
  .map(function (image) { return image.clip(STUDY_AREA_SALT_LAKE); });

// L8 has a different band layout (B1 is coastal aerosol, B2 blue, ...). We rename
// to a stable [B7..B1] order so the SI expression below sees the same bands
// regardless of sensor.
var landsat8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
  .map(maskLandsat8)
  .map(function (image) { return image.clip(STUDY_AREA_SALT_LAKE); })
  .map(function (image) {
    return image.rename(['B0', 'B1', 'B2', 'B3', 'B4', 'B5',
                         'B6', 'B7', 'B8', 'B9', 'B10', 'B11']);
  })
  .select(['B7', 'B6', 'B5', 'B4', 'B3', 'B2', 'B1']);

var mergedLandsat = ee.ImageCollection(
  landsat5.merge(landsat7).merge(landsat8)
).filterDate(START_DATE, END_DATE);

// ------------------------------------------------------------------
// 3. Composite into a regular monthly time series
// ------------------------------------------------------------------
// Median over each calendar month rejects residual cloud and edge effects
// and produces a temporally regular collection that is easy to chart.
var buildMonthlyComposites = function (collection, startDate, count, interval, units) {
  var sequence = ee.List.sequence(0, ee.Number(count).subtract(1));
  var origin = ee.Date(startDate);
  return ee.ImageCollection(sequence.map(function (i) {
    var monthStart = origin.advance(ee.Number(interval).multiply(i), units);
    var monthEnd = origin.advance(ee.Number(interval).multiply(ee.Number(i).add(1)), units);
    return collection.filterDate(monthStart, monthEnd).median()
      .set('system:time_start', monthStart.millis())
      .set('system:time_end', monthEnd.millis());
  }));
};

var monthlyComposites = buildMonthlyComposites(
  mergedLandsat, START_DATE, COMPOSITE_COUNT, COMPOSITE_INTERVAL_MONTHS, 'month'
);

// Drop empty composites (months with no Landsat overpasses produce images
// that have lost all bands after merge). We tag each composite with its
// remaining band count and keep only the fully-banded ones.
var fullyBandedComposites = monthlyComposites
  .map(function (image) { return image.set('bandCount', image.bandNames().length()); })
  .filter(ee.Filter.eq('bandCount', 7));

// ------------------------------------------------------------------
// 4. Compute the salinity index SI = (Green + Red) / 2
// ------------------------------------------------------------------
// Landsat SR values are stored as int16 scaled by 1e4. We divide by 10000
// so SI is in proper reflectance units (0-1) before binning.
var salinityIndexCollection = fullyBandedComposites.map(function (image) {
  var reflectance = image.divide(10000.0);
  var salinityIndex = reflectance.expression(
    '(GREEN + RED) / 2',
    {
      'GREEN': reflectance.select('B2'),
      'RED':   reflectance.select('B3')
    }
  ).rename('SI');
  return salinityIndex.copyProperties(image, ['system:time_start', 'system:time_end']);
});

// ------------------------------------------------------------------
// 5. For each SI bin, convert per-pixel membership to area (km^2)
//    and chart it as a time series over the salt-lake polygon
// ------------------------------------------------------------------
SALINITY_INDEX_BINS.forEach(function (bin) {
  var binAreaCollection = salinityIndexCollection.map(function (image) {
    var inBin = image.gte(bin.min).and(image.lt(bin.max));
    var areaKm2 = inBin.multiply(ee.Image.pixelArea()).divide(1e6);
    return areaKm2.copyProperties(image, ['system:time_start', 'system:time_end']);
  });

  var saltCrustAreaChart = ui.Chart.image.seriesByRegion({
    imageCollection: binAreaCollection,
    scale: 30,
    regions: STUDY_AREA_SALT_LAKE,
    reducer: ee.Reducer.sum(),
    xProperty: 'system:time_start'
  })
  .setChartType('LineChart')
  .setOptions({
    title:    'Salt-crust area, SI bin ' + bin.label,
    hAxis:    { title: 'Date' },
    vAxis:    { title: 'Area (km^2)' },
    lineWidth: 1,
    pointSize: 3,
    series:   { 0: { color: 'black' } }
  });

  print(saltCrustAreaChart);
});
