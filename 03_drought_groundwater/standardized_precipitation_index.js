/*
 * Standardized Precipitation Index (SPI) from CHIRPS for 1, 3, 6, 12, 24, 48-month windows
 *
 * Datasets:    UCSB-CHG/CHIRPS/DAILY (daily precipitation, 0.05 deg, 1981-present)
 * Region:      User-defined AOI (FeatureCollection)
 * Output:      Image collection of SPI maps + AOI-mean time-series chart
 *              + optional GeoTIFF export of SPI rasters per time step
 *
 * Method:      McKee et al. (1993). For each rolling window of W months the
 *              cumulative precipitation P_w is computed at every pixel. The
 *              long-term mean (mu) and standard deviation (sigma) of P_w are
 *              derived from the full CHIRPS archive, conditioned on the
 *              calendar position of the window (so that, e.g., June-SPI-3 is
 *              standardized against historical Apr-May-Jun totals only).
 *              SPI is then (P_w - mu) / sigma. Negative SPI = drier than
 *              normal; positive SPI = wetter than normal.
 *
 * Usage:
 *   1. Paste into the GEE Code Editor.
 *   2. Adjust AOI, START_DATE, END_DATE and TIME_WINDOW_MONTHS at the top.
 *   3. Click Run; SPI map and chart appear in the console; CSV/GeoTIFF
 *      exports are queued in the Tasks tab if EXPORT_RASTERS = true.
 *
 * Caveat:      A normal distribution is assumed for the standardization step.
 *              Precipitation is more strictly modeled with a gamma fit; SPI
 *              values here should be read as a first-order estimator.
 */

// ---------------------------------------------------------------------------
// USER CONFIGURATION
// ---------------------------------------------------------------------------

// Area of interest -- replace with your own FeatureCollection asset or polygon.
var AOI = ee.FeatureCollection('users/<YOUR_GEE_USERNAME>/Iran_Snoptics');

// Time window (in months) over which precipitation is accumulated before
// standardization. Supported values: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 24, 48.
var TIME_WINDOW_MONTHS = 3;

// Native CHIRPS resolution is ~5550 m at the equator. For large AOIs you can
// coarsen this to speed up reductions.
var RESOLUTION_METERS = 5550;

// Optional GeoTIFF export window. Set EXPORT_RASTERS = true to queue tasks.
var EXPORT_RASTERS    = false;
var EXPORT_START_DATE = '2019-01-01';
var EXPORT_END_DATE   = '2023-12-31';
var EXPORT_FOLDER     = 'GEE_SPI_Exports';

// SPI visualization: red = dry, blue = wet.
var SPI_VIS = {
  bands:   ['SPI'],
  min:     -4,
  max:      4,
  palette: ['d53e4f', 'fc8d59', 'fee08b', 'ffffbf', 'e6f598', '99d594', '3288bd']
};

// ---------------------------------------------------------------------------
// LOAD CHIRPS AND ESTABLISH THE TIME AXIS
// ---------------------------------------------------------------------------

var CHIRPS = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY');

var firstImage  = ee.Date(ee.List(CHIRPS.get('date_range')).get(0));
var latestImage = ee.Date(
  CHIRPS.limit(1, 'system:time_start', false).first().get('system:time_start')
);

print('CHIRPS available from', firstImage,  'to', latestImage);
print('SPI window (months):',  TIME_WINDOW_MONTHS);

// ---------------------------------------------------------------------------
// BUILD A LIST OF END-OF-WINDOW DATES
// ---------------------------------------------------------------------------
// We walk backwards from the most recent CHIRPS image in steps of
// TIME_WINDOW_MONTHS, then sort ascending. Each entry marks the END of a
// rolling accumulation window.
//
// Example for TIME_WINDOW_MONTHS = 3 and latestImage = 2024-05-31:
//   window 1: 2024-03-01 -> 2024-05-31
//   window 2: 2023-12-01 -> 2024-02-29
//   window 3: 2023-09-01 -> 2023-11-30
//   ... etc.

var totalSteps = latestImage
  .difference(firstImage, 'month')
  .divide(ee.Number(TIME_WINDOW_MONTHS));

var stepIndices = ee.List.sequence(0, totalSteps);

var windowEndDates = stepIndices.map(function (stepIndex) {
  var monthsBack = ee.Number(stepIndex).multiply(TIME_WINDOW_MONTHS).multiply(-1);
  // Advance one day past the latest CHIRPS image so the inclusive sum captures it.
  return latestImage.advance(1, 'day').advance(monthsBack, 'month');
}).sort();

// ---------------------------------------------------------------------------
// AGGREGATE CHIRPS INTO ROLLING-WINDOW PRECIPITATION TOTALS
// ---------------------------------------------------------------------------
// For each end date, sum the CHIRPS daily rasters over the preceding
// TIME_WINDOW_MONTHS months. We only retain windows that span the full
// requested duration (early windows can be partially clipped by firstImage).

var cumulativePrecipPerWindow = ee.ImageCollection.fromImages(
  windowEndDates.map(function (windowEnd) {
    var endDate   = ee.Date(windowEnd);
    var startDate = endDate.advance(-TIME_WINDOW_MONTHS, 'month');

    var dailyInWindow = CHIRPS
      .filterDate(startDate, endDate)
      .map(function (img) { return img.clip(AOI); });

    var summedPrecip = dailyInWindow.sum().set({
      'system:time_start': startDate.millis(),
      'system:time_end':   endDate.millis(),
      'window_months':     TIME_WINDOW_MONTHS,
      'image_count':       dailyInWindow.size()
    });

    // Skip incomplete windows at the start of the record.
    var observedMonths = endDate.difference(startDate, 'month').round();
    return ee.Algorithms.If(
      observedMonths.gte(TIME_WINDOW_MONTHS),
      summedPrecip
    );
  })
);

var precipitationCollection = ee.ImageCollection(
  cumulativePrecipPerWindow.copyProperties(CHIRPS)
);

// ---------------------------------------------------------------------------
// COMPUTE PER-WINDOW LONG-TERM MEAN AND STANDARD DEVIATION
// ---------------------------------------------------------------------------
// For sub-annual SPI (< 12 months) we standardize against all historical
// windows that share the same day-of-year footprint -- i.e. June-SPI-3 is
// only compared against other Apr-May-Jun blocks. For SPI windows >= 12
// months the calendar position is irrelevant and we standardize against
// the whole archive.

function standardizeSubAnnual (precipCol) {
  return precipCol.map(function (windowImage) {
    var startDOY = ee.Date(windowImage.get('system:time_start')).getRelative('day', 'year');
    var endDOY   = ee.Date(windowImage.get('system:time_end')).getRelative('day', 'year');

    var matchingWindows = precipCol
      .filter(ee.Filter.calendarRange(startDOY, endDOY, 'day_of_year'))
      .reduce(ee.Reducer.stdDev().combine(ee.Reducer.mean(), null, true));

    var spi = windowImage.expression(
      '(precip - mu) / sigma',
      {
        precip: windowImage.select('precipitation'),
        mu:     matchingWindows.select('precipitation_mean'),
        sigma:  matchingWindows.select('precipitation_stdDev')
      }
    ).rename('SPI');

    return windowImage.addBands(matchingWindows).addBands(spi);
  });
}

function standardizeMultiAnnual (precipCol) {
  var stats = precipCol.reduce(ee.Reducer.stdDev().combine(ee.Reducer.mean(), null, true));

  return precipCol.map(function (windowImage) {
    var spi = windowImage.expression(
      '(precip - mu) / sigma',
      {
        precip: windowImage.select('precipitation'),
        mu:     stats.select('precipitation_mean'),
        sigma:  stats.select('precipitation_stdDev')
      }
    ).rename('SPI');

    return windowImage.addBands(stats).addBands(spi);
  });
}

var THRESHOLD_MONTHS = 12;

var spiCollection = ee.ImageCollection(
  ee.Algorithms.If(
    ee.Number(TIME_WINDOW_MONTHS).gte(THRESHOLD_MONTHS),
    standardizeMultiAnnual(precipitationCollection),
    standardizeSubAnnual(precipitationCollection)
  )
);

// ---------------------------------------------------------------------------
// MAP DISPLAY
// ---------------------------------------------------------------------------

Map.centerObject(AOI, 6);

var mostRecentSpi = spiCollection.limit(1, 'system:time_start', false).first();
Map.addLayer(
  mostRecentSpi.select('SPI'),
  SPI_VIS,
  'SPI-' + TIME_WINDOW_MONTHS + ' (most recent)'
);

// ---------------------------------------------------------------------------
// TIME-SERIES CHARTS
// ---------------------------------------------------------------------------
// Average SPI and average cumulative precipitation across the AOI for every
// rolling window, plotted against the window start time.

var aoiLabeled = AOI.map(function (feature) {
  return feature
    .set('label_spi',    'SPI-' + TIME_WINDOW_MONTHS)
    .set('label_precip', TIME_WINDOW_MONTHS + '-month precip [mm]');
});

var precipChart = ui.Chart.image.seriesByRegion(
  precipitationCollection, aoiLabeled, ee.Reducer.mean(),
  'precipitation', RESOLUTION_METERS, 'system:time_start', 'label_precip'
).setOptions({
  title:  TIME_WINDOW_MONTHS + '-month cumulative precipitation (CHIRPS, AOI mean)',
  vAxis:  { title: 'Precipitation [mm]' },
  hAxis:  { title: 'Year' }
});
print(precipChart);

var spiChart = ui.Chart.image.seriesByRegion(
  spiCollection, aoiLabeled, ee.Reducer.mean(),
  'SPI', RESOLUTION_METERS, 'system:time_start', 'label_spi'
).setOptions({
  title:  'SPI-' + TIME_WINDOW_MONTHS + ' time series (CHIRPS, AOI mean)',
  vAxis:  { title: 'SPI [-]' },
  hAxis:  { title: 'Year' }
});
print(spiChart);

// ---------------------------------------------------------------------------
// OPTIONAL RASTER EXPORTS
// ---------------------------------------------------------------------------
// Each SPI image is exported as a single-band GeoTIFF, named by its window
// start date. The exports go to Drive under EXPORT_FOLDER.

if (EXPORT_RASTERS) {
  var exportCollection = spiCollection
    .select(['SPI'])
    .filterDate(EXPORT_START_DATE, EXPORT_END_DATE);

  var exportList = exportCollection.toList(exportCollection.size());
  var exportSize = exportCollection.size().getInfo();

  for (var i = 0; i < exportSize; i++) {
    var spiImage = ee.Image(exportList.get(i));
    var dateTag  = ee.Date(spiImage.get('system:time_start')).format('YYYY-MM-dd').getInfo();

    Export.image.toDrive({
      image:       spiImage,
      description: 'SPI_' + TIME_WINDOW_MONTHS + 'mo_' + dateTag,
      folder:      EXPORT_FOLDER,
      region:      AOI.geometry().bounds(),
      scale:       RESOLUTION_METERS,
      maxPixels:   1e13
    });
  }
}
