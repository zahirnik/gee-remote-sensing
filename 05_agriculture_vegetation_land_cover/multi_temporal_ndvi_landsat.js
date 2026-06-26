/*
 * Multi-Temporal NDVI from Harmonised Landsat 5 / 7 / 8 (1990s - present)
 *
 * Datasets:    LANDSAT/LT05/C01/T1   (1984-2012, TM)
 *              LANDSAT/LE07/C01/T1   (1999-2022, ETM+)
 *              LANDSAT/LC08/C01/T1   (2013-present, OLI)
 * Region:      Configurable AOI (default: Lake Urmia basin, NW Iran)
 * Output:      Monthly NDVI ImageCollection, optional annual maxNDVI composites,
 *              and a time-series chart over the AOI.
 *
 * Usage:
 *   1. Paste into the GEE Code Editor.
 *   2. Set STUDY_AREA, START_YEAR and END_YEAR.
 *   3. Click Run; monthly composites are produced, the first is shown on the map
 *      and a long-term NDVI chart is printed to the Console.
 */

// -------------------------------------------------------------------------
// 1. User configuration
// -------------------------------------------------------------------------
var STUDY_AREA          = ee.Geometry.Rectangle([44.5, 36.5, 46.5, 38.5]); // Lake Urmia basin
var START_YEAR          = 1990;
var END_YEAR            = 2023;
var CHART_REDUCER_SCALE = 1000;   // Metres; coarse scale keeps the chart responsive.
var EXPORT_SCALE_METRES = 30;     // Native Landsat resolution.

// -------------------------------------------------------------------------
// 2. Harmonised, sensor-aware NDVI extractor
//    Band mapping:  TM/ETM+  red=B3, NIR=B4
//                   OLI       red=B4, NIR=B5
// -------------------------------------------------------------------------
function landsatNdvi(image, redBand, nirBand) {
  return image.normalizedDifference([nirBand, redBand]).rename('NDVI')
              .copyProperties(image, ['system:time_start']);
}

function buildHarmonisedLandsatCollection(startDate, endDate, aoi) {
  var landsat5 = ee.ImageCollection('LANDSAT/LT05/C01/T1')
                   .filterDate(startDate, endDate)
                   .filterBounds(aoi)
                   .map(function (img) { return landsatNdvi(img, 'B3', 'B4'); });

  var landsat7 = ee.ImageCollection('LANDSAT/LE07/C01/T1')
                   .filterDate(startDate, endDate)
                   .filterBounds(aoi)
                   .map(function (img) { return landsatNdvi(img, 'B3', 'B4'); });

  var landsat8 = ee.ImageCollection('LANDSAT/LC08/C01/T1')
                   .filterDate(startDate, endDate)
                   .filterBounds(aoi)
                   .map(function (img) { return landsatNdvi(img, 'B4', 'B5'); });

  return landsat5.merge(landsat7).merge(landsat8).sort('system:time_start');
}

var harmonisedNdviCollection = buildHarmonisedLandsatCollection(
  ee.Date.fromYMD(START_YEAR, 1, 1),
  ee.Date.fromYMD(END_YEAR + 1, 1, 1),
  STUDY_AREA
);

// -------------------------------------------------------------------------
// 3. Monthly compositing (median is robust against residual cloud / shadow)
// -------------------------------------------------------------------------
function monthlyComposites(collection, startYear, endYear) {
  var years  = ee.List.sequence(startYear, endYear);
  var months = ee.List.sequence(1, 12);

  var composites = years.map(function (y) {
    return months.map(function (m) {
      var start = ee.Date.fromYMD(y, m, 1);
      var end   = start.advance(1, 'month');
      var monthly = collection.filterDate(start, end).median().rename('NDVI');
      return monthly.set({
        'system:time_start': start.millis(),
        'system:time_end':   end.millis(),
        'year':              y,
        'month':             m
      });
    });
  }).flatten();

  return ee.ImageCollection.fromImages(composites);
}

var monthlyNdviCollection = monthlyComposites(harmonisedNdviCollection,
                                              START_YEAR, END_YEAR);

// -------------------------------------------------------------------------
// 4. Annual maximum-NDVI composite (good for peak-greenness mapping)
// -------------------------------------------------------------------------
var annualMaxNdvi = ee.ImageCollection(
  ee.List.sequence(START_YEAR, END_YEAR).map(function (y) {
    var start = ee.Date.fromYMD(y, 1, 1);
    var end   = start.advance(1, 'year');
    return harmonisedNdviCollection.filterDate(start, end).max().rename('NDVI_max')
             .set('system:time_start', start.millis())
             .set('year', y);
  })
);

// -------------------------------------------------------------------------
// 5. Visualise + time-series chart
// -------------------------------------------------------------------------
var ndviPalette = ['#a50026', '#f46d43', '#fee08b', '#d9ef8b', '#66bd63', '#006837'];
var firstMonth  = monthlyNdviCollection.first();

Map.centerObject(STUDY_AREA, 7);
Map.addLayer(firstMonth, {min: -0.2, max: 0.8, palette: ndviPalette},
             'NDVI - earliest monthly composite');

var ndviTimeSeriesChart = ui.Chart.image.series({
  imageCollection: monthlyNdviCollection,
  region:          STUDY_AREA,
  reducer:         ee.Reducer.mean(),
  scale:           CHART_REDUCER_SCALE,
  xProperty:       'system:time_start'
}).setOptions({
  title: 'Mean NDVI - monthly composites, ' + START_YEAR + ' to ' + END_YEAR,
  hAxis: {title: 'Date'},
  vAxis: {title: 'NDVI'},
  lineWidth: 1.5,
  pointSize: 2,
  series: {0: {color: '#1b7837'}}
});

print(ndviTimeSeriesChart);

// -------------------------------------------------------------------------
// 6. Export annual maxNDVI stack (one task; bands are years)
// -------------------------------------------------------------------------
var annualMaxStack = annualMaxNdvi.toBands();

Export.image.toDrive({
  image:        annualMaxStack,
  description:  'annual_maxNDVI_' + START_YEAR + '_' + END_YEAR,
  folder:       'GEE_exports',
  region:       STUDY_AREA,
  scale:        EXPORT_SCALE_METRES,
  maxPixels:    1e13
});
