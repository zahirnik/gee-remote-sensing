/*
 * 20-year monthly Surface Urban Heat Island (SUHI) time series
 *
 * Datasets:    MODIS/061/MOD11A1 (daily 1 km Land Surface Temperature, Terra)
 *              MODIS/006/MCD12Q1 (yearly 500 m IGBP land cover)
 *              FAO/GAUL/2015/level2 (administrative boundaries, level 2)
 * Region:      Any FAO GAUL level-2 city by name (e.g. 'Tehran', 'Greater London',
 *              'Los Angeles', 'Berlin', 'Cairo').
 * Output:      (1) FeatureCollection 'monthlyUhiTimeSeries' of monthly
 *                  SUHI (urban-mean LST minus rural-mean LST, Kelvin).
 *              (2) Map layers of land cover, urban mask and rural mask.
 *              (3) CSV export to Drive.
 *
 * SUHI definition:
 *   SUHI = mean(LST | urban land cover) - mean(LST | rural land cover)
 *   For each month we average all daily MOD11A1 acquisitions, mask the
 *   urban / rural classes from the same year's MCD12Q1 land cover, and
 *   take the mean LST in each class within the city admin boundary.
 *
 * Usage:
 *   1. Paste this script into the GEE Code Editor.
 *   2. Set CITY_NAME and the date range constants at the top.
 *   3. Click Run; chart and CSV export task are produced.
 */

// ------------------------------------------------------------------
// 0. User-editable constants
// ------------------------------------------------------------------
var CITY_NAME       = 'Tehran';                   // FAO GAUL ADM2_NAME
var START_YEAR      = 2001;                       // MCD12Q1 starts in 2001
var END_YEAR        = 2020;                       // inclusive
var MONTHLY_COUNT   = (END_YEAR - START_YEAR + 1) * 12;   // 240 months for 2001-2020

// MODIS LST is delivered in DN; multiply by 0.02 to recover Kelvin.
var LST_SCALE_FACTOR        = 0.02;

// IGBP land-cover codes from MCD12Q1 LC_Type1:
//   13 = Urban and Built-Up Lands
//   11 = Permanent Wetlands
//   15 = Snow and Ice
//   17 = Water Bodies
// We treat class 13 as 'urban'. The 'rural' reference excludes water,
// wetlands, snow/ice AND urban itself so the SUHI difference is not
// diluted by non-land or built-up surfaces.
var URBAN_LANDCOVER_CLASS   = 13;
var EXCLUDED_RURAL_CLASSES  = [11, 13, 15, 17];

var REDUCTION_SCALE_M       = 500;                // resample LST to LC native 500 m
var EXPORT_FOLDER           = 'GEE_SUHI';
var EXPORT_FILE_PREFIX      = 'SUHI_monthly_';

var IGBP_LANDCOVER_VIS = {
  min: 1, max: 17,
  palette: ['05450a', '086a10', '54a708', '78d203', '009900',
            'c6b044', 'dcd159', 'dade48', 'fbff13', 'b6ff05',
            '27ff87', 'c24f44', 'a5a5a5', 'ff6d4c', '69fff8',
            'f9ffa4', '1c0dff']
};

// ------------------------------------------------------------------
// 1. Resolve the city admin boundary from FAO GAUL level-2
// ------------------------------------------------------------------
var faoGaulLevel2 = ee.FeatureCollection('FAO/GAUL/2015/level2');
var selectedCityAdminBoundary = faoGaulLevel2
  .filter(ee.Filter.eq('ADM2_NAME', CITY_NAME));

Map.centerObject(selectedCityAdminBoundary);
Map.addLayer(selectedCityAdminBoundary.style({color: 'black', fillColor: '00000000', width: 2}),
             {}, CITY_NAME + ' admin boundary');

// ------------------------------------------------------------------
// 2. Build per-year urban / rural masks from MCD12Q1
// ------------------------------------------------------------------
var landCoverCollection = ee.ImageCollection('MODIS/006/MCD12Q1')
  .filter(ee.Filter.calendarRange(START_YEAR, END_YEAR, 'year'))
  .map(function (image) {
    return image.select('LC_Type1').clip(selectedCityAdminBoundary);
  });

// Pre-compute a 'Urban' and 'Non-Urban' binary band per yearly land-cover
// image so they can be matched to monthly LST by acquisition year below.
var landCoverWithMasks = landCoverCollection.map(function (image) {
  var urbanMask = image.eq(URBAN_LANDCOVER_CLASS).rename('Urban');

  // Build rural mask: 1 where pixel is NOT in any excluded class.
  var ruralMask = ee.Image.constant(1).rename('Non-Urban');
  EXCLUDED_RURAL_CLASSES.forEach(function (excludedClass) {
    ruralMask = ruralMask.and(image.neq(excludedClass));
  });
  ruralMask = ruralMask.rename('Non-Urban');

  return image.addBands(urbanMask).addBands(ruralMask);
});

Map.addLayer(landCoverCollection.first(), IGBP_LANDCOVER_VIS,
             'IGBP land cover ' + START_YEAR);
Map.addLayer(landCoverWithMasks.first().select('Urban').selfMask(),
             {palette: ['red']}, 'Urban mask ' + START_YEAR);

// ------------------------------------------------------------------
// 3. Build monthly composites of daily MOD11A1 LST
// ------------------------------------------------------------------
// MOD11A1 is a daily product at 1 km. For SUHI we want stable monthly
// values, so we collapse the daily collection into a monthly median.
// Median is robust to cloud-contaminated frames sneaking through QC.
var dailyLst = ee.ImageCollection('MODIS/061/MOD11A1')
  .filterDate(START_YEAR + '-01-01', (END_YEAR + 1) + '-01-01');

function buildMonthlyCollection(collection, startDateString, monthCount) {
  var monthIndexList = ee.List.sequence(0, monthCount - 1);
  var originDate = ee.Date(startDateString);

  return ee.ImageCollection(monthIndexList.map(function (monthIndex) {
    var startDate = originDate.advance(ee.Number(monthIndex), 'month');
    var endDate   = startDate.advance(1, 'month');
    return collection.filterDate(startDate, endDate).median()
      .set('system:time_start', startDate.millis())
      .set('system:time_end',   endDate.millis());
  }));
}

var monthlyLstRaw = buildMonthlyCollection(
  dailyLst, START_YEAR + '-01-01', MONTHLY_COUNT
);

// Convert the monthly composites from DN to Kelvin and resample to the
// 500 m land-cover grid so reduceRegion can pair LST and LC pixel-for-pixel.
var monthlyLst = monthlyLstRaw.map(function (image) {
  return image.select('LST_Day_1km')
    .multiply(LST_SCALE_FACTOR)
    .set('date', ee.Date(image.get('system:time_start')))
    .reproject(image.projection(), null, REDUCTION_SCALE_M);
});

// ------------------------------------------------------------------
// 4. Pair each monthly LST image with that year's urban/rural masks
// ------------------------------------------------------------------
var monthlyLstWithLandCover = monthlyLst.map(function (image) {
  var year = ee.Date(image.get('date')).get('year');
  var yearlyLandCover = landCoverWithMasks
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .first();
  return image.addBands(yearlyLandCover);
});

// ------------------------------------------------------------------
// 5. Compute monthly SUHI = mean(urban LST) - mean(rural LST)
// ------------------------------------------------------------------
var monthlyUhiTimeSeries = monthlyLstWithLandCover.map(function (image) {
  var urbanLstMean = image.select('LST_Day_1km')
    .mask(image.select('Urban'))
    .reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: selectedCityAdminBoundary,
      scale: REDUCTION_SCALE_M,
      maxPixels: 1e10
    });

  var ruralLstMean = image.select('LST_Day_1km')
    .mask(image.select('Non-Urban'))
    .reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: selectedCityAdminBoundary,
      scale: REDUCTION_SCALE_M,
      maxPixels: 1e10
    });

  var suhiKelvin = ee.Number(urbanLstMean.get('LST_Day_1km'))
    .subtract(ee.Number(ruralLstMean.get('LST_Day_1km')));

  return ee.Feature(null, {
    suhi:       suhiKelvin,
    urban_lst:  urbanLstMean.get('LST_Day_1km'),
    rural_lst:  ruralLstMean.get('LST_Day_1km'),
    date:       ee.Date(image.get('date')).format('yyyy-MM')
  });
});

// ------------------------------------------------------------------
// 6. Chart and CSV export
// ------------------------------------------------------------------
var suhiChart = ui.Chart.feature.byFeature(monthlyUhiTimeSeries, 'date', ['suhi'])
  .setChartType('LineChart')
  .setOptions({
    title: 'Monthly SUHI - ' + CITY_NAME + ' (' + START_YEAR + '-' + END_YEAR + ')',
    hAxis: {title: 'Date'},
    vAxis: {title: 'SUHI (urban - rural LST, K)'},
    lineWidth: 1,
    pointSize: 2
  });
print(suhiChart);

Export.table.toDrive({
  collection:     monthlyUhiTimeSeries,
  description:    'Export_SUHI_monthly_' + CITY_NAME,
  folder:         EXPORT_FOLDER,
  fileNamePrefix: EXPORT_FILE_PREFIX + CITY_NAME,
  fileFormat:     'CSV'
});
