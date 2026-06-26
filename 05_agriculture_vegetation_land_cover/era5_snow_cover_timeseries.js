/*
 * ERA5-Land Snow Water Equivalent and Snow Cover Time Series
 *
 * Datasets:    ECMWF/ERA5_LAND/HOURLY  (snow_depth_water_equivalent, snow_cover)
 *              Resolution ~0.1 deg (~9 km), 1950-present.
 * Region:      Configurable AOI (default: Alborz Mountains, northern Iran)
 * Output:      Monthly SWE and snow-cover-fraction ImageCollections,
 *              regional mean charts, and a Drive export of monthly SWE bands.
 *
 * Usage:
 *   1. Paste into the GEE Code Editor.
 *   2. Set STUDY_AREA, START_DATE and END_DATE.
 *   3. Click Run; monthly stacks are built, two charts print to the Console
 *      and a snow-water-equivalent export task is queued.
 */

// -------------------------------------------------------------------------
// 1. User configuration
// -------------------------------------------------------------------------
var STUDY_AREA   = ee.Geometry.Rectangle([50.0, 35.5, 53.5, 37.0]);  // Alborz Mountains
var START_DATE   = '1981-01-01';
var END_DATE     = '2020-12-31';
var REDUCER_SCALE_METRES = 9000;     // Native ERA5-Land grid
var EXPORT_SCALE_M       = 9000;
var TARGET_CRS           = 'EPSG:4326';

// -------------------------------------------------------------------------
// 2. Load hourly ERA5-Land snow variables
//    snow_depth_water_equivalent (m of water) and snow_cover (%) are the two
//    canonical snow diagnostics in ERA5-Land.
// -------------------------------------------------------------------------
var era5LandSnow = ee.ImageCollection('ECMWF/ERA5_LAND/HOURLY')
                     .filterDate(START_DATE, END_DATE)
                     .filterBounds(STUDY_AREA)
                     .select(['snow_depth_water_equivalent', 'snow_cover']);

// -------------------------------------------------------------------------
// 3. Monthly aggregation
//    For SWE we take the monthly mean (an end-of-month snapshot is also reasonable);
//    for snow cover we keep the mean fraction across the month.
// -------------------------------------------------------------------------
function monthlySnowComposite(collection, startDateString, endDateString) {
  var start = ee.Date(startDateString);
  var end   = ee.Date(endDateString);
  var monthsCount = end.difference(start, 'month').round();

  var sequence = ee.List.sequence(0, monthsCount.subtract(1));

  return ee.ImageCollection(sequence.map(function (i) {
    var monthStart = start.advance(ee.Number(i),         'month');
    var monthEnd   = start.advance(ee.Number(i).add(1),  'month');
    return collection.filterDate(monthStart, monthEnd).mean()
             .set({
               'system:time_start': monthStart.millis(),
               'system:time_end':   monthEnd.millis(),
               'year':              monthStart.get('year'),
               'month':             monthStart.get('month')
             });
  }));
}

var monthlySnowStack = monthlySnowComposite(era5LandSnow, START_DATE, END_DATE);

// Split the dual-band stack into two single-variable ImageCollections.
var monthlySnowWaterEquivalent = monthlySnowStack
                                   .select('snow_depth_water_equivalent')
                                   .map(function (img) {
                                     return img.rename('SWE_m');
                                   });

var monthlySnowCoverFraction = monthlySnowStack
                                 .select('snow_cover')
                                 .map(function (img) {
                                   return img.rename('snow_cover_pct');
                                 });

// -------------------------------------------------------------------------
// 4. Long-term mean and standard deviation maps
// -------------------------------------------------------------------------
var sweMean   = monthlySnowWaterEquivalent.mean().clip(STUDY_AREA);
var sweStdDev = monthlySnowWaterEquivalent.reduce(ee.Reducer.stdDev()).clip(STUDY_AREA);
var snowCoverClimatology = monthlySnowCoverFraction.mean().clip(STUDY_AREA);

// -------------------------------------------------------------------------
// 5. Regional time-series charts
// -------------------------------------------------------------------------
var sweChart = ui.Chart.image.series({
  imageCollection: monthlySnowWaterEquivalent,
  region:          STUDY_AREA,
  reducer:         ee.Reducer.mean(),
  scale:           REDUCER_SCALE_METRES,
  xProperty:       'system:time_start'
}).setOptions({
  title: 'ERA5-Land monthly snow water equivalent (m)',
  hAxis: {title: 'Date'},
  vAxis: {title: 'SWE (m water)'},
  series: {0: {color: '#08519c', lineWidth: 1.2, pointSize: 1.5}}
});
print(sweChart);

var snowCoverChart = ui.Chart.image.series({
  imageCollection: monthlySnowCoverFraction,
  region:          STUDY_AREA,
  reducer:         ee.Reducer.mean(),
  scale:           REDUCER_SCALE_METRES,
  xProperty:       'system:time_start'
}).setOptions({
  title: 'ERA5-Land monthly snow cover fraction (%)',
  hAxis: {title: 'Date'},
  vAxis: {title: 'Snow cover (%)'},
  series: {0: {color: '#54278f', lineWidth: 1.2, pointSize: 1.5}}
});
print(snowCoverChart);

// -------------------------------------------------------------------------
// 6. Map visualisation
// -------------------------------------------------------------------------
var sweVis = {min: 0, max: 0.3,
              palette: ['#f7fbff', '#9ecae1', '#4292c6', '#08519c', '#08306b']};
var snowCoverVis = {min: 0, max: 100,
                    palette: ['#ffffff', '#bdbdbd', '#6baed6', '#08306b']};

Map.centerObject(STUDY_AREA, 8);
Map.addLayer(sweMean,                 sweVis,        'Mean SWE (m water) - climatology');
Map.addLayer(sweStdDev,
             {min: 0, max: 0.2,
              palette: ['#fff5eb', '#fdae6b', '#d94801']},
             'SWE inter-annual variability (stdDev)');
Map.addLayer(snowCoverClimatology,    snowCoverVis,  'Snow-cover climatology (%)');

// -------------------------------------------------------------------------
// 7. Export monthly SWE stack to Drive
// -------------------------------------------------------------------------
Export.image.toDrive({
  image:       monthlySnowWaterEquivalent.toBands(),
  description: 'monthly_SWE_' + START_DATE.substring(0, 4) +
               '_' + END_DATE.substring(0, 4),
  folder:      'GEE_exports',
  region:      STUDY_AREA,
  scale:       EXPORT_SCALE_M,
  crs:         TARGET_CRS,
  maxPixels:   1e13
});
