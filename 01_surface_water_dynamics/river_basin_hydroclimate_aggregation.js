/*
 * Basin-scale hydroclimate aggregation: precipitation, air temperature, snow
 * cover and snow-water-equivalent zonal statistics over HydroSHEDS basins.
 *
 * Datasets:    UCSB-CHG/CHIRPS/DAILY                 (precipitation, mm/day)
 *              ECMWF/ERA5_LAND/MONTHLY                (2-m air temperature)
 *              MODIS/006/MOD10A1                      (NDSI snow cover)
 *              NASA/FLDAS/NOAH01/C/GL/M/V001          (snow water equivalent)
 *              WWF/HydroSHEDS/v1/Basins/hybas_3       (level-3 basin polygons)
 *              MODIS/MOD44W/MOD44W_005_2000_02_24     (permanent-water mask)
 * Region:      User-defined polygon clipping a sub-set of HydroSHEDS basins;
 *              default polygon covers the upper Columbia / Snake River
 *              headwaters in the US Pacific Northwest.
 * Output:      Four time-series charts (precip, temperature, NDSI snow cover,
 *              FLDAS SWE) for the selected basins, plus an Export-to-Drive
 *              task that writes a long-format CSV of all four variables per
 *              basin per month.
 *
 * Usage:
 *   1. Paste this script into the GEE Code Editor.
 *   2. Replace STUDY_AREA_POLYGON / BASIN_LEVEL constants with your region.
 *   3. Adjust DATE_RANGE constants.
 *   4. Click Run; charts render in the Console and the Tasks panel will
 *      show one Export-to-Drive task to start.
 */

// ------------------------------------------------------------------
// 0. User-editable constants
// ------------------------------------------------------------------
var STUDY_AREA_POLYGON = ee.Geometry.Polygon(
  [[[-122.006640625, 46.72833847933665],
    [-118.139453125, 42.86364317514532],
    [-114.623828125, 43.44071377076342],
    [-112.60234375,  46.48681864416886],
    [-118.139453125, 48.67986178009012]]]);

var BASIN_ASSET   = 'WWF/HydroSHEDS/v1/Basins/hybas_3';
var BASIN_ID_FIELD = 'HYBAS_ID';

var START_DATE      = '2002-03-01';
var END_DATE        = '2020-01-01';
var COMPOSITE_MONTHS = 1;                // monthly composites

// ------------------------------------------------------------------
// 1. Load basin polygons that intersect the study area
// ------------------------------------------------------------------
var selectedBasins = ee.FeatureCollection(BASIN_ASSET)
  .filterBounds(STUDY_AREA_POLYGON);

Map.centerObject(STUDY_AREA_POLYGON, 6);
Map.addLayer(selectedBasins, { color: '1f78b4' }, 'Selected HydroSHEDS basins');

// ------------------------------------------------------------------
// 2. Generic monthly-mean / monthly-sum compositor
// ------------------------------------------------------------------
// The same windowing logic is reused for every hydroclimate variable so the
// resulting collections share the system:time_start axis and can be safely
// joined later (one per-basin row per month).
var buildMonthlyComposites = function (collection, startDate, count, reducer) {
  var sequence = ee.List.sequence(0, ee.Number(count).subtract(1));
  var origin = ee.Date(startDate);
  return ee.ImageCollection(sequence.map(function (i) {
    var monthStart = origin.advance(i, 'month');
    var monthEnd = origin.advance(ee.Number(i).add(1), 'month');
    var monthlyImage = (reducer === 'sum')
      ? collection.filterDate(monthStart, monthEnd).sum()
      : collection.filterDate(monthStart, monthEnd).mean();
    return monthlyImage
      .set('system:time_start', monthStart.millis())
      .set('system:time_end',   monthEnd.millis());
  }));
};

var totalMonths = ee.Date(END_DATE).difference(ee.Date(START_DATE), 'month').round();

// ------------------------------------------------------------------
// 3. Precipitation: CHIRPS daily -> monthly totals (mm/month)
// ------------------------------------------------------------------
var dailyChirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .select('precipitation')
  .filterBounds(STUDY_AREA_POLYGON);

var monthlyPrecipitation = buildMonthlyComposites(
  dailyChirps, START_DATE, totalMonths, 'sum'
).map(function (image) { return image.rename('precipitation_mm'); });

// ------------------------------------------------------------------
// 4. Air temperature: ERA5-Land monthly 2 m temperature (K -> degC)
// ------------------------------------------------------------------
var era5Temperature = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY')
  .select('temperature_2m')
  .filterBounds(STUDY_AREA_POLYGON)
  .map(function (image) {
    // Convert kelvin to Celsius for downstream interpretability.
    return image.subtract(273.15)
      .rename('temperature_c')
      .copyProperties(image, ['system:time_start', 'system:time_end']);
  });

var monthlyTemperature = buildMonthlyComposites(
  era5Temperature, START_DATE, totalMonths, 'mean'
);

// ------------------------------------------------------------------
// 5. Snow cover: MOD10A1 NDSI -> per-pixel snow-covered area (km^2)
// ------------------------------------------------------------------
// We mask out the MOD44W permanent-water pixels first, because NDSI over
// open water is meaningless and can spuriously flip across freeze/thaw.
var permanentWaterMask = ee.Image('MODIS/MOD44W/MOD44W_005_2000_02_24')
  .select('water_mask')
  .clip(STUDY_AREA_POLYGON);

var dailyMod10 = ee.ImageCollection('MODIS/006/MOD10A1')
  .filterBounds(STUDY_AREA_POLYGON)
  .filterDate(START_DATE, END_DATE)
  .select('NDSI_Snow_Cover')
  .map(function (image) {
    var clipped = image.clip(STUDY_AREA_POLYGON);
    // Replace open-water pixels with 0 (no snow) instead of letting them
    // dominate the NDSI threshold.
    var landOnly = permanentWaterMask.where(permanentWaterMask.eq(0), clipped);
    // NDSI >= 10 (0-100 scale) is the canonical 'snow present' threshold.
    var snowPresent = landOnly.gte(10);
    var snowAreaKm2 = landOnly.updateMask(snowPresent)
      .multiply(ee.Image.pixelArea().divide(1e8));     // pixels in 100 km^2 units
    return snowAreaKm2
      .rename('snow_cover_km2')
      .copyProperties(image, ['system:time_start', 'system:time_end']);
  });

var monthlySnowCover = buildMonthlyComposites(
  dailyMod10, START_DATE, totalMonths, 'mean'
);

// ------------------------------------------------------------------
// 6. Snow water equivalent: FLDAS monthly snow cover instantaneous (kg/m^2)
// ------------------------------------------------------------------
var fldasSnow = ee.ImageCollection('NASA/FLDAS/NOAH01/C/GL/M/V001')
  .filterDate(START_DATE, END_DATE)
  .filterBounds(STUDY_AREA_POLYGON)
  .select('SnowCover_inst')
  .map(function (image) {
    // Multiply by pixel area and divide by 1e6 to get km^2 snow-covered area
    // (assuming SnowCover_inst is a 0-1 fractional cover).
    var areaKm2 = image.multiply(ee.Image.pixelArea()).divide(1e6);
    return areaKm2
      .rename('fldas_snow_km2')
      .copyProperties(image, ['system:time_start', 'system:time_end']);
  });

// ------------------------------------------------------------------
// 7. Time-series charts of each variable over the study area
// ------------------------------------------------------------------
print(ui.Chart.image.series(monthlyPrecipitation, STUDY_AREA_POLYGON,
      ee.Reducer.mean(), 5000, 'system:time_start')
      .setOptions({ title: 'Mean monthly CHIRPS precipitation (mm)' }));

print(ui.Chart.image.series(monthlyTemperature, STUDY_AREA_POLYGON,
      ee.Reducer.mean(), 11000, 'system:time_start')
      .setOptions({ title: 'Mean monthly ERA5-Land 2 m temperature (deg C)' }));

print(ui.Chart.image.series(monthlySnowCover, STUDY_AREA_POLYGON,
      ee.Reducer.sum(), 500, 'system:time_start')
      .setOptions({ title: 'Total monthly MOD10A1 snow-covered area (100 km^2)' }));

print(ui.Chart.image.series(fldasSnow, STUDY_AREA_POLYGON,
      ee.Reducer.sum(), 4500, 'system:time_start')
      .setOptions({ title: 'Total monthly FLDAS snow-covered area (km^2)' }));

// ------------------------------------------------------------------
// 8. Per-basin zonal statistics, exported as a long-format CSV
// ------------------------------------------------------------------
// We tag each variable image with its variable name and reduce it over every
// selected basin; the resulting rows are flattened into a single CSV that is
// easy to pivot in pandas / R.
var reduceVariableOverBasins = function (variableCollection, variableName, scaleM) {
  return variableCollection.map(function (image) {
    return image.reduceRegions({
      collection: selectedBasins.select([BASIN_ID_FIELD]),
      reducer: ee.Reducer.mean(),
      scale: scaleM
    }).map(function (feature) {
      return feature
        .set('variable', variableName)
        .set('imageDate', image.date().format('YYYY-MM-dd'))
        .setGeometry(null);
    });
  }).flatten();
};

var precipitationRows = reduceVariableOverBasins(monthlyPrecipitation, 'precip_mm', 5000);
var temperatureRows   = reduceVariableOverBasins(monthlyTemperature,   'temp_c',    11000);
var snowCoverRows     = reduceVariableOverBasins(monthlySnowCover,     'snow_km2',  500);
var fldasSnowRows     = reduceVariableOverBasins(fldasSnow,            'fldas_snow_km2', 4500);

var hydroclimateLongTable = precipitationRows
  .merge(temperatureRows)
  .merge(snowCoverRows)
  .merge(fldasSnowRows);

Export.table.toDrive({
  collection:  hydroclimateLongTable,
  description: 'river_basin_hydroclimate_long_table',
  fileFormat:  'CSV'
});
