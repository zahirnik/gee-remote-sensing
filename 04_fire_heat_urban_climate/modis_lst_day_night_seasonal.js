/*
 * MODIS LST day vs night seasonal patterns (Terra MOD11A1 + Aqua MYD11A1)
 *
 * Datasets:    MODIS/006/MOD11A1 (Terra daily 1 km LST, 10:30 / 22:30 local)
 *              MODIS/006/MYD11A1 (Aqua  daily 1 km LST, 13:30 / 01:30 local)
 *              An asset of station / sample points (one row per location with a
 *              'station_ID' property) for time-series extraction.
 * Region:      Any user-defined ROI polygon. Default example covers a large
 *              study area roughly over the Iranian plateau.
 * Output:      Four CSV exports to Drive:
 *                Terra day, Terra night, Aqua day, Aqua night.
 *              Each CSV is a tall table of (station_ID, imageId, LST_mean,
 *              view_time, QC) - ready to pivot in pandas / R for seasonal,
 *              diurnal and inter-platform analyses.
 *
 * Why both Terra and Aqua, both day and night?
 *   - Terra and Aqua sample at different local solar times, so combining them
 *     gives 4 daily LST samples (~01:30, 10:30, 13:30, 22:30) which lets you
 *     reconstruct the diurnal temperature cycle.
 *   - Day vs night LST differ by ~10-30 K depending on surface and season;
 *     splitting them is essential for seasonal-amplitude / SUHI work.
 *
 * Usage:
 *   1. Paste this script into the GEE Code Editor.
 *   2. Set STATION_POINTS_ASSET, ROI, START_DATE, END_DATE.
 *   3. Click Run; four tasks will appear in the Tasks tab.
 */

// ------------------------------------------------------------------
// 0. User-editable constants
// ------------------------------------------------------------------
var STATION_POINTS_ASSET = 'users/<YOUR_GEE_USERNAME>/LST_POINTS';   // points with 'station_ID' property
var ROI = ee.Geometry.Polygon(
  [[[40.20253837333929, 41.8367090611045],
    [40.20253837333929, 22.715242934985714],
    [68.1517571233393,  22.715242934985714],
    [68.1517571233393,  41.8367090611045]]], null, false
);

var START_DATE = '2000-02-24';   // Terra MODIS launched late 1999; data begin Feb 2000
var END_DATE   = '2021-12-31';

var STATION_ID_PROPERTY  = 'station_ID';
var REDUCTION_SCALE_M    = 1000;             // MOD11A1 native resolution
var LST_SCALE_FACTOR     = 0.02;             // DN -> Kelvin

// ------------------------------------------------------------------
// 1. Load station points
// ------------------------------------------------------------------
var stationPoints = ee.FeatureCollection(STATION_POINTS_ASSET);

// ------------------------------------------------------------------
// 2. Helpers
// ------------------------------------------------------------------
// Build a per-station tall table from a MODIS LST collection for either the
// 'LST_Day_1km' or 'LST_Night_1km' band, with its companion view-time and QC bands.
function extractStationTimeSeries(collection, lstBand, viewTimeBand, qcBand) {
  return collection.map(function (image) {
    var lstKelvin = image.select(lstBand).multiply(LST_SCALE_FACTOR).rename('LST_K');
    var enriched  = lstKelvin
      .addBands(image.select(viewTimeBand))
      .addBands(image.select(qcBand));

    return enriched.reduceRegions({
      collection: stationPoints.select([STATION_ID_PROPERTY]),
      reducer:    ee.Reducer.mean(),
      scale:      REDUCTION_SCALE_M
    }).map(function (feature) {
      return feature
        .set('imageId',    image.id())
        .set('time_start', image.date().format('yyyy-MM-dd'));
    });
  }).flatten();
}

function buildLstCollection(collectionId, lstBand, viewTimeBand, qcBand) {
  return ee.ImageCollection(collectionId)
    .filterBounds(ROI)
    .filterDate(START_DATE, END_DATE)
    .select([lstBand, viewTimeBand, qcBand]);
}

// ------------------------------------------------------------------
// 3. Terra (MOD11A1) - day and night
// ------------------------------------------------------------------
var terraDayCollection = buildLstCollection(
  'MODIS/006/MOD11A1', 'LST_Day_1km',   'Day_view_time',   'QC_Day'
);
var terraNightCollection = buildLstCollection(
  'MODIS/006/MOD11A1', 'LST_Night_1km', 'Night_view_time', 'QC_Night'
);

var terraDayTimeSeries = extractStationTimeSeries(
  terraDayCollection,   'LST_Day_1km',   'Day_view_time',   'QC_Day'
);
var terraNightTimeSeries = extractStationTimeSeries(
  terraNightCollection, 'LST_Night_1km', 'Night_view_time', 'QC_Night'
);

// ------------------------------------------------------------------
// 4. Aqua (MYD11A1) - day and night
// ------------------------------------------------------------------
var aquaDayCollection = buildLstCollection(
  'MODIS/006/MYD11A1', 'LST_Day_1km',   'Day_view_time',   'QC_Day'
);
var aquaNightCollection = buildLstCollection(
  'MODIS/006/MYD11A1', 'LST_Night_1km', 'Night_view_time', 'QC_Night'
);

var aquaDayTimeSeries = extractStationTimeSeries(
  aquaDayCollection,   'LST_Day_1km',   'Day_view_time',   'QC_Day'
);
var aquaNightTimeSeries = extractStationTimeSeries(
  aquaNightCollection, 'LST_Night_1km', 'Night_view_time', 'QC_Night'
);

// ------------------------------------------------------------------
// 5. Map preview - example day vs night layers from the first acquisition
// ------------------------------------------------------------------
Map.centerObject(ROI, 5);
Map.addLayer(stationPoints, {color: 'black'}, 'Station points');

var lstDayVis   = {min: 280, max: 320, palette: ['blue', 'white', 'red']};
var lstNightVis = {min: 260, max: 300, palette: ['blue', 'white', 'red']};

Map.addLayer(
  terraDayCollection.first().select('LST_Day_1km').multiply(LST_SCALE_FACTOR),
  lstDayVis, 'Terra day LST (K) - first scene'
);
Map.addLayer(
  terraNightCollection.first().select('LST_Night_1km').multiply(LST_SCALE_FACTOR),
  lstNightVis, 'Terra night LST (K) - first scene'
);

// ------------------------------------------------------------------
// 6. Drive exports - one CSV per platform-and-overpass combination
// ------------------------------------------------------------------
Export.table.toDrive({
  collection: terraDayTimeSeries,
  description: 'Export_TerraDay_LST',
  folder: 'Terra_Day',
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: terraNightTimeSeries,
  description: 'Export_TerraNight_LST',
  folder: 'Terra_Night',
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: aquaDayTimeSeries,
  description: 'Export_AquaDay_LST',
  folder: 'Aqua_Day',
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: aquaNightTimeSeries,
  description: 'Export_AquaNight_LST',
  folder: 'Aqua_Night',
  fileFormat: 'CSV'
});
