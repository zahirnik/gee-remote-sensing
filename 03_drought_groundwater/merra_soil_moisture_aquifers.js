/*
 * MERRA-2 root-zone soil moisture monitoring over aquifer polygons
 *
 * Datasets:    NASA/GSFC/MERRA/lnd/2 (MERRA-2 land surface diagnostics, hourly,
 *                                      ~55 km, 1980-present)
 *              Bands used: RZMC (root-zone soil moisture content, m3/m3)
 *                          GWETROOT optional (root-zone wetness, dimensionless)
 * Region:      Aquifer FeatureCollection (default: MENA aquifers)
 * Output:      Monthly RZMC image collection, per-aquifer CSV time series and
 *              an AOI-mean time-series chart
 *
 * Method:      Hourly MERRA-2 RZMC is averaged into calendar months over the
 *              user-defined period. Each monthly image is then reduced over
 *              every aquifer polygon to obtain a long-term soil moisture
 *              record useful for groundwater recharge and agricultural
 *              drought studies (low RZMC suppresses deep percolation).
 *
 * Usage:
 *   1. Paste into the GEE Code Editor.
 *   2. Adjust AOI, START_DATE, END_DATE and N_MONTHS at the top.
 *   3. Click Run; chart appears in the console; per-polygon CSV export task
 *      is queued in the Tasks tab.
 */

// ---------------------------------------------------------------------------
// USER CONFIGURATION
// ---------------------------------------------------------------------------

var AOI            = ee.FeatureCollection('projects/bamboo-creek-269221/assets/Aquifers_Drought');
var POLYGON_ID     = 'OBJECTID';

var START_DATE     = '2004-01-01';
var END_DATE       = '2019-01-01';
var N_MONTHS       = 180;

// MERRA-2 land grid is ~0.5 deg lon x ~0.625 deg lat (~55 km).
var SCALE_METERS   = 55000;
var EXPORT_FOLDER  = 'GEE_Soil_Moisture';

// ---------------------------------------------------------------------------
// LOAD MERRA-2 ROOT-ZONE SOIL MOISTURE
// ---------------------------------------------------------------------------

var merraLand = ee.ImageCollection('NASA/GSFC/MERRA/lnd/2')
  .filterDate(START_DATE, END_DATE);

var rootZoneMoisture = merraLand.select('RZMC');   // m3/m3

// ---------------------------------------------------------------------------
// HELPER: MONTHLY MEAN COMPOSITE
// ---------------------------------------------------------------------------

function monthlyMean (collection, startDate, monthCount) {
  var origin = ee.Date(startDate);
  var months = ee.List.sequence(0, ee.Number(monthCount).subtract(1));

  return ee.ImageCollection(months.map(function (i) {
    var monthStart = origin.advance(ee.Number(i),         'month');
    var monthEnd   = origin.advance(ee.Number(i).add(1),  'month');

    return collection.filterDate(monthStart, monthEnd).mean().set({
      'system:time_start': monthStart.millis(),
      'system:time_end':   monthEnd.millis()
    });
  }));
}

var rzmcMonthly = monthlyMean(rootZoneMoisture, START_DATE, N_MONTHS);

// ---------------------------------------------------------------------------
// MAP DISPLAY
// ---------------------------------------------------------------------------

Map.centerObject(AOI, 5);
Map.addLayer(
  rzmcMonthly.mean(),
  { min: 0, max: 0.5, palette: ['ffffd9', 'a1dab4', '41b6c4', '2c7fb8', '253494'] },
  'MERRA-2 root-zone soil moisture (period mean)'
);
Map.addLayer(AOI, { color: 'red' }, 'Aquifer polygons', false);

// ---------------------------------------------------------------------------
// AOI-MEAN TIME-SERIES CHART
// ---------------------------------------------------------------------------

var rzmcChart = ui.Chart.image.seriesByRegion({
  imageCollection: rzmcMonthly,
  regions:         AOI,
  reducer:         ee.Reducer.mean(),
  scale:           SCALE_METERS,
  xProperty:       'system:time_start',
  seriesProperty:  'RZMC'
}).setChartType('LineChart').setOptions({
  title:     'MERRA-2 root-zone soil moisture -- monthly mean over aquifers',
  hAxis:     { title: 'Date' },
  vAxis:     { title: 'RZMC [m3/m3]' },
  lineWidth: 1,
  pointSize: 2
});
print(rzmcChart);

// ---------------------------------------------------------------------------
// PER-POLYGON TRIPLET EXPORT
// ---------------------------------------------------------------------------

var triplets = rzmcMonthly.map(function (image) {
  return image.select('RZMC').reduceRegions({
    collection: AOI.select([POLYGON_ID]),
    reducer:    ee.Reducer.mean(),
    scale:      SCALE_METERS
  }).map(function (feature) {
    return feature
      .set('imageId',   image.id())
      .set('imageDate', image.date())
      .setGeometry(null);
  });
}).flatten();

Export.table.toDrive({
  collection:  triplets,
  description: 'merra2_rzmc_aquifers',
  folder:      EXPORT_FOLDER,
  fileFormat:  'CSV'
});
