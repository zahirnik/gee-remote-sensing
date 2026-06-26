/*
 * MODIS ET / PET drought monitoring and Evapotranspiration Stress Index (ESI)
 *
 * Datasets:    MODIS/006/MOD16A2 (8-day actual ET and potential ET, 500 m, 2001-present)
 * Region:      Any basin / aquifer FeatureCollection (default: Iran baseflow basins)
 * Output:      Monthly ET, PET and ESI image collections + per-polygon CSV
 *              (one row per polygon per month) for downstream drought analysis
 *
 * Method:      ESI = ET / PET. Values below the long-term mean indicate that
 *              actual evapotranspiration is falling short of atmospheric
 *              demand -- a fingerprint of vegetation water stress and the
 *              early signal of agricultural drought (Anderson et al. 2011).
 *              MOD16A2 is published with a 0.1 scale factor that is removed
 *              before any arithmetic.
 *
 * Usage:
 *   1. Paste into the GEE Code Editor.
 *   2. Adjust AOI, START_DATE, END_DATE and N_MONTHS at the top.
 *   3. Click Run; the per-polygon CSV export task is queued in the Tasks tab.
 */

// ---------------------------------------------------------------------------
// USER CONFIGURATION
// ---------------------------------------------------------------------------

var AOI            = ee.FeatureCollection('users/parizi555/IranBaseflow/400Baseflows');
var POLYGON_ID     = 'OBJECTID';

var START_DATE     = '2001-04-01';
var END_DATE       = '2021-01-01';
var N_MONTHS       = 240;             // number of monthly composites to build

var SCALE_METERS   = 500;
var MOD16_SCALE    = 0.1;             // MOD16A2 stored scale factor
var EXPORT_FOLDER  = 'GEE_ET_Drought';

// ---------------------------------------------------------------------------
// LOAD MOD16A2 ACTUAL ET AND POTENTIAL ET
// ---------------------------------------------------------------------------

var mod16 = ee.ImageCollection('MODIS/006/MOD16A2')
  .filterDate(START_DATE, END_DATE);

var etRaw  = mod16.select('ET');
var petRaw = mod16.select('PET');

// ---------------------------------------------------------------------------
// HELPER: MONTHLY COMPOSITE BUILDER
// ---------------------------------------------------------------------------
// MOD16A2 ships at 8-day cadence; we re-aggregate to calendar months by
// summing the 8-day totals that fall inside each month.

function monthlyComposite (collection, startDate, monthCount) {
  var origin = ee.Date(startDate);
  var months = ee.List.sequence(0, ee.Number(monthCount).subtract(1));

  return ee.ImageCollection(months.map(function (i) {
    var monthStart = origin.advance(ee.Number(i),         'month');
    var monthEnd   = origin.advance(ee.Number(i).add(1),  'month');

    return collection.filterDate(monthStart, monthEnd).sum().set({
      'system:time_start': monthStart.millis(),
      'system:time_end':   monthEnd.millis()
    });
  }));
}

// ---------------------------------------------------------------------------
// BUILD MONTHLY ET, PET AND ESI COLLECTIONS
// ---------------------------------------------------------------------------
// Empty months produce zero-band images; we drop those, then apply the
// scale factor and compute ESI = ET / PET.

function dropEmpty (collection) {
  return collection
    .map(function (image) {
      return image.set('bandCount', image.bandNames().length());
    })
    .filter(ee.Filter.eq('bandCount', 1));
}

var etMonthly = dropEmpty(monthlyComposite(etRaw, START_DATE, N_MONTHS))
  .map(function (img) {
    return img.multiply(MOD16_SCALE).rename('ET')
      .copyProperties(img, ['system:time_start', 'system:time_end']);
  });

var petMonthly = dropEmpty(monthlyComposite(petRaw, START_DATE, N_MONTHS))
  .map(function (img) {
    return img.multiply(MOD16_SCALE).rename('PET')
      .copyProperties(img, ['system:time_start', 'system:time_end']);
  });

// Inner join on system:time_start so we only divide ET by the matching PET.
var timeFilter = ee.Filter.equals({
  leftField:  'system:time_start',
  rightField: 'system:time_start'
});
var etPetJoined = ee.ImageCollection(
  ee.Join.inner().apply(etMonthly, petMonthly, timeFilter)
);

var esiCollection = etPetJoined.map(function (pair) {
  var etImage  = ee.Image(pair.get('primary'));
  var petImage = ee.Image(pair.get('secondary'));

  // Guard against PET == 0 (rare in deserts at night).
  var esi = etImage.divide(petImage.max(0.001)).rename('ESI');

  return etImage
    .addBands(petImage)
    .addBands(esi)
    .copyProperties(etImage, ['system:time_start', 'system:time_end']);
});

// ---------------------------------------------------------------------------
// MAP DISPLAY: LONG-TERM MEAN ESI AND THE MOST RECENT MONTH
// ---------------------------------------------------------------------------

var esiClimatology  = esiCollection.select('ESI').mean();
var esiMostRecent   = esiCollection.limit(1, 'system:time_start', false).first().select('ESI');

Map.centerObject(AOI, 6);
Map.addLayer(esiClimatology, { min: 0, max: 1, palette: ['brown', 'yellow', 'green'] },
             'ESI climatology (mean)');
Map.addLayer(esiMostRecent,  { min: 0, max: 1, palette: ['brown', 'yellow', 'green'] },
             'ESI most recent month');

// ---------------------------------------------------------------------------
// TIME-SERIES CHART (AOI MEAN)
// ---------------------------------------------------------------------------

var esiChart = ui.Chart.image.seriesByRegion({
  imageCollection: esiCollection.select(['ET', 'PET', 'ESI']),
  regions:         AOI,
  reducer:         ee.Reducer.mean(),
  scale:           SCALE_METERS,
  xProperty:       'system:time_start'
}).setOptions({
  title:  'Monthly ET, PET and ESI -- AOI mean',
  hAxis:  { title: 'Date' },
  vAxis:  { title: 'ET / PET [mm]   |   ESI [-]' }
});
print(esiChart);

// ---------------------------------------------------------------------------
// PER-POLYGON TRIPLET EXPORT
// ---------------------------------------------------------------------------

var triplets = esiCollection.map(function (image) {
  return image.reduceRegions({
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
  description: 'modis_et_pet_esi_polygons',
  folder:      EXPORT_FOLDER,
  fileFormat:  'CSV'
});
