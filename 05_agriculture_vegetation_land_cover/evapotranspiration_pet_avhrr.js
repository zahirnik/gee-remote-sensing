/*
 * Potential and Actual Evapotranspiration from MODIS and AVHRR NDVI baselines
 *
 * Datasets:    MODIS/006/MOD16A2     (8-day PET / ET, 500 m)
 *              NOAA/CDR/AVHRR/NDVI/V5 (long-term NDVI baseline, 1981-present, 0.05 deg)
 *              IDAHO_EPSCOR/TERRACLIMATE (monthly actual ET, 1958-present, 4 km)
 * Region:      Configurable AOI (default: MENA grasslands)
 * Output:      Monthly PET and ET ImageCollections, time-series chart, and a
 *              Drive export of monthly PET grids.
 *
 * Usage:
 *   1. Paste into the GEE Code Editor.
 *   2. Set STUDY_AREA and the monitoring window (START_DATE / END_DATE).
 *   3. Click Run; monthly PET / ET stacks are built and a long-term chart is shown.
 */

// -------------------------------------------------------------------------
// 1. User configuration
// -------------------------------------------------------------------------
var STUDY_AREA       = ee.Geometry.Rectangle([20, 15, 60, 40]);  // MENA grasslands
var START_DATE       = '2001-01-01';
var END_DATE         = '2020-01-01';
var AVHRR_START_DATE = '1981-07-01';
var AVHRR_END_DATE   = '2020-01-01';
var MONTHS_IN_RANGE  = 12 * 19;       // PET window: 2001-2019 = 228 months
var EXPORT_SCALE_M   = 1000;          // Resample target for export
var TARGET_CRS       = 'EPSG:4326';

// -------------------------------------------------------------------------
// 2. Generic monthly compositor (reusable across PET / ET / NDVI inputs)
//    `reducerType` controls aggregation: 'sum' for fluxes, 'mean' for indices.
// -------------------------------------------------------------------------
function monthlyComposite(collection, startDateString, monthCount, reducerType) {
  var origin   = ee.Date(startDateString);
  var sequence = ee.List.sequence(0, monthCount - 1);

  return ee.ImageCollection(sequence.map(function (i) {
    var start = origin.advance(ee.Number(i),           'month');
    var end   = origin.advance(ee.Number(i).add(1),    'month');
    var slice = collection.filterDate(start, end);
    var reduced = reducerType === 'sum' ? slice.sum() : slice.mean();
    return reduced.set({
      'system:time_start': start.millis(),
      'system:time_end':   end.millis()
    });
  }));
}

// -------------------------------------------------------------------------
// 3. MOD16A2 PET (Potential Evapotranspiration)
//    MOD16A2 reports 8-day cumulative PET in 0.1 mm; we monthly-sum the cubes
//    then rescale to mm.
// -------------------------------------------------------------------------
var modisPotentialEt = ee.ImageCollection('MODIS/006/MOD16A2')
                         .filterBounds(STUDY_AREA)
                         .filterDate(START_DATE, END_DATE)
                         .select('PET');

var monthlyPetRaw  = monthlyComposite(modisPotentialEt, START_DATE,
                                      MONTHS_IN_RANGE, 'sum');

var monthlyPetMm = monthlyPetRaw.map(function (img) {
  return img.multiply(0.1).rename('PET_mm')
            .copyProperties(img, ['system:time_start', 'system:time_end']);
});

// -------------------------------------------------------------------------
// 4. TerraClimate actual evapotranspiration (AET) - long-record back-stop
// -------------------------------------------------------------------------
var terraclimateActualEt = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')
                             .select('aet')
                             .filterDate(START_DATE, END_DATE)
                             .map(function (img) {
                               return img.clip(STUDY_AREA)
                                         .copyProperties(img,
                                           ['system:time_start', 'system:time_end']);
                             });

// -------------------------------------------------------------------------
// 5. AVHRR NDVI baseline - vegetation greenness context for ET interpretation
// -------------------------------------------------------------------------
var avhrrNdviScaled = ee.ImageCollection('NOAA/CDR/AVHRR/NDVI/V5')
                        .filterDate(AVHRR_START_DATE, AVHRR_END_DATE)
                        .select('NDVI')
                        .map(function (img) {
                          return img.clip(STUDY_AREA).multiply(0.0001)
                                    .rename('NDVI')
                                    .copyProperties(img,
                                      ['system:time_start', 'system:time_end']);
                        });

var monthsOfAvhrrRecord = 462;        // 1981-07 through 2019-12
var monthlyAvhrrNdvi = monthlyComposite(avhrrNdviScaled, AVHRR_START_DATE,
                                        monthsOfAvhrrRecord, 'mean');

// -------------------------------------------------------------------------
// 6. PET / NDVI charts over the AOI
// -------------------------------------------------------------------------
var petChart = ui.Chart.image.series({
  imageCollection: monthlyPetMm,
  region:          STUDY_AREA,
  reducer:         ee.Reducer.mean(),
  scale:           EXPORT_SCALE_M,
  xProperty:       'system:time_start'
}).setOptions({
  title: 'Monthly PET (mm) - MOD16A2',
  hAxis: {title: 'Date'},
  vAxis: {title: 'PET (mm / month)'},
  series: {0: {color: '#d95f02', lineWidth: 1.5, pointSize: 2}}
});
print(petChart);

var ndviChart = ui.Chart.image.series({
  imageCollection: monthlyAvhrrNdvi,
  region:          STUDY_AREA,
  reducer:         ee.Reducer.mean(),
  scale:           5000,
  xProperty:       'system:time_start'
}).setOptions({
  title: 'AVHRR CDR monthly NDVI - long-term greenness baseline',
  hAxis: {title: 'Date'},
  vAxis: {title: 'NDVI'},
  series: {0: {color: '#1b9e77', lineWidth: 1.2, pointSize: 1.5}}
});
print(ndviChart);

// -------------------------------------------------------------------------
// 7. Visualise the latest PET grid and AVHRR climatology mean
// -------------------------------------------------------------------------
Map.centerObject(STUDY_AREA, 4);
Map.addLayer(monthlyPetMm.first(),
             {min: 0, max: 200, palette: ['#ffffe5', '#fed98e', '#fe9929', '#cc4c02']},
             'PET - first monthly composite');
Map.addLayer(monthlyAvhrrNdvi.mean(),
             {min: -0.1, max: 0.7,
              palette: ['#a50026', '#fee08b', '#66bd63', '#006837']},
             'AVHRR NDVI climatology');

// -------------------------------------------------------------------------
// 8. Export monthly PET stack
// -------------------------------------------------------------------------
Export.image.toDrive({
  image:       monthlyPetMm.toBands(),
  description: 'monthly_PET_mm_' + START_DATE.substring(0, 4) +
               '_' + END_DATE.substring(0, 4),
  folder:      'GEE_exports',
  region:      STUDY_AREA,
  scale:       EXPORT_SCALE_M,
  crs:         TARGET_CRS,
  maxPixels:   1e13
});
