/*
 * Multi-lake water-balance components: precipitation (CHIRPS), ET (MOD16) and
 * soil moisture (SMAP) aggregated per hydrological year over many lake basins.
 *
 * Datasets:    UCSB-CHG/CHIRPS/DAILY   (precipitation, mm/day, 5 km)
 *              MODIS/NTSG/MOD16A2/105  (8-day actual ET, kg/m^2/8day = mm/8day)
 *              NASA/SMAP/SPL4SMGP/007  (3-hourly surface + root-zone soil moisture)
 *              JAXA/ALOS/AW3D30/V3_2   (30 m DSM used to mask the lake water body
 *                                       out of basin-scale aggregates)
 * Region:      User-supplied lake/basin FeatureCollection with per-feature
 *              attributes 'lakeId' and 'lakeElevationM'.
 * Output:      CSV with one row per basin and one column per year for each of
 *              three water-balance components (precipitation, ET, soil moisture).
 *
 * Usage:
 *   1. Paste this script into the GEE Code Editor.
 *   2. Set BASIN_COLLECTION, BASIN_ID_PROPERTY and BASIN_ELEVATION_PROPERTY.
 *   3. Adjust the hydrological-year window if you do not use Aug-1 -> Aug-1.
 *   4. Click Run; three CSV export tasks will appear in the Tasks panel.
 */

// ------------------------------------------------------------------
// 0. User-editable constants
// ------------------------------------------------------------------
var BASIN_COLLECTION         = 'users/<YOUR_USERNAME>/<YOUR_BASIN_POLYGONS>';
var BASIN_ID_PROPERTY        = 'Hylak_id';      // unique id per basin
var BASIN_ELEVATION_PROPERTY = 'LakeEleve';     // lake water-surface elevation (m)
var LAKE_ELEV_BUFFER_M       = 5;               // mask out pixels within +/-5 m
                                                // of the lake water surface so the
                                                // open-water area does not bias
                                                // basin precipitation / ET / soil
                                                // moisture means
var START_YEAR               = 2001;
var END_YEAR                 = 2020;            // inclusive
// Hydrological year boundaries. Aug-1 -> next Aug-1 follows common practice
// for arid endorheic basins (snowpack mostly drained by late summer).
var HYDRO_YEAR_MONTH         = 8;
var HYDRO_YEAR_DAY           = 1;

// ------------------------------------------------------------------
// 1. Load basin polygons and a DSM used to define the lake-water mask
// ------------------------------------------------------------------
var basinCollection = ee.FeatureCollection(BASIN_COLLECTION);

var dsmElevation = ee.ImageCollection('JAXA/ALOS/AW3D30/V3_2')
  .select('DSM')
  .mosaic();

// ------------------------------------------------------------------
// 2. Helper to aggregate a daily/sub-daily collection into annual images
// ------------------------------------------------------------------
// reduceFn must be one of 'sum' (precipitation, ET) or 'mean' (soil moisture).
var buildAnnualCollection = function (sourceCollection, reduceFn) {
  var years = ee.List.sequence(START_YEAR, END_YEAR);
  return ee.ImageCollection.fromImages(years.map(function (year) {
    var startDate = ee.Date.fromYMD(year, HYDRO_YEAR_MONTH, HYDRO_YEAR_DAY);
    var endDate = ee.Date.fromYMD(ee.Number(year).add(1), HYDRO_YEAR_MONTH, HYDRO_YEAR_DAY);
    var yearlyImage = (reduceFn === 'sum')
      ? sourceCollection.filterDate(startDate, endDate).sum()
      : sourceCollection.filterDate(startDate, endDate).mean();
    return yearlyImage
      .set('year', year)
      .set('system:time_start', startDate.millis());
  }));
};

// ------------------------------------------------------------------
// 3. Build the three component collections
// ------------------------------------------------------------------
// CHIRPS daily precipitation (mm/day) summed -> annual precipitation (mm)
var dailyChirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY').select('precipitation');
var annualPrecipitation = buildAnnualCollection(dailyChirps, 'sum');

// MOD16 8-day ET in 0.1 mm/8day -> rescale to mm and sum to annual ET (mm/yr).
var modisEt8Day = ee.ImageCollection('MODIS/NTSG/MOD16A2/105')
  .select('ET')
  .map(function (image) {
    return image.multiply(0.1)
      .copyProperties(image, ['system:time_start', 'system:time_end']);
  });
var annualEvapotranspiration = buildAnnualCollection(modisEt8Day, 'sum');

// SMAP root-zone soil moisture, 3-hourly (m^3/m^3) -> annual mean.
var smapSoilMoisture = ee.ImageCollection('NASA/SMAP/SPL4SMGP/007')
  .select('sm_rootzone');
var annualSoilMoisture = buildAnnualCollection(smapSoilMoisture, 'mean');

// ------------------------------------------------------------------
// 4. Stack annual images into a single multi-band image per component
// ------------------------------------------------------------------
// Doing one reduceRegions() per component (instead of one per year) keeps
// the EE compute graph short and avoids hitting the 'too many concurrent
// aggregations' error on collections with thousands of basins.
var stackAnnualToBands = function (annualCollection) {
  var years = ee.List.sequence(START_YEAR, END_YEAR);
  var imageList = annualCollection.toList(annualCollection.size());
  return ee.Image(years.iterate(function (year, accumulator) {
    var index = ee.Number(year).subtract(START_YEAR);
    var bandName = ee.Number(year).format('%d');
    return ee.Image(accumulator)
      .addBands(ee.Image(imageList.get(index)).rename(bandName));
  }, ee.Image().select()));
};

var precipitationStack    = stackAnnualToBands(annualPrecipitation);
var evapotranspirationStack = stackAnnualToBands(annualEvapotranspiration);
var soilMoistureStack     = stackAnnualToBands(annualSoilMoisture);

// ------------------------------------------------------------------
// 5. Build a per-basin extractor that excludes the lake water surface
// ------------------------------------------------------------------
// We keep only pixels whose elevation is > (lakeElevation - bufferM); below
// that threshold we are on the lake itself, which we do not want to include
// in basin-scale water-balance forcings (precip falling on the lake is part
// of direct flux, not catchment input; ET over the lake is open-water ET).
var makePerBasinAggregator = function (stackedImage, reducer, aggregator, factorScale) {
  return function (basinFeature) {
    var basinElevation = ee.Number(basinFeature.get(BASIN_ELEVATION_PROPERTY));
    var elevationClipped = dsmElevation.clip(basinFeature);
    var landMask = elevationClipped.gt(basinElevation.subtract(LAKE_ELEV_BUFFER_M));

    // reduceToVectors converts the contiguous land-mask region back to a
    // FeatureCollection that reduceRegions can iterate over efficiently.
    var landPolygons = landMask.updateMask(landMask).reduceToVectors({
      reducer: ee.Reducer.countEvery(),
      geometry: basinFeature.geometry(),
      scale: 500,
      maxPixels: 1e12
    });

    var perYearStats = stackedImage.reduceRegions({
      collection: landPolygons,
      reducer: reducer,
      scale: factorScale
    });

    var years = ee.List.sequence(START_YEAR, END_YEAR);
    var output = ee.Feature(null);
    output = ee.Feature(years.iterate(function (year, accumulator) {
      var bandName = ee.Number(year).format('%d');
      var value = perYearStats[aggregator](bandName);
      return ee.Feature(accumulator).set(bandName, value);
    }, output));

    return output.copyProperties(
      basinFeature,
      [BASIN_ID_PROPERTY, BASIN_ELEVATION_PROPERTY]
    );
  };
};

// ------------------------------------------------------------------
// 6. Apply the aggregator to each component and export to Drive
// ------------------------------------------------------------------
// Precipitation: cumulative mm over the catchment -> use sum reducer.
var precipitationByBasin = basinCollection.map(makePerBasinAggregator(
  precipitationStack, ee.Reducer.mean(), 'aggregate_mean',
  dailyChirps.first().projection().nominalScale()
));

// ET: cumulative mm over the catchment land area -> mean over basin pixels.
var etByBasin = basinCollection.map(makePerBasinAggregator(
  evapotranspirationStack, ee.Reducer.mean(), 'aggregate_mean',
  modisEt8Day.first().projection().nominalScale()
));

// Soil moisture: dimensionless m^3/m^3 -> mean over basin pixels.
var soilMoistureByBasin = basinCollection.map(makePerBasinAggregator(
  soilMoistureStack, ee.Reducer.mean(), 'aggregate_mean',
  smapSoilMoisture.first().projection().nominalScale()
));

Export.table.toDrive({
  collection: precipitationByBasin,
  description: 'water_balance_precipitation_chirps',
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: etByBasin,
  description: 'water_balance_et_mod16',
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: soilMoistureByBasin,
  description: 'water_balance_soil_moisture_smap',
  fileFormat: 'CSV'
});
