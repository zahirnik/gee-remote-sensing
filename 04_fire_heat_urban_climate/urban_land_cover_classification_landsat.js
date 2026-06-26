/*
 * Supervised Landsat land-cover classification for urban areas (Isfahan, Iran)
 *
 * Datasets:    LANDSAT/LT05/C01/T1_TOA, LANDSAT/LE07/C01/T1_TOA,
 *              LANDSAT/LC08/C01/T1_TOA (Top-of-Atmosphere, single-date scenes)
 *              USGS/SRTMGL1_003 (30 m elevation)
 *              Landsat SMW LST module (Ermida et al. 2020):
 *              users/sofiaermida/landsat_smw_lst:modules/Landsat_LST.js
 *              Training / test polygons as user assets per epoch
 * Region:      Isfahan metropolitan area (Iran). Pipeline is generic and runs
 *              for any city with matching training assets.
 * Output:      (1) Per-epoch Landsat land-cover classification raster.
 *              (2) Validation error matrix, overall / kappa / producers /
 *                  consumers accuracy.
 *              (3) Variable-importance chart (random-forest feature ranking).
 *
 * Feature stack used for classification:
 *   - Reflectance bands (Landsat 5/7/8 TOA, renamed by sensor module)
 *   - 21 spectral indices (NDVI, NDBI, NDWI, MNDWI, NBR, SAVI, EVI, IBI, ...)
 *   - SRTM elevation, slope, aspect
 *   - Landsat-derived Land Surface Temperature (Ermida SMW)
 *
 * Usage:
 *   1. Paste this script into the GEE Code Editor.
 *   2. Replace LANDSAT_SCENE_ID, STUDY_AREA_ASSET, TRAINING_ASSET,
 *      TEST_ASSET and matching LST inputs for the epoch you want to classify.
 *   3. Click Run; classification layer appears on the map, validation
 *      metrics and importance chart print to the console.
 */

// ------------------------------------------------------------------
// 0. User-editable constants - swap this block per epoch (1985, 1993, ...)
// ------------------------------------------------------------------
var STUDY_AREA_ASSET  = 'users/<your_account>/IsfahanLCC/Border_GEE';
var TRAINING_ASSET    = 'users/<your_account>/IsfahanLCC/1985_Train';
var TEST_ASSET        = 'users/<your_account>/IsfahanLCC/1985_Test';

var LANDSAT_SCENE_ID  = 'LANDSAT/LT05/C01/T1_TOA/LT05_164037_19850802';
var LST_SATELLITE     = 'L5';            // 'L5' / 'L7' / 'L8' for the Ermida SMW module
var LST_START_DATE    = '1985-08-01';
var LST_END_DATE      = '1985-08-03';
var LST_USE_NDVI      = true;

var CLASS_PROPERTY    = 'LC_Class';      // attribute on training polygons (integer class label)
var SAMPLE_POINTS     = 100000;          // total stratified-sample budget
var TRAIN_TEST_SPLIT  = 0.3;             // test-budget fraction relative to SAMPLE_POINTS

var RANDOM_FOREST_TREES = 50;
var REDUCTION_SCALE_M   = 30;            // use 15 m for L7/L8 TOA pan-sharpened products
var OUTPUT_PROJECTION   = 'EPSG:32639';  // UTM zone 39N - swap to match your AOI
var OUTPUT_ASSET_NAME   = 'classified_1985';
var CHART_TITLE         = 'Variable importance: Landsat 5 1985/08/02';

var CLASSIFICATION_VIS = {
  min: 1, max: 8,
  palette: ['956300', 'dcef63', 'c31300', 'fff5d6',
            '0046c7', '00785a', '009900', '008b00']
};

// ------------------------------------------------------------------
// 1. Load study area and training / test polygons
// ------------------------------------------------------------------
var studyArea       = ee.FeatureCollection(STUDY_AREA_ASSET);
var trainingPolygons = ee.FeatureCollection(TRAINING_ASSET);
var testingPolygons  = ee.FeatureCollection(TEST_ASSET);

Map.centerObject(studyArea);

// ------------------------------------------------------------------
// 2. Landsat Surface Temperature via the SMW community module
// ------------------------------------------------------------------
// The single-channel SMW (Statistical Mono-Window) algorithm of
// Ermida et al. (2020) wraps emissivity, atmospheric correction and
// the LST inversion behind a clean collection() helper.
var landsatLstModule = require('users/sofiaermida/landsat_smw_lst:modules/Landsat_LST.js');

var landsatLstCollection = landsatLstModule.collection(
  LST_SATELLITE, LST_START_DATE, LST_END_DATE, studyArea, LST_USE_NDVI
);

var lstBand = ee.Image(landsatLstCollection.first())
  .select('LST')
  .clip(studyArea)
  .rename('LST')
  .reproject(OUTPUT_PROJECTION, null, REDUCTION_SCALE_M);

// ------------------------------------------------------------------
// 3. Landsat indices module (NDVI, NDBI, NDWI, NBR, IBI, EBBI, ...)
// ------------------------------------------------------------------
// Local repo module that returns the scene's reflectance bands plus
// ~21 spectral indices in a single multi-band image.
var landsatIndicesModule = require('users/<your_account>/LCC_Isfahan:LANDSAT_INDICES');
var indicesCollection = landsatIndicesModule.collection(LANDSAT_SCENE_ID, studyArea);

// ------------------------------------------------------------------
// 4. Terrain features from SRTM (elevation, slope, aspect)
// ------------------------------------------------------------------
// Topography drives micro-climate and tends to be a top-ranked predictor
// when the AOI spans elevation gradients (Isfahan: ~1500-3500 m).
var srtm = ee.Image('USGS/SRTMGL1_003');
var elevation = srtm.select('elevation')
  .clip(studyArea)
  .rename('DEM')
  .reproject(OUTPUT_PROJECTION, null, REDUCTION_SCALE_M);
var aspect = ee.Terrain.aspect(elevation).rename('ASPECT');
var slope  = ee.Terrain.slope(elevation).rename('SLOPE');

// ------------------------------------------------------------------
// 5. Assemble the full feature stack
// ------------------------------------------------------------------
var classificationInputs = indicesCollection
  .addBands(elevation)
  .addBands(aspect)
  .addBands(slope)
  .addBands(lstBand);

// ------------------------------------------------------------------
// 6. Stratified sampling
// ------------------------------------------------------------------
// We rasterize the polygon class labels into an extra band so that
// stratifiedSample() can draw samples *per class* in proportion to each
// class's polygon area - this avoids over-sampling large rare classes.
function buildStratifiedSample(polygons, sampleBudget) {
  var labelBand = ee.Image().int()
    .paint(polygons, CLASS_PROPERTY)
    .rename(CLASS_PROPERTY)
    .toInt();

  var stackForSampling = classificationInputs.addBands(labelBand);
  var classValues = polygons.aggregate_array(CLASS_PROPERTY).distinct();
  var totalArea = polygons.geometry().area();

  // Allocate sample count per class proportional to that class's polygon area.
  var classSampleCounts = classValues.map(function (classValue) {
    return polygons.filter(ee.Filter.equals(CLASS_PROPERTY, classValue))
      .geometry().area()
      .divide(totalArea)
      .multiply(sampleBudget)
      .ceil();
  });

  return stackForSampling.stratifiedSample({
    numPoints:   sampleBudget,
    classBand:   CLASS_PROPERTY,
    region:      polygons,
    scale:       REDUCTION_SCALE_M,
    dropNulls:   true,
    geometries:  true,
    classValues: classValues,
    classPoints: classSampleCounts
  });
}

var trainingSamples = buildStratifiedSample(trainingPolygons, SAMPLE_POINTS);
var testingSamples  = buildStratifiedSample(testingPolygons, TRAIN_TEST_SPLIT * SAMPLE_POINTS);

// ------------------------------------------------------------------
// 7. Train Random Forest classifier
// ------------------------------------------------------------------
var randomForestClassifier = ee.Classifier.smileRandomForest(RANDOM_FOREST_TREES).train({
  features:      trainingSamples,
  classProperty: CLASS_PROPERTY
});

var classifiedLandCover = classificationInputs.classify(randomForestClassifier);

// Out-of-bag error: how the trees would have scored on samples they did
// not see during training. Cheap proxy for generalisation error.
var classifierExplanation = randomForestClassifier.explain();
var outOfBagError = ee.Dictionary(classifierExplanation).get('outOfBagErrorEstimate');
print('Out-of-bag error estimate:', outOfBagError);

Map.addLayer(classifiedLandCover, CLASSIFICATION_VIS, 'Classified land cover');

// ------------------------------------------------------------------
// 8. Validation - error matrix and standard accuracy metrics
// ------------------------------------------------------------------
var validated = testingSamples.classify(randomForestClassifier);
var errorMatrix = validated.errorMatrix(CLASS_PROPERTY, 'classification');

print('Validation error matrix:',     errorMatrix);
print('Validation overall accuracy:', errorMatrix.accuracy());
print('Validation kappa accuracy:',   errorMatrix.kappa());
print('Validation consumers accuracy:', errorMatrix.consumersAccuracy());
print('Validation producers accuracy:', errorMatrix.producersAccuracy());

// ------------------------------------------------------------------
// 9. Variable-importance chart (random-forest feature ranking)
// ------------------------------------------------------------------
var importance = ee.Feature(null, ee.Dictionary(classifierExplanation).get('importance'));
print('Raw variable importance:', importance);

function sortDictionaryByValue(rawDictionary) {
  var sortable = [];
  for (var key in rawDictionary) {
    if (rawDictionary.hasOwnProperty(key)) {
      sortable.push([key, rawDictionary[key]]);
    }
  }
  sortable.sort(function (a, b) { return a[1] - b[1]; });

  var sortedKeys = [];
  var sortedValues = [];
  for (var i = 0; i < sortable.length; i += 1) {
    sortedKeys.push(sortable[i][0]);
    sortedValues.push(sortable[i][1]);
  }
  return [sortedKeys, sortedValues];
}

var sortedImportance = sortDictionaryByValue(importance.toDictionary().getInfo());
var importanceChart = ui.Chart.array.values(
    sortedImportance[1].reverse(), 0, sortedImportance[0].reverse()
  )
  .setChartType('ColumnChart')
  .setOptions({
    title: CHART_TITLE,
    legend: {position: 'none'},
    hAxis:  {title: 'Predictor', textStyle: {fontSize: 9}},
    vAxis:  {title: 'Importance'}
  });
print(importanceChart);

// ------------------------------------------------------------------
// 10. Optional Drive / Asset export of the classified raster
// ------------------------------------------------------------------
// Uncomment to push the classified raster to the Assets store.
//
// Export.image.toAsset({
//   image:       classifiedLandCover,
//   description: OUTPUT_ASSET_NAME,
//   assetId:     OUTPUT_ASSET_NAME,
//   region:      studyArea.geometry(),
//   scale:       REDUCTION_SCALE_M
// });
