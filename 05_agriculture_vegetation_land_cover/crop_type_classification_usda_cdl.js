/*
 * Crop-Type Classification with USDA NASS CDL Labels and Sentinel-2 Features
 *
 * Datasets:    USDA/NASS/CDL                  (annual crop labels, 30 m)
 *              COPERNICUS/S2_SR_HARMONIZED    (Sentinel-2 surface reflectance)
 *              MODIS/MOD09GA_006_NDVI         (auxiliary NDVI baseline)
 * Region:      Configurable AOI (default: central Illinois corn / soy belt)
 * Output:      Classified crop-type raster (Random Forest), county-level crop-area
 *              CSV export, and validation confusion matrix in the Console.
 *
 * Usage:
 *   1. Paste into the GEE Code Editor.
 *   2. Set STUDY_AREA and TARGET_YEAR; both must fall within CDL coverage (US, 1997-).
 *   3. Click Run; the classified map and county-area export task appear.
 */

// -------------------------------------------------------------------------
// 1. User configuration
// -------------------------------------------------------------------------
var STUDY_AREA            = ee.Geometry.Rectangle([-91.5, 39.5, -87.5, 42.5]); // Illinois
var TARGET_YEAR           = 2022;
var GROWING_SEASON_START  = '05-01';   // MM-DD
var GROWING_SEASON_END    = '09-30';
var CLOUD_PROB_THRESHOLD  = 20;        // QA60-derived cloud probability cutoff (%)
var NUM_TREES_RF          = 100;       // 100 trees: standard ceiling for CDL-scale problems
var TRAIN_FRACTION        = 0.7;       // 70 / 30 split is a defensible default
var SAMPLE_POINTS_PER_CLASS = 500;     // Stratified sample size per CDL crop class
var INPUT_BANDS           = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12',
                             'NDVI', 'NDWI', 'NBR'];
var CDL_CROP_CLASSES = {
  1:  'corn',
  2:  'cotton',
  4:  'sorghum',
  5:  'soybeans',
  21: 'barley',
  22: 'durum_wheat',
  23: 'spring_wheat',
  24: 'winter_wheat',
  28: 'oats',
  36: 'alfalfa'
};

// -------------------------------------------------------------------------
// 2. Load CDL crop-label raster for the target year
// -------------------------------------------------------------------------
var cdlCropImage = ee.ImageCollection('USDA/NASS/CDL')
                     .filterDate(TARGET_YEAR + '-01-01', TARGET_YEAR + '-12-31')
                     .first()
                     .select('cropland')
                     .clip(STUDY_AREA);

// Restrict to the crop classes we want to model so the classifier doesn't waste
// trees on background labels (water, urban, forest, etc.).
var targetClassValues = ee.List(Object.keys(CDL_CROP_CLASSES).map(Number));
var cdlMaskedToTargets = cdlCropImage.updateMask(
  cdlCropImage.remap(targetClassValues, targetClassValues.map(function () { return 1; })).eq(1)
);

// -------------------------------------------------------------------------
// 3. Build a Sentinel-2 growing-season median composite + indices
// -------------------------------------------------------------------------
function maskSentinel2Clouds(image) {
  // S2_SR_HARMONIZED encodes cloud / cirrus in QA60: bits 10 (clouds) and 11 (cirrus).
  var qa = image.select('QA60');
  var cloudBitMask  = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var clearMask = qa.bitwiseAnd(cloudBitMask).eq(0)
                    .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(clearMask).divide(10000);
}

function addSpectralIndices(image) {
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var ndwi = image.normalizedDifference(['B3', 'B8']).rename('NDWI');
  var nbr  = image.normalizedDifference(['B8', 'B12']).rename('NBR');
  return image.addBands([ndvi, ndwi, nbr]);
}

var sentinelGrowingSeason = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(STUDY_AREA)
  .filterDate(TARGET_YEAR + '-' + GROWING_SEASON_START,
              TARGET_YEAR + '-' + GROWING_SEASON_END)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_PROB_THRESHOLD))
  .map(maskSentinel2Clouds);

var featureStack = sentinelGrowingSeason
                     .map(addSpectralIndices)
                     .median()
                     .select(INPUT_BANDS)
                     .clip(STUDY_AREA);

// -------------------------------------------------------------------------
// 4. Stratified sampling of CDL labels over the feature stack
//    Stratified sampling avoids dominance by the most abundant class (corn / soy).
// -------------------------------------------------------------------------
var stackedForSampling = featureStack.addBands(cdlMaskedToTargets.rename('crop_class'));

var sampledFeatures = stackedForSampling.stratifiedSample({
  numPoints:     SAMPLE_POINTS_PER_CLASS,
  classBand:     'crop_class',
  region:        STUDY_AREA,
  scale:         30,
  seed:          42,
  geometries:    false,
  dropNulls:     true,
  classValues:   targetClassValues,
  classPoints:   targetClassValues.map(function () { return SAMPLE_POINTS_PER_CLASS; })
});

// -------------------------------------------------------------------------
// 5. Train / test split
// -------------------------------------------------------------------------
var sampledWithRandom = sampledFeatures.randomColumn('split', 7);
var trainingFeatureCollection   = sampledWithRandom.filter(ee.Filter.lt('split', TRAIN_FRACTION));
var validationFeatureCollection = sampledWithRandom.filter(ee.Filter.gte('split', TRAIN_FRACTION));

// -------------------------------------------------------------------------
// 6. Random Forest classifier
//    NUM_TREES_RF=100 is a sensible default for ~10 crop classes; accuracy
//    typically plateaus between 80 and 150 trees on Sentinel-2 features.
// -------------------------------------------------------------------------
var randomForestClassifier = ee.Classifier.smileRandomForest(NUM_TREES_RF)
  .train({
    features:        trainingFeatureCollection,
    classProperty:   'crop_class',
    inputProperties: INPUT_BANDS
  });

// -------------------------------------------------------------------------
// 7. Predict and validate
// -------------------------------------------------------------------------
var classifiedCropMap = featureStack.classify(randomForestClassifier);

var validated = validationFeatureCollection.classify(randomForestClassifier);
var validationErrorMatrix = validated.errorMatrix('crop_class', 'classification');
print('Validation error matrix:',  validationErrorMatrix);
print('Validation overall accuracy:', validationErrorMatrix.accuracy());
print('Validation kappa:',           validationErrorMatrix.kappa());

// -------------------------------------------------------------------------
// 8. County-level crop-area summaries (per-crop area in km^2)
// -------------------------------------------------------------------------
var usCounties = ee.FeatureCollection('TIGER/2018/Counties')
                   .filterBounds(STUDY_AREA);

var perCropBinaryStack = ee.Image.cat(
  Object.keys(CDL_CROP_CLASSES).map(function (cls) {
    var crop = CDL_CROP_CLASSES[cls];
    return cdlCropImage.eq(Number(cls)).rename(crop);
  })
);

var perCropAreaKm2 = perCropBinaryStack
                       .multiply(ee.Image.pixelArea())
                       .divide(1e6);

var countyCropAreaTable = perCropAreaKm2.reduceRegions({
  collection: usCounties.select(['NAME', 'GEOID']),
  reducer:    ee.Reducer.sum(),
  scale:      30
});

Export.table.toDrive({
  collection:  countyCropAreaTable,
  description: 'county_crop_area_km2_' + TARGET_YEAR,
  folder:      'GEE_exports',
  fileFormat:  'CSV'
});

// -------------------------------------------------------------------------
// 9. Map visualisation
// -------------------------------------------------------------------------
Map.centerObject(STUDY_AREA, 7);
Map.addLayer(featureStack, {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.3},
             'Sentinel-2 growing-season median (RGB)');
Map.addLayer(classifiedCropMap,
             {min: 1, max: 36,
              palette: ['#ffd400', '#ff2626', '#a50f15', '#267300', '#deb887',
                        '#c4a484', '#6a51a3', '#225ea8', '#9ecae1', '#a1d99b']},
             'Random Forest crop-type classification');
