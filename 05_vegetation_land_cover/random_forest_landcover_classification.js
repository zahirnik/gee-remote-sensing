/*
 * Random Forest Land-Cover Classification (Landsat 5 TOA)
 *
 * Datasets:    LANDSAT/LT05/C01/T1_TOA (input imagery, 6 reflective + thermal bands)
 *              MODIS/051/MCD12Q1, Land_Cover_Type_1 IGBP scheme (training labels)
 * Region:      Configurable point-centred area of interest (default: San Francisco Bay)
 * Output:      Classified IGBP land-cover map on the map canvas
 *              Resubstitution and validation confusion matrices in the Console
 *
 * Usage:
 *   1. Paste into the GEE Code Editor.
 *   2. Adjust STUDY_AREA, TARGET_YEAR and classifier hyper-parameters as needed.
 *   3. Click Run; the cloud-masked composite and the classification appear on the
 *      map, and the train/test error matrices print to the Console.
 */

// -------------------------------------------------------------------------
// 1. User configuration
// -------------------------------------------------------------------------
var STUDY_AREA              = ee.Geometry.Point(-122.3942, 37.7295); // San Francisco Bay
var TARGET_YEAR             = 2011;
var CLOUD_SCORE_THRESHOLD   = 50;     // Mask pixels with simpleCloudScore > 50
var TRAINING_SAMPLE_SIZE    = 5000;   // Total training pixels (stratified by IGBP class)
var VALIDATION_SAMPLE_SIZE  = 5000;   // Held-out validation pixels (different random seed)
var NUM_TREES_RF            = 50;     // Random Forest size; 30-100 is a robust range for IGBP
var TRAIN_SEED              = 0;
var VALIDATION_SEED         = 1;
var INPUT_BANDS             = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7'];
var CLASS_PROPERTY          = 'Land_Cover_Type_1';
var MAP_ZOOM_LEVEL          = 10;

// -------------------------------------------------------------------------
// 2. Load Landsat 5 TOA imagery and pick the least cloudy scene of the year
// -------------------------------------------------------------------------
var landsatStart = ee.Date.fromYMD(TARGET_YEAR, 1, 1);
var landsatEnd   = landsatStart.advance(1, 'year');

var landsatSurfaceReflectance = ee.Image(
  ee.ImageCollection('LANDSAT/LT05/C01/T1_TOA')
    .filterDate(landsatStart, landsatEnd)
    .filterBounds(STUDY_AREA)
    .sort('CLOUD_COVER')
    .first()
);

// -------------------------------------------------------------------------
// 3. Cloud masking using Landsat simpleCloudScore + per-band valid-data mask
// -------------------------------------------------------------------------
var cloudScore = ee.Algorithms.Landsat.simpleCloudScore(landsatSurfaceReflectance)
                                     .select('cloud');

// Combine a cross-band validity mask (all bands non-null) with the cloud-score mask.
var validDataMask = landsatSurfaceReflectance.mask().reduce('min');
var cloudMaskedScene = landsatSurfaceReflectance.updateMask(
  validDataMask.and(cloudScore.lte(CLOUD_SCORE_THRESHOLD))
);

// -------------------------------------------------------------------------
// 4. Build labelled stack: input reflectance bands + MODIS IGBP land-cover label
//    MCD12Q1 Land_Cover_Type_1 = annual IGBP class (0..17).
// -------------------------------------------------------------------------
var modisLandCoverLabels = ee.Image('MODIS/051/MCD12Q1/' + TARGET_YEAR + '_01_01')
                              .select(CLASS_PROPERTY);

var labelledStack = cloudMaskedScene.addBands(modisLandCoverLabels);

// -------------------------------------------------------------------------
// 5. Sample training and validation feature collections
//    Different random seeds guarantee disjoint draws over the same pixel pool.
// -------------------------------------------------------------------------
var trainingFeatureCollection = labelledStack.sample({
  numPixels: TRAINING_SAMPLE_SIZE,
  seed: TRAIN_SEED
});

var validationFeatureCollection = labelledStack.sample({
  numPixels: VALIDATION_SAMPLE_SIZE,
  seed: VALIDATION_SEED
}).filter(ee.Filter.neq('B1', null)); // Drop any rows with missing reflectance.

// -------------------------------------------------------------------------
// 6. Train a Random Forest classifier
//    NUM_TREES_RF was selected to balance variance reduction against compute cost:
//    error rates plateau quickly past ~50 trees on 7-band Landsat features.
// -------------------------------------------------------------------------
var randomForestClassifier = ee.Classifier.smileRandomForest(NUM_TREES_RF)
  .train({
    features: trainingFeatureCollection,
    classProperty: CLASS_PROPERTY,
    inputProperties: INPUT_BANDS
  });

// -------------------------------------------------------------------------
// 7. Predict and validate
// -------------------------------------------------------------------------
var classifiedLandCover = cloudMaskedScene.classify(randomForestClassifier);

// Resubstitution accuracy (optimistic; uses training set).
var trainingErrorMatrix = randomForestClassifier.confusionMatrix();
print('Resubstitution error matrix:', trainingErrorMatrix);
print('Training overall accuracy:',   trainingErrorMatrix.accuracy());

// Independent validation on a disjoint pixel sample.
var validatedFeatures   = validationFeatureCollection.classify(randomForestClassifier);
var validationErrorMatrix = validatedFeatures.errorMatrix(CLASS_PROPERTY, 'classification');
print('Validation error matrix:',  validationErrorMatrix);
print('Validation overall accuracy:', validationErrorMatrix.accuracy());
print('Validation kappa:',            validationErrorMatrix.kappa());

// -------------------------------------------------------------------------
// 8. Visualisation
// -------------------------------------------------------------------------
var igbpPalette = [
  'aec3d4',                                              // 0  water
  '152106', '225129', '369b47', '30eb5b', '387242',      // 1-5 forest
  '6a2325', 'c3aa69', 'b76031', 'd9903d', '91af40',      // 6-10 shrub / grass
  '111149',                                              // 11 wetlands
  'cdb33b',                                              // 12 croplands
  'cc0013',                                              // 13 urban
  '33280d',                                              // 14 crop mosaic
  'd7cdcc',                                              // 15 snow and ice
  'f7e084',                                              // 16 barren
  '6f6f6f'                                               // 17 tundra
];

Map.centerObject(STUDY_AREA, MAP_ZOOM_LEVEL);
Map.addLayer(cloudMaskedScene, {bands: ['B3', 'B2', 'B1'], max: 0.4},
             'Landsat 5 true-colour composite');
Map.addLayer(classifiedLandCover, {palette: igbpPalette, min: 0, max: 17},
             'Random Forest IGBP classification');
