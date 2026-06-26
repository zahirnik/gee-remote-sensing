/*
 * Crop-type classification with AlphaEarth foundation embeddings + Random Forest
 *
 * Datasets:    GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL (AlphaEarth annual embeddings),
 *              USDA/NASS/CDL (USDA Cropland Data Layer — pixel-level crop labels)
 * Region:      A small AOI in Story County, Iowa (USA Midwest corn / soy belt)
 * Output:      Classified crop-type map for the AOI, training- and validation-
 *              accuracy in the console, confusion matrix, and an EE Export task
 *              that writes the classified raster to Drive.
 *
 * Usage:
 *   1. Paste this script into the GEE Code Editor.
 *   2. Adjust the constants in the CONFIG block below if you want a different
 *      study area or year. The defaults (Story County, Iowa, 2022) run in ~30 s
 *      and require no extra data.
 *   3. Click Run. The Map shows the classified raster; the Console prints
 *      train/test accuracy + a confusion matrix; an Export task is queued
 *      that you can confirm under the Tasks tab to save the GeoTIFF to Drive.
 *
 * What this script demonstrates
 * -----------------------------
 *   - How to load AlphaEarth's 64-dimensional learned pixel embeddings
 *     (`A00`..`A63`) for a chosen year and AOI.
 *   - How to use those embeddings *directly* as features for a Random Forest
 *     classifier, with USDA Cropland Data Layer as labels.
 *   - How to do stratified training / validation sampling and report a clean
 *     accuracy + confusion matrix instead of a single number.
 *
 * Why this is interesting
 * -----------------------
 *   AlphaEarth is a foundation model trained on multi-year, multi-sensor
 *   Earth-observation data. Each pixel is summarised as a 64-dim embedding
 *   that captures spectral, temporal, and spatial context. As features for
 *   a downstream classifier, embeddings typically need far fewer labelled
 *   samples than raw spectral bands to reach the same accuracy.
 */


// ---------------------------------------------------------------------------
// CONFIG — edit these constants to retarget the study.
// ---------------------------------------------------------------------------

// Study AOI: a 0.4 deg x 0.3 deg rectangle around Story County, Iowa (USA).
// You can also pass any ee.Geometry / ee.FeatureCollection here.
var STUDY_AREA = ee.Geometry.Rectangle([-93.7, 41.85, -93.3, 42.15]);

// Year to classify. AlphaEarth covers 2017..2024 inclusive at time of writing.
var TARGET_YEAR = 2022;

// USDA CDL classes to keep (everything else is collapsed into 'OTHER' = 0).
// Codes come from https://www.nass.usda.gov/Research_and_Science/Cropland/
//   1   = Corn
//   5   = Soybeans
//   24  = Winter Wheat
// We renumber the kept classes to 1, 2, 3 so the classifier output is compact.
var CROP_CODES_OF_INTEREST = {
  1: 1,   // Corn        -> class 1
  5: 2,   // Soybeans    -> class 2
  24: 3   // Winter Wheat-> class 3
};

// Per-class training / validation point counts. We sample BOTH classes
// independently to avoid the dominant class swamping the classifier.
var POINTS_PER_CLASS_TRAINING   = 500;
var POINTS_PER_CLASS_VALIDATION = 200;

// Random Forest hyperparameters.
var NUM_TREES_RF = 100;
var BAG_FRACTION = 0.5;
var SEED         = 42;

// Export resolution (m). AlphaEarth is 10 m; we keep that to honour it.
var EXPORT_SCALE_METERS = 10;

// Folder name inside your Google Drive where the Export task will save.
var DRIVE_FOLDER_NAME = 'GEE_exports';


// ---------------------------------------------------------------------------
// 1. Load the AlphaEarth annual embedding for the target year and AOI.
// ---------------------------------------------------------------------------

// AlphaEarth is published as one tile per year per geographic tile, so we
// filter by date AND by bounds, then mosaic the tiles into a single image.
var startDate = ee.Date.fromYMD(TARGET_YEAR, 1, 1);
var endDate   = startDate.advance(1, 'year');

var alphaEarthEmbedding = ee.ImageCollection('GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL')
  .filterBounds(STUDY_AREA)
  .filterDate(startDate, endDate)
  .mosaic()
  .clip(STUDY_AREA);

// Embedding bands are named A00 .. A63. Cache the list once because we re-use
// it as the feature list for the classifier.
var EMBEDDING_BAND_NAMES = alphaEarthEmbedding.bandNames();
print('Number of AlphaEarth embedding bands:', EMBEDDING_BAND_NAMES.size());


// ---------------------------------------------------------------------------
// 2. Load USDA CDL labels for the same year and re-map kept classes.
// ---------------------------------------------------------------------------

var rawCdlImage = ee.ImageCollection('USDA/NASS/CDL')
  .filterDate(startDate, endDate)
  .first()                                     // single image per year
  .select('cropland')
  .clip(STUDY_AREA);

// Re-map: kept codes become 1..N; everything else becomes 0 (= "other").
// We build the lookup table from the CROP_CODES_OF_INTEREST dict above.
var sourceCodes = Object.keys(CROP_CODES_OF_INTEREST).map(function (k) {
  return parseInt(k, 10);
});
var targetCodes = sourceCodes.map(function (k) {
  return CROP_CODES_OF_INTEREST[k];
});
var cropLabelImage = rawCdlImage
  .remap(sourceCodes, targetCodes, 0)          // unmapped -> 0
  .rename('cropClass');


// ---------------------------------------------------------------------------
// 3. Stratified sampling of training + validation pixels.
// ---------------------------------------------------------------------------

// We build a single sample set per class (training + validation lumped), then
// split it 70 / 30 with a deterministic random column for repeatability.
var sampledFeatures = alphaEarthEmbedding
  .addBands(cropLabelImage)
  .stratifiedSample({
    numPoints: POINTS_PER_CLASS_TRAINING + POINTS_PER_CLASS_VALIDATION,
    classBand: 'cropClass',
    region:    STUDY_AREA,
    scale:     EXPORT_SCALE_METERS,
    seed:      SEED,
    geometries: false                          // we only need the values
  })
  // Deterministic 0..1 random column so the same code reproduces the split.
  .randomColumn('split_seed', SEED);

// 70 % of the sampled pixels -> training; 30 % -> validation.
var SPLIT_THRESHOLD = 0.7;
var trainingFeatures   = sampledFeatures.filter(ee.Filter.lt('split_seed', SPLIT_THRESHOLD));
var validationFeatures = sampledFeatures.filter(ee.Filter.gte('split_seed', SPLIT_THRESHOLD));

print('Training pixels:',   trainingFeatures.size());
print('Validation pixels:', validationFeatures.size());


// ---------------------------------------------------------------------------
// 4. Train a Random Forest on the 64 AlphaEarth bands.
// ---------------------------------------------------------------------------

var randomForestClassifier = ee.Classifier.smileRandomForest({
    numberOfTrees: NUM_TREES_RF,
    bagFraction:   BAG_FRACTION,
    seed:          SEED
  })
  .train({
    features:       trainingFeatures,
    classProperty:  'cropClass',
    inputProperties: EMBEDDING_BAND_NAMES
  });

// Classifier "explain" gives feature importances. Useful for sanity-checking
// that the model is actually using a diverse set of embedding dimensions.
var classifierExplain = randomForestClassifier.explain();
print('Classifier explanation:', classifierExplain);


// ---------------------------------------------------------------------------
// 5. Apply the classifier and evaluate.
// ---------------------------------------------------------------------------

// (a) Pixel-level classification of the AOI.
var classifiedCropImage = alphaEarthEmbedding
  .classify(randomForestClassifier)
  .rename('predictedClass');

// (b) Validation: predict on held-out features, build a confusion matrix.
var validationResults = validationFeatures.classify(randomForestClassifier);
var confusionMatrix = validationResults.errorMatrix('cropClass', 'classification');
print('Validation confusion matrix:', confusionMatrix);
print('Validation overall accuracy:', confusionMatrix.accuracy());
print('Validation kappa:',            confusionMatrix.kappa());


// ---------------------------------------------------------------------------
// 6. Visualise on the map.
// ---------------------------------------------------------------------------

// Display palette: index 0 = OTHER, 1 = corn (gold), 2 = soy (green),
// 3 = winter wheat (orange).
var CROP_DISPLAY_PALETTE = ['cccccc', 'ffd400', '2ca02c', 'd95f02'];
var CROP_DISPLAY_LABELS  = ['Other', 'Corn', 'Soy', 'Winter Wheat'];

Map.centerObject(STUDY_AREA, 11);
Map.addLayer(cropLabelImage.selfMask(),
             {min: 0, max: 3, palette: CROP_DISPLAY_PALETTE},
             'CDL labels');
Map.addLayer(classifiedCropImage.selfMask(),
             {min: 0, max: 3, palette: CROP_DISPLAY_PALETTE},
             'RF classification (AlphaEarth)');


// ---------------------------------------------------------------------------
// 7. Export the classified raster to Drive.
// ---------------------------------------------------------------------------

Export.image.toDrive({
  image:       classifiedCropImage.toUint8(),
  description: 'crop_classification_alphaearth_rf_' + TARGET_YEAR,
  folder:      DRIVE_FOLDER_NAME,
  region:      STUDY_AREA,
  scale:       EXPORT_SCALE_METERS,
  maxPixels:   1e10
});
