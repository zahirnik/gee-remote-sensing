/*
 * Drought-severity classification with AlphaEarth embedding trajectories
 * + CHIRPS SPI labels + Gradient-Boosted Trees
 *
 * Datasets:    GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL (AlphaEarth embeddings),
 *              UCSB-CHG/CHIRPS/DAILY (CHIRPS daily precipitation, 0.05 deg)
 * Region:      A small AOI in California's Central Valley (USA)
 * Output:      Predicted drought-severity class map (D0 / D1 / D2 / D3 / D4)
 *              for the target year, training / validation accuracy, confusion
 *              matrix, and a Drive Export task for the GeoTIFF.
 *
 * Usage:
 *   1. Paste this script into the GEE Code Editor.
 *   2. Edit the constants in the CONFIG block if you want a different AOI
 *      or year. Defaults run in ~1 min.
 *   3. Click Run. Map shows the predicted drought classes; Console prints
 *      class counts, accuracy, and the confusion matrix; Export task queued.
 *
 * What this script demonstrates
 * -----------------------------
 *   - Using AlphaEarth embeddings as a *multi-year* feature stack: we take
 *     the embedding for year T-1 and the change in embedding since T-3, so
 *     the classifier learns the *trajectory* not just the current state.
 *   - Generating drought-severity labels in-pure-GEE from CHIRPS rainfall
 *     using a percentile-binned SPI proxy.
 *   - Gradient-Boosted Trees (`ee.Classifier.smileGradientTreeBoost`) which
 *     usually outperforms RF for ordinal targets like drought severity.
 */


// ---------------------------------------------------------------------------
// CONFIG — edit these constants to retarget the study.
// ---------------------------------------------------------------------------

// Study AOI: a portion of California's Central Valley.
var STUDY_AREA = ee.Geometry.Rectangle([-121.5, 36.5, -120.5, 37.5]);

// The target year is the year we want to *predict* drought severity for.
// Features come from years TARGET_YEAR-1 and TARGET_YEAR-3.
var TARGET_YEAR = 2022;

// CHIRPS lookback to compute the SPI label: how many years of rainfall to use
// when defining each pixel's percentile reference distribution.
var REFERENCE_PERIOD_START_YEAR = 1990;
var REFERENCE_PERIOD_END_YEAR   = 2020;

// SPI breakpoints (expressed as cumulative-rainfall percentiles within the
// reference period) for the U.S. Drought Monitor severity classes:
//   D0 - Abnormally Dry             (<= 30 percentile)
//   D1 - Moderate Drought           (<= 20 percentile)
//   D2 - Severe Drought             (<= 10 percentile)
//   D3 - Extreme Drought            (<= 5  percentile)
//   D4 - Exceptional Drought        (<= 2  percentile)
// Pixels above the 30th percentile are labelled 0 (no drought).
var DROUGHT_CLASS_BREAKPOINTS = [2, 5, 10, 20, 30];

// Sampling.
var POINTS_PER_CLASS_TRAINING   = 400;
var POINTS_PER_CLASS_VALIDATION = 150;

// Gradient-Boosted Trees hyperparameters.
var NUM_TREES_GBT  = 200;
var SHRINKAGE_GBT  = 0.05;
var MAX_NODES_GBT  = 32;
var SAMPLING_FRAC  = 0.7;
var SEED           = 17;

// Export.
var EXPORT_SCALE_METERS = 30;          // CHIRPS is coarse anyway; 30 m is fine.
var DRIVE_FOLDER_NAME   = 'GEE_exports';


// ---------------------------------------------------------------------------
// 1. Load AlphaEarth embeddings for years TARGET_YEAR-1 and TARGET_YEAR-3.
// ---------------------------------------------------------------------------

var alphaEarthCollection = ee.ImageCollection('GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL')
  .filterBounds(STUDY_AREA);

// Helper: get one mosaiced embedding image for a given year. We rename bands
// with a prefix so we can stack two years without name collisions.
function getAnnualEmbeddingMosaic(year, bandPrefix) {
  var start = ee.Date.fromYMD(year, 1, 1);
  var end   = start.advance(1, 'year');
  var mosaic = alphaEarthCollection
    .filterDate(start, end)
    .mosaic()
    .clip(STUDY_AREA);
  // Rename A00..A63 -> <prefix>A00..<prefix>A63.
  var originalBandNames = mosaic.bandNames();
  var renamedBandNames  = originalBandNames.map(function (b) {
    return ee.String(bandPrefix).cat(ee.String(b));
  });
  return mosaic.rename(renamedBandNames);
}

var embeddingsYearMinus1 = getAnnualEmbeddingMosaic(TARGET_YEAR - 1, 'lag1_');
var embeddingsYearMinus3 = getAnnualEmbeddingMosaic(TARGET_YEAR - 3, 'lag3_');

// Δembedding = lag1 - lag3, with band-name prefix 'delta_'. We use map() on
// band indices so the subtraction lines up by ordinal position.
var bandIndices = ee.List.sequence(0, 63);
var deltaEmbeddingImage = ee.ImageCollection.fromImages(bandIndices.map(function (i) {
  i = ee.Number(i).toInt();
  var lagBandSuffix = ee.String('A').cat(i.format('%02d'));
  var lag1Band = embeddingsYearMinus1.select(ee.String('lag1_').cat(lagBandSuffix));
  var lag3Band = embeddingsYearMinus3.select(ee.String('lag3_').cat(lagBandSuffix));
  return lag1Band.subtract(lag3Band)
                 .rename(ee.String('delta_').cat(lagBandSuffix));
})).toBands();

// Strip the auto-generated band-index prefix that toBands() adds.
deltaEmbeddingImage = deltaEmbeddingImage.rename(
  deltaEmbeddingImage.bandNames().map(function (b) {
    return ee.String(b).slice(ee.String(b).index('delta_'));
  })
);


// ---------------------------------------------------------------------------
// 2. Stack lag1 + Δembedding into the final feature image (128 bands).
// ---------------------------------------------------------------------------

var featureStack = embeddingsYearMinus1.addBands(deltaEmbeddingImage);
var FEATURE_BAND_NAMES = featureStack.bandNames();
print('Number of features used:', FEATURE_BAND_NAMES.size());


// ---------------------------------------------------------------------------
// 3. Generate the drought-severity label from CHIRPS rainfall.
// ---------------------------------------------------------------------------

// Helper: cumulative water-year rainfall for a given year (Oct -> Sep).
// We use the "previous October to September" convention for SPI-12.
// We accept either a plain JS number or an ee.Number (the latter happens
// when this helper is called from inside ee.List.sequence().map(...)), so
// we wrap the input in ee.Number and use server-side .subtract(...).
function getWaterYearTotalRainfall(year) {
  var yearAsEeNumber = ee.Number(year);
  var start = ee.Date.fromYMD(yearAsEeNumber.subtract(1), 10, 1);
  var end   = ee.Date.fromYMD(yearAsEeNumber,             10, 1);
  return ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
    .filterDate(start, end)
    .filterBounds(STUDY_AREA)
    .sum()
    .clip(STUDY_AREA)
    .rename('totalRainfall');
}

// Reference distribution of total annual rainfall per pixel across the
// reference period. We compute the percentile breakpoints once and reuse
// them as label thresholds.
var referenceYears = ee.List.sequence(REFERENCE_PERIOD_START_YEAR,
                                      REFERENCE_PERIOD_END_YEAR);
var referenceCollection = ee.ImageCollection.fromImages(
  referenceYears.map(function (y) { return getWaterYearTotalRainfall(y); })
);

// Compute per-pixel percentiles corresponding to our breakpoints.
var percentileImage = referenceCollection.reduce(
  ee.Reducer.percentile(DROUGHT_CLASS_BREAKPOINTS)
);

// Total rainfall for the TARGET_YEAR.
var targetYearRainfall = getWaterYearTotalRainfall(TARGET_YEAR);

// Convert rainfall to a categorical drought-class image:
//   if rainfall <= p2  -> 5 (D4)
//   elif <= p5         -> 4 (D3)
//   elif <= p10        -> 3 (D2)
//   elif <= p20        -> 2 (D1)
//   elif <= p30        -> 1 (D0)
//   else               -> 0 (no drought)
// We compute it as a sum of binary masks, which is fast on EE.
var droughtSeverityLabel = ee.Image(0)
  .where(targetYearRainfall.lte(percentileImage.select('totalRainfall_p30')), 1)
  .where(targetYearRainfall.lte(percentileImage.select('totalRainfall_p20')), 2)
  .where(targetYearRainfall.lte(percentileImage.select('totalRainfall_p10')), 3)
  .where(targetYearRainfall.lte(percentileImage.select('totalRainfall_p5')),  4)
  .where(targetYearRainfall.lte(percentileImage.select('totalRainfall_p2')),  5)
  .rename('droughtClass')
  .clip(STUDY_AREA);


// ---------------------------------------------------------------------------
// 4. Sample training + validation pixels stratified by drought class.
// ---------------------------------------------------------------------------

var sampledFeatures = featureStack
  .addBands(droughtSeverityLabel)
  .stratifiedSample({
    numPoints: POINTS_PER_CLASS_TRAINING + POINTS_PER_CLASS_VALIDATION,
    classBand: 'droughtClass',
    region:    STUDY_AREA,
    scale:     EXPORT_SCALE_METERS,
    seed:      SEED,
    geometries: false
  })
  .randomColumn('split_seed', SEED);

var SPLIT_THRESHOLD = 0.7;
var trainingFeatures   = sampledFeatures.filter(ee.Filter.lt('split_seed', SPLIT_THRESHOLD));
var validationFeatures = sampledFeatures.filter(ee.Filter.gte('split_seed', SPLIT_THRESHOLD));

print('Sampled training pixels:',   trainingFeatures.size());
print('Sampled validation pixels:', validationFeatures.size());


// ---------------------------------------------------------------------------
// 5. Train a Gradient-Boosted Trees classifier on the 128-feature stack.
// ---------------------------------------------------------------------------

var gradientBoostedTreesClassifier = ee.Classifier.smileGradientTreeBoost({
    numberOfTrees: NUM_TREES_GBT,
    shrinkage:     SHRINKAGE_GBT,
    samplingRate:  SAMPLING_FRAC,
    maxNodes:      MAX_NODES_GBT,
    seed:          SEED
  })
  .train({
    features:        trainingFeatures,
    classProperty:   'droughtClass',
    inputProperties: FEATURE_BAND_NAMES
  });


// ---------------------------------------------------------------------------
// 6. Evaluate and apply to the AOI.
// ---------------------------------------------------------------------------

var predictedDroughtImage = featureStack
  .classify(gradientBoostedTreesClassifier)
  .rename('predictedDroughtClass');

var validationPredictions = validationFeatures.classify(gradientBoostedTreesClassifier);
var droughtConfusionMatrix = validationPredictions
  .errorMatrix('droughtClass', 'classification');

print('Drought confusion matrix:', droughtConfusionMatrix);
print('Overall accuracy:',         droughtConfusionMatrix.accuracy());
print('Kappa:',                    droughtConfusionMatrix.kappa());


// ---------------------------------------------------------------------------
// 7. Visualise on the map.
// ---------------------------------------------------------------------------

// Colour ramp ordered from "no drought" (white) to D4 (deep red).
var DROUGHT_DISPLAY_PALETTE = ['ffffff', 'fef0d9', 'fdcc8a', 'fc8d59', 'e34a33', 'b30000'];

Map.centerObject(STUDY_AREA, 9);
Map.addLayer(droughtSeverityLabel,
             {min: 0, max: 5, palette: DROUGHT_DISPLAY_PALETTE},
             'CHIRPS-derived drought label');
Map.addLayer(predictedDroughtImage,
             {min: 0, max: 5, palette: DROUGHT_DISPLAY_PALETTE},
             'GBT prediction (AlphaEarth)');


// ---------------------------------------------------------------------------
// 8. Export the predicted drought-class raster to Drive.
// ---------------------------------------------------------------------------

Export.image.toDrive({
  image:       predictedDroughtImage.toUint8(),
  description: 'drought_severity_alphaearth_gbt_' + TARGET_YEAR,
  folder:      DRIVE_FOLDER_NAME,
  region:      STUDY_AREA,
  scale:       EXPORT_SCALE_METERS,
  maxPixels:   1e10
});
