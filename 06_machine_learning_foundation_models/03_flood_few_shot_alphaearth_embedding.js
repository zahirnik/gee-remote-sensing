/*
 * Few-shot flood mapping from AlphaEarth pre / post embeddings
 *
 * Datasets:    GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL (AlphaEarth embeddings),
 *              COPERNICUS/S1_GRD (Sentinel-1 GRD, used as a SAR baseline),
 *              JRC/GSW1_4/MonthlyHistory (JRC Global Surface Water — to mask
 *              out permanent water bodies so we only score true flood pixels).
 * Region:      A small AOI around Beira, Mozambique — the impact zone of
 *              Cyclone Idai (March 2019).
 * Output:      Three flood masks (prototype-distance, few-shot RF on
 *              Δembedding, Sentinel-1 SAR change-detection baseline),
 *              their pairwise agreement, plus a Drive Export task.
 *
 * Usage:
 *   1. Paste this script into the GEE Code Editor.
 *   2. The default labelled points (FLOOD_TRAINING_POINTS and
 *      NON_FLOOD_TRAINING_POINTS) sit in known flooded / non-flooded areas
 *      from Cyclone Idai. Move them or add more in the Code Editor to
 *      experiment with the few-shot setting.
 *   3. Click Run. The Map shows the three flood-mask layers; the Console
 *      reports pixel counts + pairwise agreement; an Export task is queued.
 *
 * What this script demonstrates
 * -----------------------------
 *   - Few-shot disaster mapping using a foundation model — we use only a
 *     handful of labelled points to either compute an embedding prototype
 *     (no training) or to train a tiny Random Forest.
 *   - Δembedding = post-event - pre-event captures CHANGE rather than state,
 *     which is the right signal for an event-based hazard like flooding.
 *   - Direct comparison with a conventional Sentinel-1 SAR pre / post
 *     change-detection baseline as a sanity-check.
 */


// ---------------------------------------------------------------------------
// CONFIG — edit these constants to retarget the analysis.
// ---------------------------------------------------------------------------

// Study AOI — wider region around Beira, Mozambique. Cyclone Idai made
// landfall on 14-15 March 2019.
var STUDY_AREA = ee.Geometry.Rectangle([34.70, -20.10, 35.20, -19.60]);

// AlphaEarth is annual, so we use the year before the event as "pre" and
// the year of the event as "post". For finer temporal change you could
// substitute a different sensor (Sentinel-1) but this script's point is to
// demonstrate the foundation-embedding-based approach.
var PRE_EVENT_YEAR  = 2018;
var POST_EVENT_YEAR = 2019;

// Sentinel-1 baseline: pre/post date windows for the SAR change-detection.
var SAR_PRE_DATE_START  = '2019-02-15';
var SAR_PRE_DATE_END    = '2019-03-10';
var SAR_POST_DATE_START = '2019-03-15';
var SAR_POST_DATE_END   = '2019-04-15';

// SAR change-detection threshold: how many dB darker the post image must be
// than the pre image at a pixel to be flagged as "flooded".
var SAR_DECIBEL_DROP_THRESHOLD = 1.25;

// Slope mask (degrees): flat areas only — flooding doesn't happen on hills.
var SLOPE_MASK_DEGREES = 5;

// Cosine-similarity threshold for the prototype-distance method. Pixels
// whose Δembedding has cosine-similarity to the mean flood prototype above
// this threshold are flagged as "flooded".
var PROTOTYPE_SIMILARITY_THRESHOLD = 0.6;

// Random Forest hyperparameters for the few-shot approach.
var NUM_TREES_RF = 50;
var SEED         = 7;

// Export.
var EXPORT_SCALE_METERS = 10;
var DRIVE_FOLDER_NAME   = 'GEE_exports';


// ---------------------------------------------------------------------------
// 1. Few-shot labelled points (edit / extend these in the Code Editor).
// ---------------------------------------------------------------------------

// FLOOD_TRAINING_POINTS sit on flooded land per public CEMS damage products
// for the Idai event. NON_FLOOD_TRAINING_POINTS are well clear of the
// inundation zone. Move them around the map to see how the few-shot model
// responds.
var FLOOD_TRAINING_POINTS = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([34.79, -19.83]), {label: 1}),
  ee.Feature(ee.Geometry.Point([34.83, -19.85]), {label: 1}),
  ee.Feature(ee.Geometry.Point([34.86, -19.79]), {label: 1}),
  ee.Feature(ee.Geometry.Point([34.90, -19.82]), {label: 1}),
  ee.Feature(ee.Geometry.Point([34.94, -19.85]), {label: 1}),
  ee.Feature(ee.Geometry.Point([34.97, -19.88]), {label: 1})
]);

var NON_FLOOD_TRAINING_POINTS = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([34.75, -19.65]), {label: 0}),
  ee.Feature(ee.Geometry.Point([34.78, -19.70]), {label: 0}),
  ee.Feature(ee.Geometry.Point([35.10, -20.00]), {label: 0}),
  ee.Feature(ee.Geometry.Point([35.15, -19.95]), {label: 0}),
  ee.Feature(ee.Geometry.Point([35.05, -19.65]), {label: 0}),
  ee.Feature(ee.Geometry.Point([35.00, -19.70]), {label: 0})
]);


// ---------------------------------------------------------------------------
// 2. Load AlphaEarth pre and post embeddings, build Δembedding.
// ---------------------------------------------------------------------------

var alphaEarthCollection = ee.ImageCollection('GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL')
  .filterBounds(STUDY_AREA);

function loadAnnualEmbeddingMosaic(year) {
  var start = ee.Date.fromYMD(year, 1, 1);
  var end   = start.advance(1, 'year');
  return alphaEarthCollection
    .filterDate(start, end)
    .mosaic()
    .clip(STUDY_AREA);
}

var preEventEmbedding  = loadAnnualEmbeddingMosaic(PRE_EVENT_YEAR);
var postEventEmbedding = loadAnnualEmbeddingMosaic(POST_EVENT_YEAR);

// Δembedding = post - pre, per band. Rename to make features unambiguous.
var deltaEmbeddingImage = postEventEmbedding
  .subtract(preEventEmbedding)
  .rename(preEventEmbedding.bandNames().map(function (b) {
    return ee.String('d').cat(ee.String(b));
  }));

var DELTA_BAND_NAMES = deltaEmbeddingImage.bandNames();


// ---------------------------------------------------------------------------
// 3. Method A — prototype distance (NO training).
// ---------------------------------------------------------------------------
// We compute the mean Δembedding over the flood training points. Any pixel
// whose Δembedding has high cosine-similarity to that prototype is flagged
// as flooded. This is the closest thing to "zero-shot" we can do.

var floodPrototypeMeansDictionary = deltaEmbeddingImage
  .reduceRegions({
    collection: FLOOD_TRAINING_POINTS,
    reducer:    ee.Reducer.mean(),
    scale:      EXPORT_SCALE_METERS
  })
  .aggregate_array('first');   // placeholder — we re-pivot below

// The above doesn't quite give us a flat 64-vector. We do it explicitly:
var floodPrototypeVector = deltaEmbeddingImage
  .reduceRegion({
    reducer:  ee.Reducer.mean(),
    geometry: FLOOD_TRAINING_POINTS.geometry(),
    scale:    EXPORT_SCALE_METERS,
    maxPixels: 1e7
  });

// Build a constant image with the prototype vector tiled to every pixel.
// This gives us an image we can take the dot-product against.
var floodPrototypeImage = ee.Image.constant(
  DELTA_BAND_NAMES.map(function (bandName) {
    return ee.Number(floodPrototypeVector.get(bandName));
  })
).rename(DELTA_BAND_NAMES);

// Cosine similarity = (a . b) / (||a|| * ||b||).
function safeMagnitude(image) {
  return image.pow(2).reduce(ee.Reducer.sum()).sqrt();
}

var dotProductImage = deltaEmbeddingImage
  .multiply(floodPrototypeImage)
  .reduce(ee.Reducer.sum());

var cosineSimilarityImage = dotProductImage.divide(
  safeMagnitude(deltaEmbeddingImage).multiply(safeMagnitude(floodPrototypeImage)).add(1e-9)
);

var prototypeBasedFloodMask = cosineSimilarityImage
  .gte(PROTOTYPE_SIMILARITY_THRESHOLD)
  .rename('flood_prototype');


// ---------------------------------------------------------------------------
// 4. Method B — few-shot Random Forest on Δembedding.
// ---------------------------------------------------------------------------

var fewShotTrainingFeatures = deltaEmbeddingImage
  .sampleRegions({
    collection: FLOOD_TRAINING_POINTS.merge(NON_FLOOD_TRAINING_POINTS),
    properties: ['label'],
    scale:      EXPORT_SCALE_METERS,
    geometries: false
  });

var fewShotRfClassifier = ee.Classifier.smileRandomForest({
    numberOfTrees: NUM_TREES_RF,
    seed:          SEED
  })
  .train({
    features:        fewShotTrainingFeatures,
    classProperty:   'label',
    inputProperties: DELTA_BAND_NAMES
  });

var fewShotRfFloodMask = deltaEmbeddingImage
  .classify(fewShotRfClassifier)
  .rename('flood_few_shot_rf');


// ---------------------------------------------------------------------------
// 5. Method C — Sentinel-1 SAR pre/post change-detection baseline.
// ---------------------------------------------------------------------------

var s1Collection = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(STUDY_AREA)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .select('VH');

var sarPreMedian  = s1Collection
  .filterDate(SAR_PRE_DATE_START,  SAR_PRE_DATE_END).median();
var sarPostMedian = s1Collection
  .filterDate(SAR_POST_DATE_START, SAR_POST_DATE_END).median();

// Speckle filter: small focal mean to suppress noise before differencing.
var sarPreSmoothed  = sarPreMedian.focal_mean({radius: 50, units: 'meters'});
var sarPostSmoothed = sarPostMedian.focal_mean({radius: 50, units: 'meters'});

var sarDifferenceImage = sarPreSmoothed.subtract(sarPostSmoothed);

// Slope mask using SRTM — exclude steep terrain where SAR change is dominated
// by foreshortening rather than real surface change.
var srtmElevationImage = ee.Image('USGS/SRTMGL1_003');
var srtmSlopeImage     = ee.Terrain.slope(srtmElevationImage);
var flatTerrainMask    = srtmSlopeImage.lt(SLOPE_MASK_DEGREES);

// Permanent-water mask: exclude pixels that JRC GSW reports as water in the
// month before the event (we don't want to flag the ocean as "flooded").
var jrcMonthlyWater = ee.ImageCollection('JRC/GSW1_4/MonthlyHistory')
  .filter(ee.Filter.eq('year',  PRE_EVENT_YEAR))
  .filter(ee.Filter.eq('month', 12))
  .first();
var permanentWaterMask = ee.Image(jrcMonthlyWater).eq(2);   // class 2 = water

var sarBaselineFloodMask = sarDifferenceImage
  .gt(SAR_DECIBEL_DROP_THRESHOLD)
  .and(flatTerrainMask)
  .and(permanentWaterMask.not())
  .rename('flood_sar_baseline');


// ---------------------------------------------------------------------------
// 6. Agreement summary and visualisation.
// ---------------------------------------------------------------------------

var allMethodsAgreementImage = prototypeBasedFloodMask
  .add(fewShotRfFloodMask)
  .add(sarBaselineFloodMask)
  .rename('methods_in_agreement');

print('Flood-pixel counts by method:');
[prototypeBasedFloodMask, fewShotRfFloodMask, sarBaselineFloodMask]
  .forEach(function (mask) {
    print(mask.bandNames().get(0),
          mask.reduceRegion({
            reducer:  ee.Reducer.sum(),
            geometry: STUDY_AREA,
            scale:    EXPORT_SCALE_METERS,
            maxPixels: 1e10
          }));
  });

Map.centerObject(STUDY_AREA, 10);
Map.addLayer(prototypeBasedFloodMask.selfMask(), {palette: ['1f78b4']},
             'A: prototype-distance flood');
Map.addLayer(fewShotRfFloodMask.selfMask(),      {palette: ['33a02c']},
             'B: few-shot RF flood');
Map.addLayer(sarBaselineFloodMask.selfMask(),    {palette: ['e31a1c']},
             'C: SAR baseline flood');
Map.addLayer(allMethodsAgreementImage,
             {min: 0, max: 3, palette: ['ffffff', 'fed976', 'fd8d3c', 'bd0026']},
             'Number of methods in agreement');
Map.addLayer(FLOOD_TRAINING_POINTS,
             {color: '1f78b4'}, 'Flood training points');
Map.addLayer(NON_FLOOD_TRAINING_POINTS,
             {color: '999999'}, 'Non-flood training points');


// ---------------------------------------------------------------------------
// 7. Export the agreement image to Drive.
// ---------------------------------------------------------------------------

Export.image.toDrive({
  image:       allMethodsAgreementImage.toUint8(),
  description: 'flood_agreement_alphaearth_few_shot',
  folder:      DRIVE_FOLDER_NAME,
  region:      STUDY_AREA,
  scale:       EXPORT_SCALE_METERS,
  maxPixels:   1e10
});
