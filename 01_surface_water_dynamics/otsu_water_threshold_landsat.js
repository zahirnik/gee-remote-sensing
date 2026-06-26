/*
 * Otsu between-class-variance water classification on a Landsat 5/7/8 composite
 *
 * Datasets:    LANDSAT/LT05/C01/T1_SR, LANDSAT/LE07/C01/T1_SR, LANDSAT/LC08/C01/T1_SR
 * Region:      Any user-defined polygon (default: Pyramid Lake, NV, USA)
 * Output:      Map layers: false-colour Landsat composite, NIR-band histogram,
 *              binary water mask classified at the Otsu-optimal threshold.
 *
 * Usage:
 *   1. Paste this script into the GEE Code Editor.
 *   2. Replace STUDY_AREA polygon and DATE_RANGE constants below.
 *   3. Click Run; the classified water layer and the BSS curve chart appear
 *      in the Map and Console panels respectively.
 *
 * Method:
 *   Otsu (1979) selects the digital number that maximises the between-class
 *   variance of the histogram, i.e. the value that best separates the bimodal
 *   distribution of dark (water) and bright (non-water) NIR pixels.
 */

// ------------------------------------------------------------------
// 0. User-editable constants
// ------------------------------------------------------------------
var STUDY_AREA = ee.Geometry.Polygon(
  [[[-119.8238273926751, 40.24346421078928],
    [-119.8238273926751, 39.76589137578468],
    [-119.1371818848626, 39.76589137578468],
    [-119.1371818848626, 40.24346421078928]]], null, false);

var START_DATE   = '2000-01-01';
var END_DATE     = '2000-04-01';   // a short window keeps the histogram bimodal
var HISTOGRAM_BANDS = ['B4'];      // NIR-equivalent band (renamed below)
var HISTOGRAM_BUCKET_WIDTH = 2;
var HISTOGRAM_BUCKET_COUNT = 255;

// ------------------------------------------------------------------
// 1. Cloud / shadow / snow masking helpers for Landsat C01 SR products
// ------------------------------------------------------------------
// L4/5/7 pre-Collection-2 cloud bits: bit5=cloud, bit7=high-confidence cloud,
// bit3=cloud shadow, bit4=snow. We mark a pixel bad if it is a confident cloud
// OR a (cloud shadow + snow) combination, then drop edge pixels that are
// masked in any band.
var maskLandsat457 = function (image) {
  var qa = image.select('pixel_qa');
  var cloud = qa.bitwiseAnd(1 << 5)
    .and(qa.bitwiseAnd(1 << 7))
    .or(qa.bitwiseAnd(1 << 3))
    .and(qa.bitwiseAnd(1 << 4));
  var edgeMask = image.mask().reduce(ee.Reducer.min());
  return image.updateMask(cloud.not()).updateMask(edgeMask);
};

// L8 SR uses simpler bit semantics: bit3=cloud shadow, bit4=snow, bit5=cloud.
var maskLandsat8 = function (image) {
  var qa = image.select('pixel_qa');
  var mask = qa.bitwiseAnd(1 << 3).eq(0)
    .and(qa.bitwiseAnd(1 << 5).eq(0))
    .and(qa.bitwiseAnd(1 << 4).eq(0));
  return image.updateMask(mask);
};

// ------------------------------------------------------------------
// 2. Build a harmonised Landsat 5/7/8 collection
// ------------------------------------------------------------------
// We rename L8 bands so blue/green/red/NIR/SWIR1/SWIR2 sit on the same
// band-name spots as L5/7. That way later code does not need to branch
// on which sensor produced each scene.
var landsat5 = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')
  .select(['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'pixel_qa'])
  .map(maskLandsat457)
  .map(function (image) { return image.clip(STUDY_AREA); });

var landsat7 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
  .select(['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'pixel_qa'])
  .map(maskLandsat457)
  .map(function (image) { return image.clip(STUDY_AREA); });

var landsat8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
  .map(maskLandsat8)
  .select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'pixel_qa'])
  .map(function (image) { return image.clip(STUDY_AREA); })
  .map(function (image) {
    return image.rename(['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'pixel_qa']);
  });

var mergedLandsat = ee.ImageCollection(
  landsat5.merge(landsat7).merge(landsat8)
).filterDate(START_DATE, END_DATE);

// ------------------------------------------------------------------
// 3. Build a single representative composite for thresholding
// ------------------------------------------------------------------
// Median over the full date window is robust to residual cloud / shadow
// pixels that survived the QA mask.
var landsatComposite = mergedLandsat.median().clip(STUDY_AREA);

Map.centerObject(STUDY_AREA, 10);
Map.addLayer(
  landsatComposite,
  { bands: ['B4', 'B3', 'B2'], min: 5000, max: 15000, gamma: 1.3 },
  'Landsat false-colour composite'
);

// ------------------------------------------------------------------
// 4. Compute the NIR-band histogram over the study area
// ------------------------------------------------------------------
// The NIR band is the strongest single-band water/land discriminator on
// Landsat (water absorbs almost all NIR, vegetation reflects strongly).
var nirHistogram = landsatComposite.select(HISTOGRAM_BANDS).reduceRegion({
  reducer: ee.Reducer.histogram(HISTOGRAM_BUCKET_COUNT, HISTOGRAM_BUCKET_WIDTH)
    .combine('mean', null, true)
    .combine('variance', null, true),
  geometry: STUDY_AREA,
  scale: 30,
  bestEffort: true
});
print('NIR histogram + stats', nirHistogram);
print(Chart.image.histogram(landsatComposite.select(HISTOGRAM_BANDS), STUDY_AREA, 30));

// ------------------------------------------------------------------
// 5. Otsu's method: maximise between-class variance over the histogram
// ------------------------------------------------------------------
// Given a single-band histogram, find the threshold t that maximises the
// between-class variance BSS(t) = w0*w1*(mu0 - mu1)^2. Equivalently we
// compute, for every candidate split index i:
//    BSS(i) = N0(i)*(mu0(i)-mu)^2 + N1(i)*(mu1(i)-mu)^2
// and return the bucket mean at argmax(BSS).
var computeOtsuThreshold = function (histogramDictionary) {
  var counts = ee.Array(ee.Dictionary(histogramDictionary).get('histogram'));
  var bucketMeans = ee.Array(ee.Dictionary(histogramDictionary).get('bucketMeans'));
  var bucketCount = bucketMeans.length().get([0]);
  var totalCount = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
  var weightedSum = bucketMeans.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
  var grandMean = weightedSum.divide(totalCount);

  var splitIndices = ee.List.sequence(1, bucketCount);

  // For each candidate split i, compute between-class sum of squares (BSS).
  var betweenClassVariance = splitIndices.map(function (i) {
    var leftCounts = counts.slice(0, 0, i);
    var leftCount = leftCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var leftMeans = bucketMeans.slice(0, 0, i);
    var leftMean = leftMeans.multiply(leftCounts)
      .reduce(ee.Reducer.sum(), [0]).get([0])
      .divide(leftCount);
    var rightCount = totalCount.subtract(leftCount);
    var rightMean = weightedSum.subtract(leftCount.multiply(leftMean)).divide(rightCount);
    return leftCount.multiply(leftMean.subtract(grandMean).pow(2))
      .add(rightCount.multiply(rightMean.subtract(grandMean).pow(2)));
  });

  // Helpful diagnostic chart so the user can sanity-check the chosen split.
  print(
    'Between-class variance vs threshold',
    ui.Chart.array.values(ee.Array(betweenClassVariance), 0, bucketMeans)
  );

  // The bucket mean corresponding to the highest BSS is the Otsu threshold.
  return bucketMeans.sort(betweenClassVariance).get([-1]);
};

var otsuThreshold = computeOtsuThreshold(nirHistogram.get(HISTOGRAM_BANDS[0] + '_histogram'));
print('Otsu NIR threshold (water if value <)', otsuThreshold);

// ------------------------------------------------------------------
// 6. Apply the threshold and display the classified water mask
// ------------------------------------------------------------------
var waterMask = landsatComposite.select(HISTOGRAM_BANDS[0]).lt(otsuThreshold);
Map.addLayer(waterMask.mask(waterMask), { palette: 'blue' }, 'Water (Otsu)');
