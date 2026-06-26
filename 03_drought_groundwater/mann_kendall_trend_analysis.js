/*
 * Mann-Kendall non-parametric trend test with Sen-slope and tie-corrected variance
 *
 * Datasets:    MODIS/006/MOD13A1 (EVI / NDVI, 500 m, 16-day) by default
 *              -- the test is generic and accepts any per-pixel annual or
 *                 seasonal image collection (precipitation, ET, RUE, etc.).
 * Region:      MENA / global; any FeatureCollection works
 * Output:      Per-pixel Mann-Kendall S statistic, Sen-slope, Z-score,
 *              p-value, and a binary mask of statistically significant trends
 *              (alpha = 0.05 two-sided).
 *
 * Method:      For every pair of observations (i, j) with t_j > t_i:
 *                S       = sum( sign(x_j - x_i) )
 *                slope   = median over all pairs of (x_j - x_i) / (t_j - t_i)
 *                Var(S)  = [ n(n-1)(2n+5) - sum_g g(g-1)(2g+5) ] / 18
 *                          where the second term corrects for ties in groups
 *                          of equal observations.
 *                Z       = (S - sign(S)) / sqrt(Var(S))
 *                p-value = 1 - Phi(|Z|)
 *              See Mann (1945), Kendall (1975), Sen (1968).
 *
 * Usage:
 *   1. Paste into the GEE Code Editor.
 *   2. Swap INPUT_COLLECTION / INPUT_BAND for the series you want tested.
 *   3. Click Run; trend rasters appear on the map.
 */

// ---------------------------------------------------------------------------
// USER CONFIGURATION
// ---------------------------------------------------------------------------

var AOI               = ee.FeatureCollection('users/<YOUR_GEE_USERNAME>/MENA');

// Source collection -- swap for any per-pixel annual/seasonal stack.
var INPUT_COLLECTION  = ee.ImageCollection('MODIS/006/MOD13A1');
var INPUT_BAND        = 'EVI';

// Optional seasonal filter: keeps only images falling inside this calendar
// range (e.g. months 8-9 = late summer). Set START_MONTH = 1 and END_MONTH = 12
// for the full year.
var START_MONTH       = 8;
var END_MONTH         = 9;

var SIGNIFICANCE_ALPHA = 0.05;       // two-sided
var TREND_PALETTE      = ['red', 'white', 'green'];
var PARALLEL_SCALE     = 2;          // raise if you hit memory limits

// ---------------------------------------------------------------------------
// LOAD AND CLIP THE INPUT TIME SERIES
// ---------------------------------------------------------------------------

var inputCollection = INPUT_COLLECTION
  .select(INPUT_BAND)
  .filter(ee.Filter.calendarRange(START_MONTH, END_MONTH, 'month'))
  .map(function (img) { return img.clip(AOI); });

Map.centerObject(AOI, 5);
Map.addLayer(inputCollection.mean(), {}, 'Mean of input series', false);

// ---------------------------------------------------------------------------
// JOIN EACH IMAGE WITH ALL LATER IMAGES (the (i, j) pair space)
// ---------------------------------------------------------------------------

var afterFilter = ee.Filter.lessThan({
  leftField:  'system:time_start',
  rightField: 'system:time_start'
});

var joined = ee.ImageCollection(
  ee.Join.saveAll('after').apply({
    primary:   inputCollection,
    secondary: inputCollection,
    condition: afterFilter
  })
);

// ---------------------------------------------------------------------------
// MANN-KENDALL S STATISTIC: SUM OF SIGNS OVER ALL PAIRS
// ---------------------------------------------------------------------------

function sign (imageI, imageJ) {
  return ee.Image(imageJ).neq(imageI)
    .multiply(ee.Image(imageJ).subtract(imageI).clamp(-1, 1))
    .int();
}

var kendallStatistic = ee.ImageCollection(joined.map(function (current) {
  var laterImages = ee.ImageCollection.fromImages(current.get('after'));
  return laterImages.map(function (image) {
    // unmask(0) keeps the sum well-defined when either side is masked.
    return ee.Image(sign(current, image)).unmask(0);
  });
}).flatten()).reduce('sum', PARALLEL_SCALE);

Map.addLayer(kendallStatistic, { palette: TREND_PALETTE }, 'Mann-Kendall S');

// ---------------------------------------------------------------------------
// SEN SLOPE: MEDIAN PAIRWISE SLOPE
// ---------------------------------------------------------------------------

function pairwiseSlope (imageI, imageJ) {
  return ee.Image(imageJ).subtract(imageI)
    .divide(ee.Image(imageJ).date().difference(ee.Image(imageI).date(), 'days'))
    .rename('slope')
    .float();
}

var allSlopes = ee.ImageCollection(joined.map(function (current) {
  var laterImages = ee.ImageCollection.fromImages(current.get('after'));
  return laterImages.map(function (image) {
    return ee.Image(pairwiseSlope(current, image));
  });
}).flatten());

var sensSlope = allSlopes.reduce(ee.Reducer.median(), PARALLEL_SCALE);
Map.addLayer(sensSlope, { palette: TREND_PALETTE }, 'Sen slope (per day)');

// ---------------------------------------------------------------------------
// TIE CORRECTION FOR THE VARIANCE OF S
// ---------------------------------------------------------------------------
// If g is the size of a group of equal observations, that group reduces the
// variance of S by g(g-1)(2g+5)/18. We have to detect tied values per pixel
// and accumulate the correction.

// Mark values that appear in a tie group (else zero).
var tiedValues = inputCollection.map(function (i) {
  var matches = inputCollection.map(function (j) {
    return i.eq(j);
  }).sum();
  return i.multiply(matches.gt(1));
});

// Given an array of values per pixel, return the lengths of equal-value runs.
function groupSizes (array) {
  var length  = array.arrayLength(0);
  var indices = ee.Image([1])
    .arrayRepeat(0, length)
    .arrayAccum(0, ee.Reducer.sum())
    .toArray(1);
  var sorted = array.arraySort();
  var left   = sorted.arraySlice(0, 1);
  var right  = sorted.arraySlice(0, 0, -1);
  // 1 marks the end of every run; pad the last position with 1.
  var mask = left.neq(right).arrayCat(ee.Image(ee.Array([[1]])), 0);
  var runEnds = indices.arrayMask(mask);
  return runEnds.arraySlice(0, 1).subtract(runEnds.arraySlice(0, 0, -1));
}

// Term f(n) = n(n-1)(2n+5) from Kendall's variance formula.
function varianceTerm (image) {
  return image.expression('b() * (b() - 1) * (b() * 2 + 5)');
}

var tieGroupSizes  = groupSizes(tiedValues.toArray());
var tieCorrection  = varianceTerm(tieGroupSizes).arrayReduce('sum', [0]).arrayGet([0, 0]);
var imageCount     = joined.count();

var kendallVariance = varianceTerm(imageCount)
  .subtract(tieCorrection)
  .divide(18)
  .float();

Map.addLayer(kendallVariance, {}, 'Var(S)', false);

// ---------------------------------------------------------------------------
// Z-STATISTIC AND P-VALUE
// ---------------------------------------------------------------------------
// Z = (S - sign(S)) / sqrt(Var(S))
// The piecewise construction below handles S > 0, S < 0, and S = 0 separately
// in pure ee.Image arithmetic.

var zeroPart = kendallStatistic.multiply(kendallStatistic.eq(0));
var positive = kendallStatistic.multiply(kendallStatistic.gt(0)).subtract(1);
var negative = kendallStatistic.multiply(kendallStatistic.lt(0)).add(1);

var zStatistic = zeroPart
  .add(positive.divide(kendallVariance.sqrt()))
  .add(negative.divide(kendallVariance.sqrt()));

Map.addLayer(zStatistic, { min: -2, max: 2 }, 'Z statistic');

// Standard-normal CDF using the error function: Phi(z) = 0.5*(1 + erf(z/sqrt(2)))
function standardNormalCdf (z) {
  return ee.Image(0.5).multiply(
    ee.Image(1).add(ee.Image(z).divide(ee.Image(2).sqrt()).erf())
  );
}

var pValue = ee.Image(1).subtract(standardNormalCdf(zStatistic.abs()));
Map.addLayer(pValue, { min: 0, max: 1 }, 'p-value', false);

// Reject the null hypothesis of no trend when p <= alpha/2 (two-sided).
var significantTrends = pValue.lte(SIGNIFICANCE_ALPHA / 2);
Map.addLayer(
  significantTrends,
  { min: 0, max: 1, palette: ['white', 'black'] },
  'Significant trends (alpha=' + SIGNIFICANCE_ALPHA + ')'
);
