/*
 * Landsat thermal band -> Land Surface Temperature (LST) with NDVI emissivity
 *
 * Datasets:    LANDSAT/LC08/C01/T1_SR (Surface Reflectance, including thermal B10)
 * Region:      Any user-defined study-area polygon (default: a generic urban AOI).
 * Output:      (1) A median annual composite of NDVI, fractional vegetation cover (FV),
 *                  emissivity and LST (degC).
 *              (2) Map layers of thermal brightness, emissivity and LST.
 *
 * Algorithm:
 *   Brightness temperature (Tb)         = B10 * 0.1                 (scale factor)
 *   Fractional vegetation cover (FV)    = ((NDVI - NDVI_min) /
 *                                          (NDVI_max - NDVI_min))^2
 *   Emissivity (eps)                    = 0.004 * FV + 0.986        (Sobrino-style)
 *   LST [degC] = ( Tb / (1 + (0.00115 * Tb / 1.438) * ln(eps)) ) - 273.15
 *
 * Usage:
 *   1. Paste this script into the GEE Code Editor.
 *   2. Set STUDY_AREA, START_DATE and END_DATE.
 *   3. Click Run; LST map layer appears centered on the AOI.
 */

// ------------------------------------------------------------------
// 0. User-editable constants
// ------------------------------------------------------------------
var STUDY_AREA = ee.Geometry.Polygon(
  [[[44.37581102848051, 38.69236283768664],
    [44.37581102848051, 36.45044284416454],
    [46.89168016910551, 36.45044284416454],
    [46.89168016910551, 38.69236283768664]]], null, false
);

var START_DATE      = '2018-01-01';
var END_DATE        = '2018-12-31';
var THERMAL_SCALE_FACTOR = 0.1;      // Landsat 8 C1 SR thermal scale (B10 * 0.1 -> Kelvin)
var EMISSIVITY_A    = 0.004;          // Sobrino et al. (2008) coefficients
var EMISSIVITY_B    = 0.986;
var REDUCTION_SCALE_M = 30;           // native Landsat optical resolution

// ------------------------------------------------------------------
// 1. Cloud / cloud-shadow mask using the SR pixel_qa band
// ------------------------------------------------------------------
// Landsat 8 SR pixel_qa: bit 3 = cloud shadow, bit 5 = cloud.
// Both bits must be zero for a clear-sky pixel.
function maskLandsat8Sr(image) {
  var cloudShadowBitMask = 1 << 3;
  var cloudsBitMask      = 1 << 5;
  var qa = image.select('pixel_qa');
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
              .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask);
}

// ------------------------------------------------------------------
// 2. Load Landsat 8 SR, filter, mask and composite
// ------------------------------------------------------------------
var landsatCollection = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
  .filterDate(START_DATE, END_DATE)
  .filterBounds(STUDY_AREA)
  .map(maskLandsat8Sr);

// Median composite suppresses residual cloud/haze that survived the QA mask.
var landsatComposite = landsatCollection.median();

// ------------------------------------------------------------------
// 3. NDVI
// ------------------------------------------------------------------
// NDVI = (NIR - RED) / (NIR + RED) = (B5 - B4) / (B5 + B4)
var ndvi = landsatComposite.normalizedDifference(['B5', 'B4']).rename('NDVI');

// ------------------------------------------------------------------
// 4. Thermal brightness temperature (Kelvin)
// ------------------------------------------------------------------
// B10 is delivered in 10*K (DN). Multiplying by 0.1 yields brightness
// temperature in Kelvin. We do NOT subtract 273.15 here because the LST
// equation below operates on Kelvin.
var thermalBrightnessKelvin = landsatComposite
  .select('B10')
  .multiply(THERMAL_SCALE_FACTOR)
  .rename('Tb');

// ------------------------------------------------------------------
// 5. Fractional vegetation cover (FV)
// ------------------------------------------------------------------
// FV uses scene-wide NDVI min/max so that the linear NDVI mixing model is
// stretched between the bare-soil and fully-vegetated end-members observed
// in the AOI, rather than over the global theoretical [-1, 1] range.
var ndviStats = ndvi.reduceRegion({
  reducer: ee.Reducer.min().combine({reducer2: ee.Reducer.max(), sharedInputs: true}),
  geometry: STUDY_AREA,
  scale: REDUCTION_SCALE_M,
  maxPixels: 1e9
});

var ndviMin = ee.Number(ndviStats.get('NDVI_min'));
var ndviMax = ee.Number(ndviStats.get('NDVI_max'));

// FV = ( (NDVI - NDVI_min) / (NDVI_max - NDVI_min) )^2  (squared per Carlson & Ripley 1997).
var fractionalVegCover = ndvi.subtract(ndviMin)
  .divide(ndviMax.subtract(ndviMin))
  .pow(2)
  .rename('FV');

// ------------------------------------------------------------------
// 6. Emissivity (NDVI-threshold method, Sobrino-style)
// ------------------------------------------------------------------
// eps = a * FV + b. The small coefficient pair (0.004, 0.986) yields the
// usual physical range ~ [0.986, 0.990] for vegetated / mixed surfaces.
var emissivity = fractionalVegCover
  .multiply(EMISSIVITY_A)
  .add(EMISSIVITY_B)
  .rename('EMM');

// ------------------------------------------------------------------
// 7. Land Surface Temperature
// ------------------------------------------------------------------
// LST = Tb / (1 + (lambda * Tb / rho) * ln(eps)) - 273.15
// where lambda = 10.8 um (B10 effective wavelength) and
// rho = h*c/sigma_B ~ 1.438e-2 m K. The coefficient 0.00115 collapses
// the wavelength term so the equation can be written compactly.
var lstCelsius = thermalBrightnessKelvin.expression(
  '(Tb / (1 + (0.00115 * (Tb / 1.438)) * log(Ep))) - 273.15',
  {
    Tb: thermalBrightnessKelvin.select('Tb'),
    Ep: emissivity.select('EMM')
  }
).rename('LST');

// ------------------------------------------------------------------
// 8. Visualize results
// ------------------------------------------------------------------
Map.centerObject(STUDY_AREA, 8);

var rgbVis = {bands: ['B4', 'B3', 'B2'], min: 0, max: 3000, gamma: 1.4};
Map.addLayer(landsatComposite, rgbVis, 'Landsat 8 RGB (annual median)');

Map.addLayer(ndvi,
  {min: -0.2, max: 0.8, palette: ['blue', 'white', 'green']},
  'NDVI');

Map.addLayer(thermalBrightnessKelvin,
  {min: 290, max: 320, palette: ['blue', 'white', 'red']},
  'Thermal brightness (K)');

Map.addLayer(emissivity,
  {min: 0.986, max: 0.990},
  'Emissivity');

Map.addLayer(lstCelsius, {
  min: 15, max: 45,
  palette: [
    '040274', '0502b8', '0502ff', '307ef3', '30c8e2',
    '3ff38f', 'b5e22e', 'fff705', 'ffb613', 'ff500d', 'a71001'
  ]
}, 'LST (degC)');
