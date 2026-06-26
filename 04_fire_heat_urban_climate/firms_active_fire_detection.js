/*
 * FIRMS active-fire detection with temporal compositing
 *
 * Datasets:    FIRMS (Fire Information for Resource Management System,
 *              MODIS Terra/Aqua thermal-anomaly product, 1 km native).
 * Region:      Any user-supplied study-area FeatureCollection (e.g. a fire
 *              season AOI). Default example is set to a generic AOI placeholder.
 * Output:      (1) Map layers showing maximum brightness temperature (T21) and
 *                  detection confidence over the fire season.
 *              (2) An ImageCollection of daily FIRMS frames clipped to the AOI.
 *              (3) A Drive export of the daily clipped frames as GeoTIFFs.
 *
 * Usage:
 *   1. Paste this script into the GEE Code Editor.
 *   2. Set STUDY_AREA_ASSET, START_DATE and END_DATE at the top.
 *   3. Click Run; map layers appear, the batch export task queues in the
 *      Tasks tab.
 *
 * Notes:
 *   - FIRMS reports 'T21' as brightness temperature in Kelvin of MODIS band 21
 *     (the 4 micron thermal channel) at the active-fire pixel. Higher T21
 *     indicates a hotter / more intense fire pixel.
 *   - 'confidence' is the FIRMS detection-confidence percentage (0-100).
 *     Production hotspot maps usually keep only pixels with confidence >= 50.
 */

// ------------------------------------------------------------------
// 0. User-editable constants
// ------------------------------------------------------------------
var STUDY_AREA_ASSET    = 'users/<YOUR_GEE_USERNAME>/Australiafinal';   // AOI polygon(s)
var START_DATE          = '2019-11-01';
var END_DATE            = '2020-02-29';                        // exclusive
var CONFIDENCE_THRESHOLD = 50;                                 // % (0-100)
var EXPORT_SCALE_M      = 500;                                 // FIRMS native ~1 km, 500 keeps sub-pixel detail
var EXPORT_FOLDER       = 'GEE_FIRMS';
var EXPORT_CRS          = 'EPSG:4326';

// ------------------------------------------------------------------
// 1. Load AOI and the FIRMS collection
// ------------------------------------------------------------------
var studyArea = ee.FeatureCollection(STUDY_AREA_ASSET);

// FIRMS is a daily ImageCollection. We filter by date once and then derive
// two parallel collections - one for the T21 brightness and one for the
// detection confidence - to keep the rest of the pipeline cheap.
var firmsCollection = ee.ImageCollection('FIRMS')
  .filterDate(START_DATE, END_DATE)
  .filterBounds(studyArea);

// ------------------------------------------------------------------
// 2. Clip each daily image to the AOI and carry timestamps through
// ------------------------------------------------------------------
// We preserve 'system:time_start' / 'system:time_end' so downstream charts,
// animations or compositors can still slice the collection by date.
var clipToAoi = function (image) {
  return image
    .clip(studyArea)
    .copyProperties(image, ['system:time_start', 'system:time_end']);
};

var brightnessDailyCollection = firmsCollection
  .select('T21')
  .map(clipToAoi);

var confidenceDailyCollection = firmsCollection
  .select('confidence')
  .map(clipToAoi);

// ------------------------------------------------------------------
// 3. Build seasonal composites
// ------------------------------------------------------------------
// For active-fire detection we want the *hottest* observation per pixel over
// the season (max), not the mean - because mean dilutes short-lived fires.
var seasonalMaxBrightness = brightnessDailyCollection.max().clip(studyArea);

// For confidence we keep the max as well: any pixel that was ever flagged
// high-confidence during the season stays high-confidence in the composite.
var seasonalMaxConfidence = confidenceDailyCollection.max().clip(studyArea);

// High-confidence active-fire mask: a binary layer ready for hotspot mapping,
// burned-area estimation, or fire-frequency aggregation.
var highConfidenceFireMask = seasonalMaxConfidence
  .gte(CONFIDENCE_THRESHOLD)
  .selfMask()
  .rename('high_confidence_fire');

// ------------------------------------------------------------------
// 4. Visualize on the map
// ------------------------------------------------------------------
Map.centerObject(studyArea, 5);
Map.addLayer(studyArea.style({color: 'black', fillColor: '00000000', width: 1}),
             {}, 'Study area');

Map.addLayer(seasonalMaxBrightness, {
  min: 320, max: 400,
  palette: ['yellow', 'orange', 'red', 'darkred']
}, 'Max T21 brightness (K) - season');

Map.addLayer(seasonalMaxConfidence, {
  min: 0, max: 100,
  palette: ['white', 'yellow', 'orange', 'red']
}, 'Max detection confidence (%) - season');

Map.addLayer(highConfidenceFireMask, {
  palette: ['red']
}, 'High-confidence fire pixels');

// ------------------------------------------------------------------
// 5. Batch export of daily clipped frames
// ------------------------------------------------------------------
// Uses the community 'geetools' batch module to queue one Drive task per
// daily image (otherwise GEE forces one manual click per image).
var batch = require('users/fitoprincipe/geetools:batch');

batch.Download.ImageCollection.toDrive(
  brightnessDailyCollection,
  EXPORT_FOLDER,
  {
    scale: EXPORT_SCALE_M,
    region: studyArea.geometry(),
    crs: EXPORT_CRS,
    name: 'FIRMS_T21_{system_date}'
  }
);
