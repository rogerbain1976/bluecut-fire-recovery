// =====================================================================
// Blue Cut Fire Vegetation Recovery: Rebuilt Analysis (Google Earth Engine)
// Author: Roger Bain
// Paste into the Earth Engine Code Editor (code.earthengine.google.com)
// and click Run. Exports appear in the Tasks tab.
// =====================================================================

// ---------------------------------------------------------------------
// 1. FIRE PERIMETER (MTBS)
// Filtered by name AND location near Cajon Pass so we get the right fire.
// CHECK: the Console should print exactly 1 feature.
// ---------------------------------------------------------------------
var cajonPass = ee.Geometry.Point([-117.45, 34.29]);
var fire = ee.FeatureCollection('USFS/GTAC/MTBS/burned_area_boundaries/v1')
  .filter(ee.Filter.stringContains('Incid_Name', 'BLUE CUT'))
  .filterBounds(cajonPass.buffer(20000));
print('Fire perimeter (should be 1 feature):', fire);
var region = fire.geometry();
print('Perimeter area (acres):', region.area(1).divide(4046.86));

// ---------------------------------------------------------------------
// 2. LANDSAT 8/9 COLLECTION 2 SURFACE REFLECTANCE, PROPERLY SCALED
// This is the fix for the original project: raw SR values must be
// multiplied by 0.0000275 and offset by -0.2 before computing indices.
// Clouds, cloud shadow, cirrus, snow, and saturated pixels are masked.
// ---------------------------------------------------------------------
function prepLandsat(img) {
  var sr = img.select('SR_B.').multiply(0.0000275).add(-0.2);
  var qa = img.select('QA_PIXEL');
  // Bits 0-5: fill, dilated cloud, cirrus, cloud, cloud shadow, snow
  var clear = qa.bitwiseAnd(63).eq(0);
  var unsaturated = img.select('QA_RADSAT').eq(0);
  return img.addBands(sr, null, true)
    .updateMask(clear)
    .updateMask(unsaturated);
}

var landsat = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
  .filterBounds(region)
  .map(prepLandsat);

// ---------------------------------------------------------------------
// 3. TIME PERIODS
// Recovery comparisons use SUMMER composites so every period is compared
// in the same season as the pre-fire baseline. Month6 is a winter
// composite: use it for the trajectory chart only, not recovery percent.
// The fire started Aug 16, 2016 and was contained around Aug 23.
// ---------------------------------------------------------------------
var periods = [
  {name: 'Prefire',  start: '2016-06-01', end: '2016-08-15'},
  {name: 'Postfire', start: '2016-08-24', end: '2016-10-15'},
  {name: 'Month6',   start: '2017-01-15', end: '2017-03-15'},
  {name: 'Year1',    start: '2017-06-01', end: '2017-08-31'},
  {name: 'Year5',    start: '2021-06-01', end: '2021-08-31'},
  {name: 'Year9',    start: '2025-06-01', end: '2025-08-31'},
  {name: 'Year10',   start: '2026-06-01', end: '2026-08-31'}
];

// Median composite per period, then NDVI and NBR
var indexBands = [];
periods.forEach(function(p) {
  var col = landsat.filterDate(p.start, p.end);
  print(p.name + ' scenes used:', col.size());
  var comp = col.median();
  var ndvi = comp.normalizedDifference(['SR_B5', 'SR_B4']).rename(p.name + '_NDVI');
  var nbr  = comp.normalizedDifference(['SR_B5', 'SR_B7']).rename(p.name + '_NBR');
  indexBands.push(ndvi, nbr);
});
var indices = ee.Image.cat(indexBands);

// dNBR from properly scaled data (standard severity thresholds now apply)
var dnbr = indices.select('Prefire_NBR')
  .subtract(indices.select('Postfire_NBR'))
  .rename('dNBR');

// ---------------------------------------------------------------------
// 4. INDEPENDENT AND COVARIATE DATA
// MTBS severity: independent check on your dNBR
//   (1 unburned/low, 2 low, 3 moderate, 4 high, 5 increased greenness, 6 no data)
// NLCD 2016: pre-fire land cover, so developed areas and water can be
//   filtered out and vegetation type can be tested as a recovery factor.
// 3DEP 10 m DEM: continuous elevation, slope, and aspect.
// ---------------------------------------------------------------------
var mtbs = ee.ImageCollection('USFS/GTAC/MTBS/annual_burn_severity_mosaics/v1')
  .filterDate('2016-01-01', '2016-12-31')
  .filter(ee.Filter.stringContains('system:index', 'CONUS'))
  .first()
  .select('Severity')
  .rename('MTBS_Severity');

var nlcd = ee.ImageCollection('USGS/NLCD_RELEASES/2019_REL/NLCD')
  .filter(ee.Filter.eq('system:index', '2016'))
  .first()
  .select('landcover')
  .rename('NLCD2016');

var dem = ee.Image('USGS/3DEP/10m').select('elevation');
var terrain = ee.Image.cat([
  dem.rename('Elevation_m'),
  ee.Terrain.slope(dem).rename('Slope_deg'),
  ee.Terrain.aspect(dem).rename('Aspect_deg')
]);

// Full analysis stack. Masked values (clouds, gaps) become -9999
// so every point is kept and missing values can be handled in Python.
var stack = ee.Image.cat([indices, dnbr, mtbs, nlcd, terrain])
  .clip(region)
  .unmask(-9999);

// ---------------------------------------------------------------------
// 5. RANDOM SAMPLE POINTS
// Fixed seed makes the sample reproducible. Coordinates are kept so
// spatial autocorrelation can be tested.
// ---------------------------------------------------------------------
var points = ee.FeatureCollection.randomPoints({
  region: region,
  points: 10000,
  seed: 2016,
  maxError: 1
}).map(function(f) {
  var c = f.geometry().coordinates();
  return f.set({lon: c.get(0), lat: c.get(1)});
});

var samples = stack.sampleRegions({
  collection: points,
  scale: 30,
  tileScale: 4,
  geometries: false
});

// ---------------------------------------------------------------------
// 6. EXPORTS (run each from the Tasks tab)
// ---------------------------------------------------------------------
Export.table.toDrive({
  collection: samples,
  description: 'BlueCut_Recovery_Points',
  folder: 'BlueCut_GEE',
  fileFormat: 'CSV'
});

// GeoTIFFs for StoryMap maps, UTM Zone 11N, 30 m
periods.forEach(function(p) {
  Export.image.toDrive({
    image: indices.select(p.name + '_NDVI').clip(region).toFloat(),
    description: 'BlueCut_NDVI_' + p.name,
    folder: 'BlueCut_GEE',
    region: region,
    scale: 30,
    crs: 'EPSG:32611',
    maxPixels: 1e10
  });
});
Export.image.toDrive({
  image: dnbr.clip(region).toFloat(),
  description: 'BlueCut_dNBR',
  folder: 'BlueCut_GEE',
  region: region,
  scale: 30,
  crs: 'EPSG:32611',
  maxPixels: 1e10
});

// ---------------------------------------------------------------------
// 7. MAP DISPLAY (same color stretch for every period, so screenshots
// and swipe comparisons are honest)
// ---------------------------------------------------------------------
var ndviVis = {min: 0, max: 0.6,
  palette: ['8c510a', 'd8b365', 'f6e8c3', 'c7e9c0', '74c476', '238b45', '00441b']};
var dnbrVis = {min: -0.1, max: 0.8,
  palette: ['1a9850', 'd9ef8b', 'fee08b', 'fc8d59', 'd73027', '67001f']};

Map.centerObject(region, 11);
Map.setOptions('SATELLITE');
periods.forEach(function(p) {
  Map.addLayer(indices.select(p.name + '_NDVI').clip(region), ndviVis,
    p.name + ' NDVI', p.name === 'Prefire');
});
Map.addLayer(dnbr.clip(region), dnbrVis, 'dNBR (burn severity)', false);
Map.addLayer(ee.Image().paint(fire, 0, 2), {palette: 'ffffff'}, 'Fire perimeter');
