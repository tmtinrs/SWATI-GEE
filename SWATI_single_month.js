// ============================================================================
// Surface Water Availability-Temperature Index (SWATI)
// Google Earth Engine Implementation — Publication Version
// Near-global monthly SWATI, using fixed normalization bounds
// (mean +/- 3*sd outlier thresholds, derived once for 2000-03 -> 2025-12)
// ============================================================================
//
// ── Citation ────────────────────────────────────────────────────────────────
//   If you use this code, please cite:
//     Thai, M.-T. & Liou, Y.-A. Surface Water Availability-Temperature
//     Index (SWATI) for Global Drought Monitoring. IEEE Trans. Geosci.
//     Remote Sens. 63, 1-14, https://doi.org/10.1109/tgrs.2025.3541287
//     (2025).
//
// ── Contact ───────────────────────────────────────────────────────────────
//   Minh-Tin Thai (tmtinrs@gmail.com) and Yuei-An Liou (yueian@csrsr.ncu.edu.tw)
//   Center for Space and Remote Sensing Research, National Central University,
//   Taoyuan, Taiwan
//
// ── SWATI Formula  ───────────────────────────────────────────────────────────
//
//   SWATI = sqrt[ (sqrt(3)/3 - NNDLI)^2 + (sqrt(3)/3 - NNDVI)^2 + NLST^2 ]
//
//   Variable  Description              Normalized range
//   ────────  ───────────────────────  ────────────────
//   NNDLI     Normalized NDLI          [0, sqrt(3)/3]  — water availability
//   NNDVI     Normalized NDVI          [0, sqrt(3)/3]  — vegetation
//   NLST      Normalized LST           [0, sqrt(3)/3]  — temperature
//   sqrt(3)/3 ~= 0.57735 = ideal wet vertex in 3-D Euclidean distance space
//
//   SWATI -> 0 : wet, cool, vegetated  (no dryness)
//   SWATI -> 1 : hot, dry, bare        (extreme dryness)
//
// ── Data Sources ────────────────────────────────────────────────────────────
//   MODIS MOD09GA v061  — Surface Reflectance 500 m daily (Terra)
//   MODIS MOD11A1 v061  — LST 1 km daily (Terra)
//   MODIS MOD44W        — Permanent Water Mask 250 m (static)
//
// ── Normalization approach ────────────────────────────────────────────────────
//   Phase 1 (computed offline, once, over the full study period):
//     For each month 2000-03 -> 2025-12, compute monthly mean composites,
//     then reduceRegion(mean + stdDev). Per-month outlier thresholds:
//       upper = mean + 3*sd
//       lower = mean - 3*sd   (NDVI lower fixed at 0.0)
//     The FIXED dataset-wide bounds are then:
//       high* = max(all monthly uppers),  low* = min(all monthly lowers)
//     These fixed bounds are set below — every month is normalized
//     against the SAME scale, making SWATI values comparable across time.
//   Phase 2 (this script):
//     clamp each monthly composite to [low*, high*], then unitScale -> [0, sqrt(3)/3].
//
// ── Repository ────────────────────────────────────────────────────────────────
//   See ../README.md for the complete pipeline description, fixed-bounds
//   table, usage instructions, and citation.
// ============================================================================

// ── Study area — near-global land boundary used in this study.
//    Replace with your own registered GEE asset (see README, "Fixed
//    normalization bounds", before using a different region). ──────────────
var table = ee.FeatureCollection('projects/ee-hrsl/assets/GlobalBoundary');

// ── Target period — change YEAR and MONTH only ───────────────────────────────
var YEAR  = 2025;
var MONTH = 1;
// ─────────────────────────────────────────────────────────────────────────────

var startDate  = ee.Date.fromYMD(YEAR, MONTH, 1);
var endDate    = startDate.advance(1, 'month');
var monthStr   = (MONTH < 10 ? '0' + MONTH : '' + MONTH);
var exportName = 'SWATI_' + YEAR + monthStr;

// ============================================================================
// SECTION 1 — FIXED NORMALIZATION BOUNDS
//              (near-global study area, 2000-03 -> 2025-12)
//
//   WHAT THESE SIX NUMBERS ARE
//   They were computed ONCE for this study's near-global study area and
//   study period (310 months, 2000-03 -> 2025-12), as follows:
//     1. For every month, the mean and standard deviation (sd) of each
//        variable (NDVI, NDLI, LST) were calculated over the whole study
//        area. The monthly outlier thresholds are mean - 3*sd (lower) and
//        mean + 3*sd (upper).
//     2. The fixed bound of each variable is the most extreme threshold found
//        in any month:  low = min(all monthly lower thresholds),
//                       high = max(all monthly upper thresholds).
//   They are set here as fixed values so that every month is normalized on exactly
//   the same scale, which keeps SWATI values comparable across time.
//
//   NDVI notes:
//     • lowNDVI  = 0.0  (fixed by design, not derived statistically.
//                         NDVI <= 0 pixels, which represent non-vegetated
//                         surfaces such as water and snow/ice, are masked
//                         out before computing the mean/sd used for
//                         highNDVI, and NDVI values below 0 are clamped
//                         to 0 in Section 4.)
//     • highNDVI = max monthly (mean + 3*sd), computed the same way as
//                  NDLI's and LST's upper bounds.
//
//   ┌─────────┬──────────────────────┬──────────────────────┐
//   │ Variable│ lowBound             │ highBound            │
//   ├─────────┼──────────────────────┼──────────────────────┤
//   │ NDVI    │ 0.0000 (fixed)       │ 0.9999               │
//   │ NDLI    │ -0.2408              │ 0.2063               │
//   │ LST °C  │ -64.2579             │ 82.3896              │
//   └─────────┴──────────────────────┴──────────────────────┘
//
//   WHEN TO CHANGE THESE VALUES
//   • They are valid only for THIS study area and THIS period. If you use a
//     different region or a different study period, define the six values
//     again for your own area of interest and period (same steps as above).
//   • Use ONE set of bounds for a whole time series. Never mix SWATI maps
//     produced with different bounds. If you extend the record or change
//     the bounds, re-run every month so the entire archive stays on one scale.
//
//   WHAT TO EXPECT FROM NEAR-GLOBAL BOUNDS
//   • The bounds cover the full range of conditions found across the
//     near-global extent, so SWATI shows a strong contrast between arid /
//     semi-arid regions (high SWATI) and densely vegetated regions (low SWATI).
//   • Within a smaller region, SWATI values fall in a narrow part of the 0-1
//     range. They describe the absolute dryness level of that region (a desert
//     is always "dry"), not whether the region is drier or wetter than usual.
//
//   USING THE NEAR-GLOBAL SWATI FOR A SPECIFIC REGION (see README.md)
//   Option 1 — Rescale SWATI to 0-1 to present wetness/dryness conditions
//              within the region:
//                SWATI_scaled = (SWATI - SWATI_min) / (SWATI_max - SWATI_min)
//              SWATI_min / SWATI_max = lowest / highest SWATI of your region
//              (or its 2nd / 98th percentiles). Use the same two values for
//              every month. 0 = wettest, 1 = driest WITHIN that region; the
//              result is not comparable across regions and the original SWATI
//              class thresholds do not apply to it.
//   Option 2 — Standardize SWATI (recommended for regional analysis): detect
//              dryness/wetness anomalies relative to long-term conditions, per
//              pixel and calendar month m:
//                SWATI_z = (SWATI - mean_m) / sd_m
//              mean_m / sd_m = long-term mean / sd of SWATI for the same
//              calendar month. SWATI_z > 0 = drier than normal,
//              SWATI_z < 0 = wetter than normal.
// ============================================================================

var lowNDVI  = 0.0000;
var highNDVI = 0.9999;

var lowNDLI  = -0.2408;
var highNDLI =  0.2063;

var lowLST   = -64.2579;
var highLST  =  82.3896;

var SQRT3_3  = 0.57735;       // sqrt(3)/3 — ideal wet/cool/vegetated vertex

// ============================================================================
// SECTION 2 — MASKING FUNCTIONS
// ============================================================================

// ── 2a. MOD09GA cloud & shadow mask (state_1km QA) ───────────────────────────
//   Bits 0-1 = 00 -> clear sky
//   Bit  2   =  0 -> no cloud shadow
var applyMOD09GACloudMask = function(img) {
  var qa       = img.select('state_1km');
  var isClear  = qa.bitwiseAnd(3).eq(0);
  var noShadow = qa.bitwiseAnd(4).eq(0);
  return img.updateMask(isClear.and(noShadow));
};

// ── 2b. MOD11A1 LST quality mask (QC_Day) ────────────────────────────────────
//   Bits 0-1 = 00 (good quality) or 01 (nominal quality) -> retain
var applyLSTQC = function(img) {
  var qc   = img.select('QC_Day');
  var good = qc.bitwiseAnd(3).lte(1);
  return img.updateMask(good);
};

// ============================================================================
// SECTION 3 — LOAD & COMPOSITE INPUT VARIABLES
//   Spatial resolution harmonisation:
//     MOD09GA reflectance bands (500 m native) are resampled to 1 km using
//     bilinear interpolation via .resample('bilinear'). The actual pixel
//     grid is resolved at 1 km during export (scale=1000, crs='EPSG:4326'),
//     matching the MOD11A1 LST native resolution.
//     No explicit .reproject() is applied inside the map function, which
//     would force premature server-side materialisation and exceed GEE
//     memory limits on a near-global study area.
// ============================================================================

// ── 3a. MOD09GA — single cloud-masked collection (shared by NDVI & NDLI) ─────
var modisReflectance = ee.ImageCollection('MODIS/061/MOD09GA')
  .filter(ee.Filter.date(startDate, endDate))
  .map(applyMOD09GACloudMask);

// ── 3b. NDVI — Normalized Difference Vegetation Index ────────────────────────
//   NDVI = (NIR - Red) / (NIR + Red)
//   MOD09GA: NIR = Band 2, Red = Band 1
//   Resampled to 1 km (bilinear) to match LST grid at export scale.
var meanNDVI = modisReflectance
  .map(function(img) {
    var img_rs = img.resample('bilinear');   // bilinear; resolved at scale=1000
    var nir = img_rs.select('sur_refl_b02');
    var red = img_rs.select('sur_refl_b01');
    return nir.subtract(red)
              .divide(nir.add(red))
              .rename('NDVI');
  })
  .mean()
  .clip(table);

// ── 3c. NDLI — Normalized Difference Latent heat Index ───────────────────────
//   NDLI = (Green - Red) / (Green + Red + SWIR)
//   MOD09GA: Green = Band 4, Red = Band 1, SWIR = Band 6
//   Resampled to 1 km (bilinear) to match LST grid at export scale.
var meanNDLI = modisReflectance
  .map(function(img) {
    var img_rs = img.resample('bilinear');
    var swir  = img_rs.select('sur_refl_b06');
    var red   = img_rs.select('sur_refl_b01');
    var green = img_rs.select('sur_refl_b04');
    return green.subtract(red)
                .divide(green.add(red).add(swir))
                .rename('NDLI');
  })
  .mean()
  .clip(table);

// ── 3d. LST — Land Surface Temperature (MOD11A1, converted to °C) ────────────
//   Scale factor: 0.02 K/DN; offset: -273.15 K -> °C
var meanLST = ee.ImageCollection('MODIS/061/MOD11A1')
  .filter(ee.Filter.date(startDate, endDate))
  .map(applyLSTQC)
  .map(function(img) {
    return img.select('LST_Day_1km')
              .multiply(0.02)
              .subtract(273.15)
              .rename('LST');
  })
  .mean()
  .clip(table);

// ============================================================================
// SECTION 4 — CLAMP TO FIXED BOUNDS, NORMALIZE, COMPUTE SWATI
//   Using the fixed dataset-wide mu +/- 3*sd bounds defined above (not per-month
//   percentiles). This ensures temporal comparability across all 310 months.
// ============================================================================

// ── 4a. Clamp to remove extreme outliers ─────────────────────────────────────
meanLST  = meanLST.clamp(lowLST,   highLST);
meanNDVI = meanNDVI.clamp(lowNDVI, highNDVI);
meanNDLI = meanNDLI.clamp(lowNDLI, highNDLI);

// ── 4b. Normalize each variable to [0, sqrt(3)/3] ────────────────────────────
//   unitScale(lo, hi) maps [lo, hi] -> [0, 1]; multiply by sqrt(3)/3 -> [0, sqrt(3)/3]
var normalizedLST  = meanLST.unitScale(lowLST,   highLST).multiply(SQRT3_3);
var normalizedNDVI = meanNDVI.unitScale(lowNDVI, highNDVI).multiply(SQRT3_3);
var normalizedNDLI = meanNDLI.unitScale(lowNDLI, highNDLI).multiply(SQRT3_3);

// ── 4c. SWATI = sqrt[ (sqrt(3)/3 - NNDLI)^2 + (sqrt(3)/3 - NNDVI)^2 + NLST^2 ] ─
var swati = normalizedNDLI.subtract(SQRT3_3).pow(2)
  .add(normalizedNDVI.subtract(SQRT3_3).pow(2))
  .add(normalizedLST.pow(2))
  .sqrt()
  .rename('SWATI');

// ============================================================================
// SECTION 5 — WATER MASK (MOD44W — land pixels only)
// ============================================================================

var waterMask = ee.ImageCollection('MODIS/006/MOD44W')
  .filter(ee.Filter.date('2000-01-01', '2015-01-01'))
  .select('water_mask')
  .mosaic()
  .eq(0);   // 0 = land; 1 = permanent water -> retain land only

var swatiMasked = swati.clip(table).updateMask(waterMask);

// ============================================================================
// SECTION 6 — DISPLAY
// ============================================================================

var swatiParams = {
  min: 0, max: 1,
  palette: ['040274','040281','0502a3','0502b8','0502ce','0502e6',
            '0602ff','235cb1','307ef3','269db1','30c8e2','32d3ef',
            '3be285','3ff38f','86e26f','3ae237','b5e22e','d6e21f',
            'fff705','ffd611','ffb613','ff8b13','ff6e08','ff500d',
            'ff0000','de0101','c21301','a71001','911003']
};

Map.addLayer(swatiMasked, swatiParams, 'SWATI ' + YEAR + '-' + monthStr);
Map.centerObject(table, 2);

// ── Legend ───────────────────────────────────────────────────────────────────
var legend = ui.Panel({ style: { position: 'bottom-left', padding: '8px 15px' } });

legend.add(ui.Label({
  value: 'SWATI',
  style: { fontWeight: 'bold', fontSize: '18px', margin: '0 0 4px 0', padding: '0' }
}));

legend.add(ui.Panel({ widgets: [ui.Label('' + swatiParams.max)] }));

var legendGradient = ee.Image.pixelLonLat().select('latitude')
  .multiply((swatiParams.max - swatiParams.min) / 100.0)
  .add(swatiParams.min);

legend.add(ui.Thumbnail({
  image: legendGradient.visualize(swatiParams),
  params: { bbox: '0,0,15,100', dimensions: '15x100' },
  style: { padding: '1px', position: 'bottom-center' }
}));

legend.add(ui.Panel({ widgets: [ui.Label('' + swatiParams.min)] }));
Map.add(legend);

// ============================================================================
// SECTION 7 — EXPORT TO GOOGLE DRIVE
// ============================================================================

Export.image.toDrive({
  image       : swatiMasked,
  description : exportName,
  folder      : 'SWATI_results',
  region      : table.geometry(),
  scale       : 1000,
  crs         : 'EPSG:4326',
  maxPixels   : 1e13
});
