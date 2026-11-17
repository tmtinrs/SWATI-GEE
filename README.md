# SWATI — Surface Water Availability-Temperature Index (Google Earth Engine)

Google Earth Engine implementation of the **Surface Water Availability-Temperature
Index (SWATI)**, a 3‑D Euclidean-distance dryness/wetness (DW) index built
from MODIS NDVI, NDLI, and LST. The code produces **near-global** monthly SWATI
maps at 1 km resolution. This repository accompanies:

> Thai, M.-T. & Liou, Y.-A. Surface Water Availability-Temperature Index (SWATI)
> for Global Drought Monitoring. *IEEE Trans. Geosci. Remote Sens.* 63, 1–14,
> https://doi.org/10.1109/tgrs.2025.3541287 (2025).

If you use this code, please cite the paper above (see
[`CITATION.cff`](CITATION.cff)) and the MODIS data providers acknowledged below.

---

## Contents

```
SWATI-GEE/
├── README.md
├── LICENSE
├── CITATION.cff
└── code/
    └── SWATI_single_month.js   <- computes SWATI for one target month
```

## Method summary

SWATI integrates three MODIS-derived variables into a single dryness/wetness
index using Euclidean distance theory in normalized 3-D space:

```
SWATI = sqrt[ (√3/3 − NNDLI)² + (√3/3 − NNDVI)² + NLST² ]
```

where `NNDLI`, `NNDVI`, and `NLST` are NDLI, NDVI, and LST each normalized to
`[0, √3/3]`. SWATI ranges from 0 (wet, cool, vegetated) to 1 (hot, dry, bare).

**Inputs**:
| Variable | Source | Native resolution | Notes |
|---|---|---|---|
| NDLI, NDVI | MODIS MOD09GA.061 surface reflectance | 500 m | cloud/shadow-masked via `state_1km` QA, resampled to 1 km (bilinear) |
| LST | MODIS MOD11A1.061 | 1 km | quality-masked via `QC_Day` |
| Water mask | MODIS MOD44W (2000–2015 static) | 250 m | permanent water bodies excluded from output |

**Pipeline** (see the paper's Methods section for full detail):
1. **Monthly compositing** — daily images are cloud/QC-masked individually,
   then averaged per calendar month (per-pixel arithmetic mean).
2. **Outlier removal** — for each month, the mean and standard deviation (SD)
   of each variable are calculated over the study area, and the outlier
   thresholds are set to mean ± 3 SD. NDVI's lower bound is fixed at 0, and its
   mean and SD are computed from NDVI > 0 pixels only. The monthly thresholds
   are then collapsed into a single **fixed bound pair per variable** (max of
   all monthly uppers, min of all monthly lowers). These fixed values are
   defined in the script (`lowNDVI`, `highNDVI`, `lowNDLI`, `highNDLI`,
   `lowLST`, `highLST`), so every month is normalized on an identical scale.
3. **Normalization** — each monthly composite is clamped to the fixed
   bounds and min–max scaled to `[0, √3/3]`.
4. **SWATI computation** — normalized NDLI, NDVI, and LST are combined via
   the Euclidean-distance formula above; static water bodies are masked
   out of the final SWATI output.

## Fixed normalization bounds

**Values used in this study** (near-global study area, 2000-03 to 2025-12):

| Variable | Lower bound | Upper bound |
|---|---|---|
| NDLI | −0.2408 | 0.2063 |
| NDVI | 0.0000 (fixed) | 0.9999 |
| LST (°C) | −64.2579 | 82.3896 |

**What these values are.** They were calculated once, for this study's
near-global study area and for the period March 2000 to December 2025 (310
months). Each value is the most extreme mean ± 3 SD threshold found in any
month of that period. They are fixed so that a given NDLI, NDVI, or LST value
always maps to the same normalized value, whichever month it occurs in.

**Use them only for this study area and period.** If you work with another
region or another study period, define the six bounds again for your own area
of interest and your own period. Use one set of bounds for a whole time series,
and never mix SWATI maps produced with different bounds.

**What to expect from near-global bounds.** Because the bounds cover the whole
near-global extent, they span a very wide range of conditions (for example, LST
from about −64 °C to 82 °C). As a result:

- SWATI shows **strong contrast** between arid/semi-arid regions (high SWATI,
  dry) and densely vegetated regions (low SWATI, wet).
- Within a single smaller region, SWATI values fall in a **narrow part of the
  0–1 range**. They describe the **absolute dryness level** of that region (a
  desert is always "dry"), not whether the region is drier or wetter than usual.

### Using the near-global SWATI for a specific region

Because of this, there are two ways to work with the near-global SWATI maps
for one region. Choose the one that matches what you want to show.

**Option 1 — Rescale SWATI to 0–1 (to present wetness and dryness conditions
within a region).** Stretch the SWATI values of your region so that the narrow
range they occupy fills the full 0–1 scale:

```
SWATI_scaled = ( SWATI − SWATI_min ) / ( SWATI_max − SWATI_min )
```

where `SWATI_min` and `SWATI_max` are the lowest and highest SWATI values in
your region (for your period of interest). To reduce the effect of extreme
pixels, you can use the 2nd and 98th percentiles instead. After rescaling,
0 is the wettest and 1 is the driest condition **within that region**.

- Use the same `SWATI_min` and `SWATI_max` for every month of your analysis.
- The result shows the **relative** pattern inside the region only. It cannot be
  compared with other regions, and it is no longer the original SWATI, so do
  not apply the original SWATI dryness class thresholds to it.

**Option 2 — Standardize SWATI (recommended for regional analysis).** If you use
the near-global SWATI maps for a specific region, apply standardization to
detect dryness/wetness anomalies relative to long-term conditions. For each
pixel and each calendar month *m*:

```
SWATI_z(year, m) = ( SWATI(year, m) − mean_m ) / sd_m
```

where `mean_m` and `sd_m` are the long-term mean and standard deviation of SWATI
for that same calendar month over the study period (for example, 2000–2025).
Because SWATI increases with dryness:

- `SWATI_z > 0` means **drier** than normal for that month,
- `SWATI_z < 0` means **wetter** than normal for that month,
- the larger the absolute value, the stronger the anomaly.

In short: use **Option 1** to *display* where a region is relatively wet or dry,
and **Option 2** to *detect* whether a given month or year is drier or wetter
than usual.

## Requirements

- A Google Earth Engine account with a registered cloud project
  ([sign up](https://earthengine.google.com/)).
- A GEE feature-collection asset defining your study area boundary.

## Usage

Open `code/SWATI_single_month.js` in the [GEE Code
Editor](https://code.earthengine.google.com/):

1. Replace the study-area asset path (`table`) with your own registered GEE
   asset.
2. Set `YEAR` and `MONTH` for the month you want to compute.
3. Run. The script displays the SWATI layer on the map with a legend and
   submits an export task (Google Drive, GeoTIFF, EPSG:4326, 1 km) — start
   the task from the **Tasks** tab in the Code Editor.

**Using a different region or period.** The script already contains this
study's bounds. To use your own region or period, recompute the
six values: for every month of your period, calculate the mean and SD of NDVI
(NDVI > 0 pixels only), NDLI, and LST over your study area; take mean − 3 SD
and mean + 3 SD as the monthly thresholds; then take the minimum of all monthly
lower thresholds and the maximum of all monthly upper thresholds for each
variable (keep NDVI's lower bound at 0). Replace the six values in the
script and process every month with the same values.

### Extending the record with new months

Because the normalization bounds are fixed to the current period of record,
appending future months may introduce values outside the existing envelope.
To preserve temporal comparability, recompute the fixed bounds over the
full extended record (not just the new months), update the six values, and
re-export the **entire** archive — not just the new months — with the
updated bounds.

## Output

Each monthly export is a single-band GeoTIFF (`SWATI_YYYYMM.tif`), EPSG:4326,
1 km resolution, values in `[0, 1]`, land pixels only within the study-area
boundary (permanent water masked).

## License

Code released under the [MIT License](LICENSE). MODIS data are subject to
NASA LP DAAC's terms of use and data citation guidelines.

## Citation

If you use this code, please cite:

> Thai, M.-T. & Liou, Y.-A. Surface Water Availability-Temperature Index (SWATI)
> for Global Drought Monitoring. *IEEE Trans. Geosci. Remote Sens.* 63, 1–14,
> https://doi.org/10.1109/tgrs.2025.3541287 (2025).

See also [`CITATION.cff`](CITATION.cff).

## Contact

Minh-Tin Thai ([tmtinrs@gmail.com](mailto:tmtinrs@gmail.com)) and
Yuei-An Liou ([yueian@csrsr.ncu.edu.tw](mailto:yueian@csrsr.ncu.edu.tw))

Center for Space and Remote Sensing Research, National Central University,
Taoyuan, Taiwan

For questions about this code, open an issue on this repository or contact
the authors directly.
