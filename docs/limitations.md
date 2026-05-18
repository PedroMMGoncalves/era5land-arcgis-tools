# Limitations and caveats

What this toolbox does well, and where it should not be trusted
without further checks.

## Spatial

- **Resolution.** ERA5-Land is delivered at 0.1° (about 9 km at the
  equator, 7 km at mid-latitudes). It is a reanalysis on a regular
  grid: not appropriate for sub-kilometre microclimate work, urban
  heat-island studies, or steep-relief sites where the grid cell
  smooths real topography. For those use cases, downscale with a
  high-resolution method (e.g. WRF, CHELSA, or empirical regression)
  before treating outputs as point estimates.
- **Land-only mask.** Oceans are NoData. Coastal pixels are land if
  the cell centre is on land; this can leave island groups (e.g.
  Azores) represented by 1 or 2 pixels each. Interpret with care.
- **Bias near complex terrain and coastlines.** Reanalysis fields
  are constrained by atmospheric observations indirectly; surface
  fields inherit assumptions of the land-surface model. Validate
  against in-situ stations (IPMA in Portugal, AEMET in Spain) before
  applying corrections.
- **No reprojection.** Tool 1 outputs in EPSG:4326. Downstream tools
  assume the same grid throughout a run. Mixing reprojected inputs
  triggers the strict grid-consistency check.

## Temporal

- **Reanalysis, not observations.** ERA5-Land is a model product
  constrained by observations through the parent ERA5 reanalysis.
  Treat trends as the model's representation of climate, not as
  direct measurements.
- **No projection / scenario support.** Future climate (RCP, SSP,
  CMIP6, EURO-CORDEX) is out of scope here. A sibling repo,
  `cordex-arcgis-tools`, is planned for that.
- **No ensemble.** ERA5-Land is a single deterministic dataset. No
  uncertainty quantification is built in. For uncertainty bands,
  pair with ensemble observational products (e.g. EOBS).

## Variable-specific

- **Accumulated variables (`tp`, `ssrd`, `ssr`, `strd`, `str`, `e`,
  `pev`, `ro`).** ERA5-Land hourly stores these as 1-hour
  accumulations whose timestamp marks the END of the accumulation
  period. The companion `era5land-downloads` script handles this
  with a -1h shift before resampling. Daily totals produced this
  way match ECMWF documentation. If a user supplies daily-sum files
  produced by other means, verify the convention.
- **Wind (`u10`, `v10`).** Tool 3 builds magnitude as
  `sqrt(u10^2 + v10^2)` on the fly. This is the scalar wind speed,
  not a vector quantity. Direction is not currently exposed (planned
  as a derived metric in a future version).
- **Snow variables.** Low signal for most of mainland Portugal; use
  only when the region of interest has meaningful snow cover (e.g.
  Serra da Estrela, Pyrenees).

## Software

- **ArcGIS Pro 3.x only.** No support for ArcMap. QGIS users must
  port the algorithms manually (the engine is pure
  numpy / netCDF4 / GDAL and would translate cleanly).
- **No automated test suite.** Validation is empirical: by Tool 2
  inspection, by spot-checks against literature value ranges, and
  by comparing successive versions with raster calculator. Users
  contributing new metrics should follow the same approach.
- **Solo-developer maintenance.** Bug-reports via the GitHub issue
  tracker will get a response, but turnaround depends on availability.
  See `LICENSE` for the formal disclaimer.

## OneDrive (Windows users specifically)

Running Tool 1 / Tool 3 with output folders inside a OneDrive-synced
location can cause:

1. Real-time sync of every TIFF as it is written, hammering CPU and
   network.
2. "Files On-Demand" marking outputs as cloud-only, breaking the
   next read in Tool 3 (the manifest audit in v0.2.1 catches and
   reports this).
3. Conflicts when two processes write to the same folder.

Use a local non-synced path for Tool 1 / Tool 3 outputs (e.g.
`C:\era5land-outputs\`).

## Numerical / methodological

- **Percentile metrics (`p05`, `p10`, `p90`, `p95`)** use linear
  interpolation (NumPy default). For comparability with publications
  that used method-7 from Hyndman & Fan (1996), this matches.
- **NaN propagation.** All aggregations are NaN-aware. A pixel with
  zero valid days returns NaN, not zero. This matters for sparse
  inputs (very short year ranges, heavy masking).
- **Polygon mask rasterisation** uses `CELL_CENTER` assignment.
  Pixels whose centre falls outside the polygon are excluded, even
  if part of the cell area is covered. For very small polygons
  relative to the 0.1° cell, expect zero-pixel coverage and switch
  to a higher-resolution dataset.
