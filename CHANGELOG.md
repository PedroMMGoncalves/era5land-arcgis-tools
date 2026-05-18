# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned

- Optional: trend analysis (Theil-Sen + Mann-Kendall) per pixel.
- Optional: structured JSON logging for reproducibility.
- Optional: bioclimatic variables BIO1 to BIO19 (multi-variable, deferred to v2.0+).

## [1.0.0] - 2026-05-18

First public release. All Entrega A + B + C items complete, validated
against the Portuguese G.POT use case.

### Added

- **Tool 4: Extract regional time series** (`ExtractRegionalTimeSeries`).
  Computes zonal statistics over polygon regions for each TIFF produced
  by Tool 1 (daily values) or Tool 3 (annual indicators, climatology,
  IPMA / WMO normals). Emits a long-format CSV. Auto-detects input source
  from `manifest.csv` (presence of `metric` column => Tool 3).
  Helpers: `rasterise_regions_to_labels`, `get_or_compute_region_labels`,
  `compute_zonal_stats` (numpy bincount + ufunc.at), `discover_t4_inputs`.
- **Tool 3: G.POT-oriented temperature metrics** (`heating_season_length`,
  `cooling_season_length`). Count of days with `T_mean` below/above a
  user-tunable threshold (defaults 15.5 / 22.0 degC, aligned with
  Spinoni HDD/CDD bases). Designed as direct input for shallow geothermal
  potential mapping (Casasso & Sethi 2016, doi:10.1016/j.energy.2016.03.091).
- **Tool 3: WMO standard normals checkbox** (`compute_wmo_normals`).
  Three climatology checkboxes now exposed:
  - **"Compute climatology"** (default ON): mean over the full input
    year range -> `climatology/<variable>/<statistic>/`.
  - **"Compute IPMA normals"** (default OFF): 1981-2010 + 1991-2020 ->
    `ipma_normals/<variable>/<statistic>/`.
  - **"Compute WMO standard normals"** (default OFF): 1961-1990 +
    1981-2010 + 1991-2020 -> `wmo_normals/<variable>/<statistic>/`.
  Overlapping periods are computed twice and written to both subfolders
  when both checkboxes are enabled (intentional per-checkbox separation).
- **Tool 3: dynamic parameter enable/disable**. The advanced threshold
  and base parameters (generic, HDD/CDD/GDD base, precipitation,
  PV, radiation, wind, heating/cooling season) auto-enable only when a
  metric in their category is selected. Reduces UI noise for
  single-domain runs (e.g. only temperature). Defensive `except IndexError`
  fails open if the parameter list shape changes in the future.
- **Internal refactor**: extracted shared `_emit_period_set` helper inside
  `_emit_climatology` so IPMA-normals and WMO-normals branches share one
  implementation.
- **`CITATION.cff`** for machine-readable software citation and Zenodo
  DOI generation on release.
- **`docs/data_acquisition.md`**: dual-path download instructions
  (DestinE EDH primary, CDS API fallback) with quota notes.
- **`docs/workflow.md`**: end-to-end pipeline walkthrough from raw
  NetCDFs through Tools 1, 3 and 4.
- **`docs/limitations.md`**: spatial, temporal, variable-specific,
  software, OneDrive, and numerical caveats with mitigations.
- **README polished**: 4-tool table, sibling repository map, citation
  block, license and ArcGIS badges.

### Notes

- New Tool 3 parameters appended at the END of the return list to
  preserve all pre-existing positional indices in `execute()`. Heating
  threshold = index 21, cooling threshold = 22, compute_wmo_normals = 23.
- Tool 4 long-CSV schema: `region_id, region_label, variable, statistic,
  metric, period_kind, year, start_year, end_year, date, units` plus
  one `zonal_<stat>` column per requested statistic.
- For Tool 4, region ID field must be integer-typed (SHORT, LONG, OID);
  OID values are copied to a LONG proxy to dodge ERROR 003658 (same
  workaround as `rasterise_clip_mask`, Tool 1 design decision #7).
- See companion `era5land-downloads` repo for the data acquisition
  step; this toolbox starts at Tool 1 with annual NetCDFs already on
  disk.

## [0.3.0] - 2026-05-18

### Added

- **Tool 4: Extract regional time series** (`ExtractRegionalTimeSeries`).
  Computes zonal statistics over polygon regions for each TIFF produced
  by Tool 1 (daily values) or Tool 3 (annual indicators, climatology,
  IPMA normals). Emits a long-format CSV ready for downstream analysis
  in R / Pandas / SQL. Auto-detects input source from manifest.csv
  (presence of `metric` column => Tool 3).
- New helpers:
  - `rasterise_regions_to_labels`: polygon FC -> int label grid
    (one region id per pixel), with 64-bit OID/GRID workaround.
  - `get_or_compute_region_labels`: cached wrapper keyed by
    (grid signature, FC path, region_id_field).
  - `compute_zonal_stats`: vectorised mean/min/max/std/count per region
    via numpy bincount + ufunc.at (no Python loops over pixels).
  - `discover_t4_inputs`: unified manifest discovery returning a uniform
    item list across Tool 1 and Tool 3 sources.
  - `_safe_int`: small helper for tolerant int parsing of manifest cells.

### Notes

- Tool 4 long-CSV schema (always present columns):
  `region_id, region_label, variable, statistic, metric, period_kind,
  year, start_year, end_year, date, units`
  plus one `zonal_<stat>` column per requested statistic.
- For Tool 1 inputs `metric` is empty and `period_kind = 'daily'`;
  for Tool 3 inputs `period_kind` is one of `annual`, `climatology`,
  `ipma_normals_1981_2010`, `ipma_normals_1991_2020`.
- Region ID field must be integer-typed (SHORT/LONG/OID); OID values
  are copied to a LONG proxy to dodge ERROR 003658 (64-bit OID -> GRID).
  Same workaround as `rasterise_clip_mask` (Tool 1 design decision #7).
- Optional `region_label` column populated from a user-picked text
  field on the polygon FC (kept blank if not provided).

## [0.2.2] - 2026-05-18

### Added

- **Tool 3: two new temperature metrics** for the G.POT shallow
  geothermal potential workflow (Casasso & Sethi 2016,
  doi:10.1016/j.energy.2016.03.091):
  - `heating_season_length`: count of days with `T_mean <
    heating_threshold` (degC, default 15.5, parameter
    `heating_threshold`). Designed as direct input for the heating-
    season `tc` parameter in the G.POT empirical correlation.
  - `cooling_season_length`: count of days with `T_mean >
    cooling_threshold` (degC, default 22.0, parameter
    `cooling_threshold`). Designed for G.POT cooling-mode mapping.
- Two new UI parameters on `ComputeERA5LandIndicators`
  (`heating_threshold`, `cooling_threshold`), appended at the END of
  the parameter list to avoid shifting positional indices of
  pre-existing parameters in `execute()`.
- `docs/metrics.md` entries with formulas, parameter semantics,
  use-case rationale, and peer-reviewed reference.

### Notes

- Both metrics are temperature-specific and expect `daily_mean` input.
  They differ from existing fixed-threshold ETCCDI metrics (`frost_days`
  / `ice_days` / `summer_days` / `tropical_nights` / `hot_days`) by
  using `T_mean` (not Tmin/Tmax) and by exposing the threshold as a
  user-configurable parameter.
- Defaults align with Spinoni HDD/CDD base temperatures so the same
  threshold semantics carry across `hdd_spinoni`, `cdd_spinoni`, and
  the new season-length metrics.

## [0.2.1] - 2026-05-16

### Performance
- **Parallel TIFF reads** in Tool 3 via `ThreadPoolExecutor`. The new
  `read_paths_into_stack(paths)` helper streams reads into a
  pre-allocated stack buffer using configurable worker threads
  (default 8, override via env var `ERA5LAND_READ_WORKERS=N`; set to
  1 to force serial). GDAL is thread-safe in read-only mode and Python
  releases the GIL during the C-level I/O, so this scales near-linearly
  until disk bandwidth saturates. Empirical 4-6x speedup on multi-core
  workstations; ~80 s saved per year-stack of 365 daily TIFFs.
- **Stack pre-allocation** in `read_paths_into_stack`: streams each
  read directly into its slot in the destination buffer instead of
  building a Python list of N arrays and then `np.stack`-ing.
  Reduces peak memory during stack build from ~2x to ~1x of the
  final stack size.

### Robustness
- **Manifest filesystem audit** (`_audit_index_against_disk`):
  `build_inputs_index` now verifies every TIFF path listed in the Tool 1
  manifest exists on disk, drops missing entries with a warning, and
  errors out if all entries vanished. Prevents crashes mid-run when
  the user has deleted some Tool 1 outputs between runs or OneDrive
  has offloaded files marked as "online only".

### UX
- **Empirical ETA in Tool 3 progressor**: after the first year of any
  pipeline (Regular, DTR, Wind) completes, the per-year progressor
  label includes a rolling ETA (e.g. `Regular year 5/30: 2024 (ETA
  12m30s)`). Computed via `ComputeERA5LandIndicators._eta_label`.
- **Tool 1 default log consolidation**: without verbose, each NetCDF
  emits a single compact one-liner (`<var> <statistic> (<units>),
  <n_steps> steps <year_first>-<year_last>`) instead of 4-5 separate
  lines. Verbose mode keeps the full multi-line breakdown. Reduces
  log noise by ~75% for typical 30-year batches without losing debug
  visibility.

### Configuration
- New environment variable `ERA5LAND_READ_WORKERS` (default 8) controls
  thread count for parallel TIFF reads.

## [0.2.0] - 2026-05-16
- Tool 4: Extract regional time series (zonal mean over polygons)
- docs/data_acquisition.md (step-by-step CDS download guide)
- docs/workflow.md (end-to-end pipeline documentation)
- docs/limitations.md (caveats: resolution, reanalysis bias, etc.)
- CITATION.cff for Zenodo DOI on release
- Optional: trend analysis (Theil-Sen + Mann-Kendall) per pixel
- Optional: structured JSON logging for reproducibility

### Planned for v1.0.0
- Public release after validation against real Portuguese data
- Zenodo DOI
- Complete documentation

### Deferred to v2.0+
- BIO1-19 bioclimatic variables (multi-variable derived metrics)

## [0.2.0] - 2026-05-16

### Added
- **Tool 3: Compute climate indicators** (`ComputeERA5LandIndicators`)
  with 46 metrics across 6 categories:
  - **Generic (17):** mean, min, max, range, std, sum, median, p10,
    p90, count_above_threshold, count_below_threshold,
    seconds_above_threshold, seconds_below_threshold,
    accumulation_above_threshold, accumulation_below_threshold,
    max_consecutive_above, max_consecutive_below.
  - **Temperature-specific (11):** hdd_spinoni, cdd_spinoni, gdd,
    frost_days, ice_days, summer_days, tropical_nights, hot_days,
    gsl (WMO NH variant, cumsum windowed), txn, tnx.
  - **Precipitation-specific (5):** total_precip, wet_days,
    consecutive_dry_days, consecutive_wet_days, max_1day_precip.
  - **Radiation-specific (4):** radiation_w_per_m2_mean,
    radiation_kwh_per_m2_year, photovoltaic_potential,
    radiation_above_threshold_days.
  - **Wind-specific (5, multi-input u10+v10):** wind_speed_mean,
    wind_speed_max, wind_speed_p90, wind_power_density, calm_days.
  - **DTR (4, multi-input daily_max+daily_min):** dtr_mean, dtr_max,
    dtr_p90, dtr_std.
- Single source of truth for metrics via `METRIC_REGISTRY` spec table.
- Auto-conversion Kelvin to Celsius for temperature variables; m to mm
  for precipitation. Radiation kept in raw J/m^2/day with per-formula
  factors. Wind builds sqrt(u^2+v^2) magnitude on the fly.
- Annual outputs + optional climatology (mean of annual values) over
  the full year range.
- Optional IPMA normals emulation: 1981-2010 and 1991-2020 windows,
  with warnings on partial coverage.
- Auto-discovery of (variable, statistic, year) groups via Tool 1
  manifest.csv; filesystem walk fallback when manifest absent.
- DTR auto-discovery: when a temperature variable has both
  daily_maximum and daily_minimum present, DTR metrics process them
  pairwise without extra UI configuration.
- Wind auto-discovery: when u10 and v10 are both present at the
  selected statistic, wind metrics process them pairwise.
- Validation warnings for variable-type incompatibility and
  statistic-expectation mismatches per metric.
- Optional clip polygon support with the same mask cache as Tool 1.
- docs/metrics.md catalogue with definitions, formulas, units and
  peer-reviewed references for every metric.

### Changed
- **Tool 1 output layout** changed to include the statistic level in
  the path: `output/<variable>/<statistic>/<year>/<variable>_<statistic>_YYYY-MM-DD.tif`.
  Statistic is auto-detected from the source NetCDF filename using
  the pattern `era5land_<variable>_<short_stat>_<year>.nc` produced
  by the `era5land-downloads` companion concat step. Falls back to
  `unknown_stat` if not detectable. Prevents silent overwrites when
  processing daily_mean + daily_maximum + daily_minimum of the same
  variable in one run.
- **Tool 1 manifest.csv** gained a `statistic` column between
  `variable` and `standard_name`.

### Technical
- Same netCDF4 + GDAL fast-path philosophy reused: Tool 3 reads
  daily TIFFs via `gdal.Open` + `ReadAsArray`, writes via
  `write_geotiff_gdal`. No `arcpy.RasterToNumPyArray` overhead.
- NaN-aware implementations throughout; helper `_nan_if_all_invalid`
  distinguishes "zero by computation" from "no valid data".
- GSL implemented with NumPy cumsum-based sliding-window detection
  rather than per-pixel loop.
- Tool 3 output structure aligned with the planned Tool 4 input:
  per-(variable, statistic, metric) sub-trees in `annual/` simplify
  zonal-mean time series extraction.
- Cosmetic NumPy `RuntimeWarning: All-NaN slice encountered` (emitted
  by nanmean/nanmin/etc. on fully-masked pixels like ocean) suppressed
  at module load via `warnings.filterwarnings`.
- Cross-pipeline grid validation: Tool 3 verifies that DTR and wind
  pipeline TIFFs share the geotransform of the primary sample. If not,
  the affected pipeline is skipped with a warning rather than producing
  outputs with wrong georeferencing.
- Regular pipeline failure (missing (variable, statistic) key or empty
  year range) degrades to a warning and skips the regular metrics
  rather than aborting the run; DTR and wind pipelines still execute
  if their inputs resolved.
- Radiation metric specs declare `expects_statistic="daily_sum"` and
  docs/metrics.md documents the CDS daily_mean accumulation gotcha
  explicitly. Users should source radiation via
  `derived-era5-land-daily-aggregations` (daily_sum), not
  `derived-era5-land-daily-statistics` (daily_mean), to get correct
  energy units.

## [0.1.0] - 2026-05-15

### Added
- Initial skeleton release (Entrega A).
- Tool 1: Convert ERA5-Land NetCDF to GeoTIFF (variable-agnostic).
- Tool 2: Inspect ERA5-Land NetCDF (debug).
- netCDF4 + GDAL fast path conversion engine.
- Auto-detection of CF-compliant time / lat / lon / data variables.
- Multi-NetCDF input handling with grid consistency validation.
- Per-variable, per-year subfolder organisation of outputs.
- Optional clip polygon with bbox slicing and polygon mask cache.
- DEFLATE-compressed GeoTIFF output.
- Manifest CSV with full metadata (variable, units, dates, source NetCDF).
- Verbose logging option.
- README, LICENSE (Apache-2.0), .gitignore.

### Technical
- Toolbox uses default `arcgispro-py3` environment, no pip installs.
- Dependencies: arcpy, numpy, netCDF4, osgeo.gdal (all bundled with ArcGIS Pro 3.x).
- Mask rasterisation uses temporary SHORT field workaround for 64-bit OIDs
  (same pattern as ecde-arcgis-tools v1.0.1).
