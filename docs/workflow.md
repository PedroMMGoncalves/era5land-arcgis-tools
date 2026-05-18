# End-to-end workflow

This document walks through the full pipeline from raw downloads to
publication-ready outputs. Four tools in sequence.

## Pipeline overview

```
[Source]                  [Tool]                    [Output]
─────────                 ──────                    ─────────
EDH or CDS                                          downloads/annual/
  Zarr / API   ─────►  era5land-downloads  ─────►    <var>/<stat>/
  (companion repo)                                   era5land_<var>_<stat>_<year>.nc

                                                    output_t1/
                       Tool 1: Convert    ─────►    <var>/<stat>/<year>/
                       NetCDF -> GeoTIFF            era5land_<var>_<stat>_YYYY-MM-DD.tif

                                                    output_t3/
                       Tool 3: Compute    ─────►    annual/<var>/<stat>/<metric>/
                       climate indicators            <metric>_<year>.tif
                                                    climatology/<var>/<stat>/
                                                     <metric>_climatology_<period>.tif

                                                    output_t4.csv
                       Tool 4: Extract    ─────►    region_id, region_label,
                       regional time series          variable, statistic, metric,
                                                     period_kind, year, ...,
                                                     zonal_mean, [zonal_min, ...]
```

Tool 2 (Inspect) is a debug-only side path; not part of the production
flow.

## Step 1: Download annual NetCDFs

See [`data_acquisition.md`](data_acquisition.md). End state:

```
downloads/annual/2m_temperature/daily_mean/
   era5land_2m_temperature_mean_1991.nc
   era5land_2m_temperature_mean_1992.nc
   ...
   era5land_2m_temperature_mean_2020.nc
```

## Step 2: Tool 1 (Convert NetCDF -> GeoTIFF)

Inputs:

- NetCDF files (multi-select; drag and drop the `downloads/annual/`
  tree or pick individual files).
- Output folder (NOT inside OneDrive; see `limitations.md`).
- Clip polygon (optional, Portugal continental + islands recommended
  for IPMA-comparable outputs).

Outputs:

- Per-day TIFFs in `<output>/<variable>/<statistic>/<year>/`.
- `manifest.csv` listing every TIFF with its variable, statistic,
  year, date, units, standard_name, long_name, source_nc.

Time: about 10 to 20 minutes for a 30-year single-variable batch on
an SSD; longer with multiple variables.

Validation: open one TIFF per year in ArcGIS Pro. For `2m_temperature
daily_mean`, expect values around 275 to 305 K (subtract 273.15 for
degC). NoData over oceans (ERA5-Land is land-only).

## Step 3: Tool 3 (Compute climate indicators)

Inputs:

- Tool 1 output folder (the one with `manifest.csv`).
- Variable + statistic selected from the auto-populated dropdowns.
- Metrics multi-select from the registry (see
  [`metrics.md`](metrics.md) for the full catalogue).
- Optional clip polygon (must match the Tool 1 polygon if used).
- Climatology and IPMA-normals toggles.

Outputs:

- `annual/<variable>/<statistic>/<metric>/<metric>_<year>.tif`
- `climatology/<variable>/<statistic>/<metric>_climatology_<period>.tif`
- `ipma_normals/...` if enabled.
- `manifest.csv` with `metric`, `period_kind`, `start_year`,
  `end_year` columns added relative to Tool 1.

Time: about 5 to 15 minutes per metric per (variable, statistic)
combination on a multi-core workstation with the default parallel
TIFF reader (8 workers).

Unit conversion happens automatically:

- Kelvin temperature variables converted to degC before metric runs.
- Precipitation (`tp`) converted from m to mm.
- Radiation kept raw (J/m^2/day); per-metric factors handle the
  conversion to W/m^2 or kWh/m^2 as appropriate.
- Wind multi-input metrics build `sqrt(u10^2 + v10^2)` in m/s on
  the fly when both `u10` and `v10` are available at the selected
  statistic.

## Step 4: Tool 4 (Extract regional time series)

Inputs:

- Tool 1 OR Tool 3 output folder (auto-detected via the manifest's
  schema; presence of a `metric` column => Tool 3).
- Polygon FC with one polygon per region (NUTS, municipality,
  basin, custom zones).
- Integer field for region ID (SHORT, LONG, or OBJECTID).
- Optional text field for human-readable region label.
- Output CSV path.
- Optional variable / statistic / metric / year filters.
- Multi-select zonal stats: mean (default), min, max, std, count.

Output: long-format CSV with one row per (region, variable, statistic,
metric, period_kind, time_key). Ready for `read_csv` in Pandas or R.

Time: a few minutes for typical inputs (hundreds of TIFFs across
dozens of regions).

## Reproducibility

The end-to-end run is fully captured by:

1. The downloader configuration (`download_era5land_edh.py` source
   with `DOWNLOADS`, `YEAR_RANGE`, `AREA` block visible at top).
2. The Tool 1 output `manifest.csv` (every TIFF with source NetCDF).
3. The Tool 3 output `manifest.csv` (every metric with parameters
   implicit in the metric registry entry; values logged at run time).
4. The Tool 4 output CSV (every zonal datum with region identifier).

For citation in a methods section, reference this toolbox's Zenodo
DOI and the ERA5-Land dataset DOI (see README).
