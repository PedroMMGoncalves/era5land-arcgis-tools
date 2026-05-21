# era5land-arcgis-tools

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)
[![ArcGIS Pro](https://img.shields.io/badge/ArcGIS_Pro-3.x-green.svg)](https://www.esri.com/arcgis/products/arcgis-pro/)
[![Status](https://img.shields.io/badge/status-v1.0.0-brightgreen.svg)](CHANGELOG.md)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20275805.svg)](https://doi.org/10.5281/zenodo.20275805)

ArcGIS Pro Python toolbox for processing ERA5-Land daily NetCDF data
into GIS-ready GeoTIFFs and derived climate indicators, with regional
zonal-statistics extraction to long-format CSV.

## Contents

- [Why this exists](#why-this-exists)
- [What is in the toolbox](#what-is-in-the-toolbox)
- [Companion downloader](#companion-downloader)
- [Quickstart](#quickstart)
- [Performance](#performance)
- [Worked example](#worked-example-hdd-climatology-for-portugal)
- [Documentation](#documentation)
- [FAQ](#faq)
- [Dependencies](#dependencies)
- [Citation](#citation)
- [License](#license)
- [Acknowledgments](#acknowledgments)
- [Sibling repositories](#sibling-repositories)

## Why this exists

ERA5-Land is a global 0.1 deg reanalysis dataset (1950 to present)
widely used in climatology, energy, agriculture and hydrology
research. Most climate-indicator toolchains target Python or R users
(xclim, climdex.pcic, ICCLIM, CDO). This toolbox fills the equivalent
niche for ArcGIS Pro users, with two design priorities:

- **Zero pip installs.** Runs on the default `arcgispro-py3` Python
  environment. Drag the `.pyt` into a project and you are done.
- **Variable-agnostic.** Works with any CF-compliant ERA5-Land variable
  out of the box (temperature, precipitation, radiation, wind, soil
  moisture, vegetation indices, snow). New metrics extend a single
  registry without touching the engine.

## What is in the toolbox

| Tool | Role |
| --- | --- |
| **1. Convert NetCDF to GeoTIFF** | Per-day TIFFs from annual NetCDFs. Variable-agnostic, netCDF4 + GDAL fast path (~100x faster than `arcpy.md.MakeNetCDFRasterLayer`). Optional clip polygon. Manifest CSV. |
| **2. Inspect NetCDF** | Debug-only: prints dims, variables, units, decoded time axis. Run before Tool 1 to verify file structure. |
| **3. Compute climate indicators** | 48 metrics across 6 categories (full list in [`docs/metrics.md`](docs/metrics.md)). Annual outputs + full-range climatology + optional IPMA normals + optional WMO standard normals. Auto-unit conversion (K to deg C, m to mm). |
| **4. Extract regional time series** | Zonal mean (and optional min, max, std, count) over polygon regions, for Tool 1 OR Tool 3 outputs (auto-detected). Long-format CSV ready for Pandas / R / SQL. |

Tool 3 metric breakdown:

- **Generic (17):** mean, min, max, sum, p05, p10, p50, p90, p95, p99,
  range, std, etc. Applicable to any variable.
- **Temperature (13):** Spinoni HDD / CDD / GDD, ETCCDI frost_days,
  summer_days, tropical_nights, hot_days, TXn, TNx, etc.; G.POT
  heating / cooling season lengths.
- **Precipitation (5):** annual total, R10mm, R20mm, RX1day, CDD wet.
- **DTR (4):** diurnal temperature range, mean / max / std / extreme.
- **Radiation (4):** annual total, mean daily, peak, photovoltaic
  potential proxy.
- **Wind (5):** mean / max / p95 / hours above threshold / Weibull
  scale.

Climatology windows: full year range (always) + IPMA normals
(1981-2010, 1991-2020) + WMO standard normals (1961-1990, 1981-2010,
1991-2020).

## Companion downloader

[`era5land-downloads`](https://github.com/PedroMMGoncalves/era5land-downloads)
(sibling repo) handles data acquisition from the DestinE Earth Data
Hub (recommended, no queue, 500k req/month free quota) or the
Copernicus CDS API (fallback). See
[`docs/data_acquisition.md`](docs/data_acquisition.md).

## Quickstart

1. **Acquire data** (~30-60 min via EDH). Use
   [`era5land-downloads`](https://github.com/PedroMMGoncalves/era5land-downloads)
   to download annual NetCDFs into a non-OneDrive folder
   (e.g. `C:\era5land\downloads\`).
2. **Add toolbox** (<1 min). In ArcGIS Pro: Catalog pane -> Toolboxes
   -> right-click -> Add Toolbox -> select `era5land_tools.pyt`. Four
   tools appear under "ERA5-Land ArcGIS Tools".
3. **Convert** (~5-15 min per variable per 30 years). Run Tool 1 over
   the downloaded NetCDFs.
4. **Indicators** (~5-30 min per metric per variable). Run Tool 3
   selecting metrics relevant to your study (HDD / CDD /
   heating_season_length / etc.).
5. **Regional extraction** (~1-5 min per polygon FC). Run Tool 4
   with your polygon FC to get a long-format CSV by region.

End-to-end walkthrough: [`docs/workflow.md`](docs/workflow.md).

## Performance

Design choices that drive runtime:

- **Tool 1 fast path.** `netCDF4` + GDAL direct writes instead of
  `arcpy.md.MakeNetCDFRasterLayer`. ~100x speedup empirically on the
  same hardware (no per-slice arcpy overhead, no temp layer
  materialisation).
- **Tool 1 mask cache.** Clip polygon is rasterised once per unique
  grid signature and reused across all subsequent files. For a
  multi-year, multi-variable batch on a fixed AOI, all files after the
  first use a cached mask.
- **Tool 3 parallel TIFF reader.** A `ThreadPoolExecutor` loads N
  daily TIFFs into a pre-allocated `(N, lat, lon)` stack. Default 8
  workers; override with env var `ERA5LAND_READ_WORKERS=N` (set to 1
  to force serial). 4-6x speedup on 8-core machines with SSD-backed
  Tool 1 output. Bottleneck is disk bandwidth.
- **Tool 3 manifest audit.** At startup, drops manifest entries whose
  `.tif` files are missing on disk (OneDrive offload, user deletion).
  Prevents mid-pipeline crashes; surfaces the issue immediately.

Real-world reference: a full 30-year (1991-2020), 4-statistic, single
variable batch over Portugal sub-region (~14 deg x 11 deg, ~140 x 110
pixels at 0.1 deg) takes Tool 1 about 10 minutes and Tool 3 about
10 minutes per metric. The companion downloader took 30 hours of
wall clock to fetch the 76-year (1950-2025) 22-variable source data
that feeds these tools.

## Worked example: HDD climatology for Portugal

Goal: produce per-year and 30-year-climatology rasters of Spinoni
heating degree days (HDD, base 15.5 deg C) for 1991-2020 over
Portugal, plus a long-format CSV with HDD by NUTS-3 region.

1. **Get data** via the
   [downloader](https://github.com/PedroMMGoncalves/era5land-downloads):
   `2m_temperature` daily_mean for 1991-2020 over the Portugal area
   bounding box. Output annual NetCDFs land in
   `C:\era5land\downloads\annual\2m_temperature\daily_mean\`.
2. **Add the toolbox** to your ArcGIS Pro project (Catalog ->
   Toolboxes -> Add Toolbox -> `era5land_tools.pyt`).
3. **Tool 1 -> Convert NetCDF to GeoTIFF.** Input: the 30 annual
   NetCDFs. Output folder: `C:\era5land\outputs\tifs\`. Optional clip:
   your Portugal polygon FC. ~10 min. Produces ~11000 daily TIFFs
   plus a `manifest.csv`.
4. **Tool 3 -> Compute climate indicators.** Inputs folder:
   `C:\era5land\outputs\tifs\` (Tool 1 manifest auto-discovered).
   Metric: `hdd_15_5` (Spinoni HDD, base 15.5 deg C). Year range
   1991-2020. Optional: tick "WMO standard normals" to also get the
   1991-2020 normal. ~10 min. Produces 30 annual HDD rasters +
   1 full-range climatology + 1 WMO normal.
5. **Tool 4 -> Extract regional time series.** Input folder: the
   Tool 3 output. Regions FC: NUTS-3 polygons. Stats: mean (default).
   Output: `hdd_by_nuts3.csv` in long format
   (`region_id, year, metric, value, unit`).
6. **Open the CSV** in Pandas / R / Excel. Sanity check: lowland
   southern NUTS-3 regions should land at ~700-1200 HDD; interior
   highland regions at ~1500-2500 HDD.

## Documentation

- [`docs/metrics.md`](docs/metrics.md) catalogue of every climate
  indicator (formula, units, expected input, peer-reviewed reference)
- [`docs/workflow.md`](docs/workflow.md) end-to-end pipeline
- [`docs/data_acquisition.md`](docs/data_acquisition.md) download
  paths (DestinE EDH primary, CDS API fallback)
- [`docs/limitations.md`](docs/limitations.md) what to trust and what
  to verify before publication
- [`CHANGELOG.md`](CHANGELOG.md) version history

## FAQ

**Why this toolbox instead of xclim / climdex.pcic / ICCLIM?** Same
underlying methodology (ETCCDI, Spinoni, etc.) but no pip / conda
required. Drag the `.pyt` into ArcGIS Pro and it runs. If you live
in Python anyway, xclim is excellent; if your GIS team works in
ArcGIS Pro and does not maintain a Python data-science stack, this
saves a lot of friction.

**Can I add a custom metric?** Yes. Tool 3's `METRIC_REGISTRY` is a
single dict at the top of the relevant section in `era5land_tools.pyt`.
Add an entry with `category`, `label`, `applicable_to`,
`expects_statistic`, `params`, `compute` (a callable on the
`(N, lat, lon)` stack), `units_out`, and `climatology` aggregator.
The UI multi-select populates automatically. No engine code to touch.

**Why Kelvin in Tool 1 output but Celsius in Tool 3 thresholds?** Tool
1 is variable-agnostic and stores RAW source units (K for temperature,
m for precipitation, J/m^2 for radiation totals) so the same converter
works for any future variable. Tool 3 applies unit conversion ONCE per
stack at load time (Kelvin -> Celsius for `KELVIN_VARIABLES`, m -> mm
for `PRECIPITATION_VARIABLES`) and all threshold values in the UI are
in the human-friendly unit (15.5 deg C for HDD, 1 mm for R1mm).

**Does this work with hourly ERA5 / ERA5-Land?** Not directly; Tool 1
expects (time, lat, lon) daily NetCDFs. The companion downloader
produces those by either (EDH) reading hourly Zarr and computing
daily stats locally, or (CDS) requesting the `daily-statistics`
dataset. Hourly data straight from `reanalysis-era5-land` would need
a pre-processing step to aggregate to daily before feeding Tool 1.

**macOS / Linux support?** Untested. ArcGIS Pro is Windows-only, so
the toolbox UI is by definition Windows. The conversion engine
(`netCDF4` + GDAL, no `arcpy` for the hot path) is portable in
principle but the `arcpy.PolygonToRaster` step (for clip polygon
rasterisation) requires `arcpy`.

**ArcGIS Pro 2.x?** Untested. The toolbox was developed against
ArcGIS Pro 3.5 and uses the default `arcgispro-py3` env from that
release (Python 3.11). ArcGIS Pro 2.x ships Python 3.7, which may
or may not have the needed numpy / netCDF4 versions; if you try it
and run into issues, open a GitHub issue with the traceback.

**How do I extend the toolbox to a new ERA5-Land variable?** Tool 1
already accepts any CF-compliant NetCDF: no code change. Tool 3 may
need a new entry in `TEMPERATURE_VARIABLES`, `PRECIPITATION_VARIABLES`,
`RADIATION_VARIABLES`, etc. (variable-type sets near the top of
`era5land_tools.pyt`) so the right unit conversion fires. Then if you
want metrics specific to that variable, add them via `METRIC_REGISTRY`
with the appropriate `applicable_to` set.

## Dependencies

- ArcGIS Pro 3.x
- `arcgispro-py3` default Python (`arcpy`, `numpy`, `netCDF4`,
  `osgeo.gdal`)
- No pip installs, no conda env clone

Optional environment variable `ERA5LAND_READ_WORKERS=N` controls the
Tool 3 parallel TIFF reader (default 8 threads). GDAL itself reads
`GDAL_NUM_THREADS=ALL_CPUS` (set at module load) to multi-thread the
DEFLATE / LZW decompression inside each TIFF read; the two layers
compose well.

## Citation

If you use this toolbox in published work, please cite both the
software and the underlying dataset (see [`CITATION.cff`](CITATION.cff)
for machine-readable form).

Software:

```text
Goncalves P (2026). era5land-arcgis-tools: ArcGIS Pro toolbox for
processing ERA5-Land daily NetCDF data into climate indicators and
regional time series. Version 1.0.0. Zenodo. DOI: 10.5281/zenodo.20275805
```

ERA5-Land dataset:

```text
Munoz-Sabater J et al. (2021). ERA5-Land: a state-of-the-art global
reanalysis dataset for land applications. Earth System Science Data
13:4349-4383. DOI: 10.5194/essd-13-4349-2021
```

Specific metrics carry their own peer-reviewed citations (Spinoni for
HDD/CDD, Casasso & Sethi 2016 for G.POT season-length inputs, Frich
et al. 2002 / Zhang et al. 2011 for ETCCDI). See
[`docs/metrics.md`](docs/metrics.md).

## License

[Apache License 2.0](LICENSE).

## Acknowledgments

- Copernicus Climate Change Service (C3S), ECMWF, and the
  [Destination Earth](https://destination-earth.eu/) initiative for
  hosting ERA5-Land and the EDH Zarr mirror.
- Sibling
  [`ecde-arcgis-tools`](https://github.com/PedroMMGoncalves/ecde-arcgis-tools)
  repo for the shared engine pattern (netCDF4 + GDAL fast path, mask
  cache, OID workaround).

## Sibling repositories

| Repo | Dataset | Status |
| --- | --- | --- |
| [`ecde-arcgis-tools`](https://github.com/PedroMMGoncalves/ecde-arcgis-tools) | ECDE pre-aggregated indicators (Spinoni HDD/CDD, 0.25 deg) | Released v1.0.1 |
| **[`era5land-arcgis-tools`](https://github.com/PedroMMGoncalves/era5land-arcgis-tools)** | ERA5-Land daily, raw (0.1 deg) | **This repo, v1.0.0** |
| [`era5land-downloads`](https://github.com/PedroMMGoncalves/era5land-downloads) | Downloader for the above (EDH or CDS) | Released v1.0.0 |
| `cordex-arcgis-tools` | EURO-CORDEX projections with RCP scenarios | Planned |
