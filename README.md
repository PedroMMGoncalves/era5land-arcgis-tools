# era5land-arcgis-tools

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)
[![ArcGIS Pro](https://img.shields.io/badge/ArcGIS_Pro-3.x-green.svg)](https://www.esri.com/arcgis/products/arcgis-pro/)
[![Status](https://img.shields.io/badge/status-v1.0.0-brightgreen.svg)](CHANGELOG.md)

ArcGIS Pro Python toolbox for processing ERA5-Land daily NetCDF data
into GIS-ready GeoTIFFs and derived climate indicators, with regional
zonal-statistics extraction to long-format CSV.

## Why this exists

ERA5-Land is a global 0.1° reanalysis dataset (1950 to present) widely
used in climatology, energy, agriculture and hydrology research. Most
climate-indicator toolchains target Python or R users (xclim,
climdex.pcic, ICCLIM, CDO). This toolbox fills the equivalent niche
for ArcGIS Pro users, with two design priorities:

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
| **3. Compute climate indicators** | 48 metrics in 6 categories: Generic (17), Temperature (13, includes Spinoni HDD/CDD/GDD, ETCCDI extremes, G.POT heating / cooling season lengths), Precipitation (5), DTR (4), Radiation (4), Wind (5). Annual outputs + full-range climatology + optional IPMA normals (1981-2010, 1991-2020) + optional WMO standard normals (1961-1990, 1981-2010, 1991-2020). Auto-unit conversion (K to degC, m to mm). |
| **4. Extract regional time series** | Zonal mean (and optional min, max, std, count) over polygon regions, for Tool 1 OR Tool 3 outputs (auto-detected). Long-format CSV ready for Pandas / R / SQL. |

The full metric catalogue with formulas, expected input statistic, and
peer-reviewed references is in [`docs/metrics.md`](docs/metrics.md).

## Companion downloader

`era5land-downloads` (sibling repo) handles data acquisition from the
DestinE Earth Data Hub (recommended, no queue) or the Copernicus CDS
API (fallback). See [`docs/data_acquisition.md`](docs/data_acquisition.md).

## Quickstart

1. **Acquire data.** Use `era5land-downloads` to download annual
   NetCDFs into a non-OneDrive folder (e.g. `C:\era5land\downloads\`).
2. **Add toolbox.** In ArcGIS Pro: Catalog pane -> Toolboxes ->
   right-click -> Add Toolbox -> select `era5land_tools.pyt`. Four
   tools appear under "ERA5-Land ArcGIS Tools".
3. **Convert.** Run Tool 1 over the downloaded NetCDFs.
4. **Indicators.** Run Tool 3 selecting metrics relevant to your study
   (HDD / CDD / heating_season_length / etc.).
5. **Regional extraction.** Run Tool 4 with your polygon FC to get a
   long-format CSV by region.

End-to-end walkthrough: [`docs/workflow.md`](docs/workflow.md).

## Documentation

- [`docs/metrics.md`](docs/metrics.md) catalogue of every climate
  indicator (formula, units, expected input, peer-reviewed reference)
- [`docs/workflow.md`](docs/workflow.md) end-to-end pipeline
- [`docs/data_acquisition.md`](docs/data_acquisition.md) download
  paths (DestinE EDH primary, CDS API fallback)
- [`docs/limitations.md`](docs/limitations.md) what to trust and what
  to verify before publication
- [`CHANGELOG.md`](CHANGELOG.md) version history

## Dependencies

- ArcGIS Pro 3.x
- `arcgispro-py3` default Python (`arcpy`, `numpy`, `netCDF4`,
  `osgeo.gdal`)
- No pip installs, no conda env clone

Optional environment variable `ERA5LAND_READ_WORKERS=N` controls the
Tool 3 parallel TIFF reader (default 8 threads).

## Citation

If you use this toolbox in published work, please cite both the
software and the underlying dataset (see [`CITATION.cff`](CITATION.cff)
for machine-readable form).

Software:

```text
Goncalves P (2026). era5land-arcgis-tools: ArcGIS Pro toolbox for
processing ERA5-Land daily NetCDF data into climate indicators and
regional time series. Version 1.0.0. Zenodo. DOI: <pending>
```

ERA5-Land dataset:

```text
Munoz-Sabater J et al. (2021). ERA5-Land: a state-of-the-art global
reanalysis dataset for land applications. Earth System Science Data
13:4349-4383. DOI: 10.5194/essd-13-4349-2021
```

Specific metrics carry their own peer-reviewed citations (Spinoni for
HDD/CDD, Casasso & Sethi 2016 for G.POT season-length inputs, Frich
et al. 2002 / Zhang et al. 2011 for ETCCDI). See `docs/metrics.md`.

## License

[Apache License 2.0](LICENSE).

## Acknowledgments

- Copernicus Climate Change Service (C3S), ECMWF, and the
  Destination Earth initiative for hosting ERA5-Land and the EDH
  Zarr mirror.
- Companion `ecde-arcgis-tools` repo for the shared engine pattern
  (netCDF4 + GDAL fast path, mask cache, OID workaround).

## Sibling repositories

| Repo | Dataset | Status |
| --- | --- | --- |
| `ecde-arcgis-tools` | ECDE pre-aggregated indicators (Spinoni HDD/CDD, 0.25°) | Released v1.0.1 |
| **`era5land-arcgis-tools`** | ERA5-Land daily, raw (0.1°) | **This repo, v1.0.0** |
| `era5land-downloads` | Downloader for the above (EDH or CDS) | Companion |
| `cordex-arcgis-tools` | EURO-CORDEX projections with RCP scenarios | Planned |
