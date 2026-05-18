# Data acquisition

The `era5land-arcgis-tools` toolbox consumes annual NetCDF files of
ERA5-Land daily variables. Two acquisition paths are supported.

## Path A: DestinE Earth Data Hub (recommended)

The DestinE Earth Data Hub (EDH) exposes the raw ERA5-Land hourly Zarr
store over plain HTTPS, with no queue. Daily statistics are computed
locally from the hourly source. Typical throughput: a 30-year
climatology of one variable for Portugal completes in 30 to 60 minutes
of wall clock.

Setup:

1. Register a free DestinE Platform account at
   <https://platform.destine.eu/> (academic affiliation in EU member
   states qualifies under the "Academia and Research" category).
2. Generate a Personal Access Token at
   <https://earthdatahub.destine.eu/account-settings>.
3. Run the companion `era5land-downloads` repo's `00_setup.bat` which
   creates a dedicated conda env, prompts for the token and writes
   `~/_netrc`, and exposes `download_era5land_edh.py` for batch
   downloads.
4. Outputs land in `downloads/annual/<variable>/<statistic>/`
   matching the layout Tool 1 expects.

Quotas as of 2026-05: 500 000 requests per month per user. A
30-year Portugal batch with the full 22-variable bundle uses on the
order of 30 000 chunk requests, well within budget.

## Path B: Copernicus Climate Data Store (CDS) API

The traditional path. The `derived-era5-land-daily-statistics`
dataset is computed on demand by the CDS service; there is no
server-side cache.

Setup:

1. Register a free CDS account at <https://cds.climate.copernicus.eu>.
2. Generate a Personal Access Token in your CDS profile.
3. Accept the dataset's Terms of Use.
4. Use the companion `era5land-downloads` repo's `download_era5land.py`
   (uses `cdsapi`, supports state persistence and retry).

Limits: the CDS forces one month per request for this dataset. A
30-year batch of one variable plus three statistics is 1080 requests.
Median throughput in 2026 is one request every 30 to 80 minutes during
European business hours; expect 2 to 5 days of wall clock. Use this
path only if EDH does not host a variable you need, or if you cannot
register for DestinE.

## Pre-existing NetCDFs from other sources

Tool 1 is variable-agnostic and reads any CF-compliant NetCDF with
`(time, lat, lon)` (or `(valid_time, latitude, longitude)`)
dimensions on a regular WGS84 grid. ECDE-derived files, custom
aggregations, or CHELSA / TerraClimate-style inputs work without
modification provided the grid is consistent across files and the
units are CF-compliant.

## Verifying the downloaded files

Before launching a large Tool 1 run, sample one input through Tool 2
(Inspect). Confirm:

- Time dimension parses (`valid_time` or `time`)
- Variable units (`K` for temperature, `m` for precipitation, etc.)
- Grid extent matches the area you requested
- No silent truncation (e.g. 8760 hourly slots present for a year)

If Tool 2 reports `unknown` for `standard_name`, that is expected for
ERA5-Land daily statistics; Tool 3 infers variable type from the
short name (`t2m`, `tp`, etc.).
