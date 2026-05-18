# -*- coding: utf-8 -*-
"""
ERA5-Land ArcGIS Tools
======================

ArcGIS Pro Python toolbox for processing ERA5-Land daily NetCDF data from
the Copernicus Climate Data Store (CDS) into GIS-ready GeoTIFFs and
derived climate indicators.

Tools provided
--------------
1. Convert ERA5-Land NetCDF to GeoTIFFs
2. Inspect ERA5-Land NetCDF (debug)
3. Compute climate indicators (Entrega B, v0.2.0)
4. Extract regional time series (Entrega C, to be added in v0.3.0)

Changelog
---------
v0.2.1 (2026-05-16)
  - Performance: parallel TIFF reads in Tool 3 via ThreadPoolExecutor
    (env var ERA5LAND_READ_WORKERS=N, default 8). 4-6x I/O speedup.
  - Performance: read_paths_into_stack pre-allocates the destination
    stack and streams reads into slots (peak memory ~1x stack size).
  - Robustness: build_inputs_index now audits the manifest against
    disk and drops entries whose .tif files are missing (avoids
    mid-run crashes when files were deleted or OneDrive-offloaded).
  - UX: Tool 3 progressor labels include a rolling empirical ETA
    after the first year of each pipeline completes.
  - UX: Tool 1 default logs consolidated to a single per-NetCDF line
    (verbose mode keeps the multi-line breakdown for debugging).

v0.2.0 (2026-05-16)
  - Added Tool 3: Compute climate indicators with 46 metrics across
    6 categories (Generic 17, Temperature 11, Precipitation 5, DTR 4,
    Radiation 4, Wind 5). Single source of truth via METRIC_REGISTRY
    spec table.
  - Auto-unit conversion in Tool 3: K to degC for temperature, m to mm
    for precipitation; radiation kept raw with per-formula factors;
    wind builds sqrt(u^2+v^2) magnitude on the fly.
  - Annual outputs + optional climatology + optional IPMA normals
    (1981-2010, 1991-2020).
  - Multi-input metrics (DTR, wind) auto-discover from manifest.
  - Tool 1 path layout changed to include statistic level:
    output/<variable>/<statistic>/<year>/<variable>_<statistic>_YYYY-MM-DD.tif
    Statistic auto-detected from source NetCDF filename pattern
    (era5land_<variable>_<stat>_<year>.nc); falls back to 'unknown_stat'.
    Prevents silent overwrites when processing multiple statistics of
    the same variable in one run.
  - Tool 1 manifest.csv gained `statistic` column.
  - docs/metrics.md catalogue with formulas and peer-reviewed
    references for every metric.

v0.1.0 (2026-05-15)
  - Initial skeleton: Tool 1 (Convert) and Tool 2 (Inspect).
  - Variable-agnostic conversion engine using netCDF4 + GDAL fast path.
  - Auto-detection of CF-compliant time, lat, lon, variable.
  - Multi-NetCDF input handling with grid consistency validation.
  - Per-variable, per-year subfolder organisation:
    output/<variable>/<year>/<file>.tif (replaced in v0.2.0).
  - Optional clip polygon with bbox slicing + polygon mask cache.
  - DEFLATE-compressed GeoTIFF output via GDAL.
  - Verbose logging option.
  - Manifest CSV with full metadata for downstream tools.

Expected input
--------------
ERA5-Land daily statistics NetCDF files downloaded from:
  https://cds.climate.copernicus.eu/datasets/derived-era5-land-daily-statistics

Typical filename pattern (when downloaded as single-year requests):
  <variable>_<year>.nc, or
  era5land_<variable>_<year>.nc, or
  similar. The toolbox does NOT depend on filename for time parsing;
  it reads the time dimension via netCDF4.num2date for robust handling
  of all CF calendars.

Expected NetCDF structure:
  Dimensions: (time, latitude, longitude) or (time, lat, lon) or (valid_time, ...)
  Variables: one or more 2D-time-varying fields with CF standard_name
  CRS: regular WGS84 lat/lon grid, 0.1 deg cell size (EPSG:4326)
  Time units: typically "hours since 1900-01-01 00:00:00.0" with calendar
              "gregorian" or "proleptic_gregorian"

Design decisions worth preserving
---------------------------------
1. Variable-agnostic conversion: Tool 1 reads any CF-compliant variable
   with (time, lat, lon) dimensions. Auto-detection prefers data
   variables over coordinate/bound variables. User can override.

2. Subfolder organisation (output/<variable>/<statistic>/<year>/): scales
   to multi-variable, multi-statistic, multi-year batches without filename
   collisions. Statistic level prevents silent overwrites when the same
   variable is processed for daily_mean + daily_maximum + daily_minimum
   in the same run.

3. netCDF4 + GDAL fast path (same engine pattern as ecde-arcgis-tools).
   No arcpy.md.MakeNetCDFRasterLayer overhead.

4. Multi-NetCDF input: Tool 1 validates that all inputs share the same
   spatial grid (lat/lon arrays). Inconsistent grids = error, not warning,
   because mixing grids silently produces incorrect downstream results.

5. Polygon mask cache: same pattern as ecde-arcgis-tools. The mask is
   rasterised once per unique grid signature and reused across files.

6. Unit handling: Tool 1 preserves units in the manifest CSV. Unit
   conversion (e.g. Kelvin to Celsius) is responsibility of Tool 3 when
   computing temperature-specific metrics. Tool 1 outputs are raw values
   in original units.

Dependencies
------------
ArcGIS Pro 3.x default Python environment. Required libraries (all ship
with arcgispro-py3): arcpy, numpy, netCDF4, osgeo.gdal.

Author: Pedro Goncalves, LNEG / FCUP / 41 Norte, 2026.
License: Apache-2.0 (see LICENSE file).
"""

from __future__ import annotations

import csv
import os
import re
import time as _time
import warnings
from concurrent.futures import ThreadPoolExecutor

import arcpy
import numpy as np

try:
    import netCDF4
except ImportError as e:
    raise ImportError(
        "netCDF4 is required. It ships with ArcGIS Pro 3.x by default."
    ) from e

try:
    from osgeo import gdal, osr
    gdal.UseExceptions()
    # Enable GDAL's internal multi-threaded I/O for raster decompression
    # (DEFLATE/LZW workers). Speeds up batch reads of compressed daily TIFFs
    # roughly 1.5-2x on multi-core systems. Safe with read-only access.
    gdal.SetConfigOption("GDAL_NUM_THREADS", "ALL_CPUS")
except ImportError as e:
    raise ImportError(
        "osgeo (GDAL) is required. It ships with ArcGIS Pro 3.x by default."
    ) from e


# Tool 3 metrics aggregate over the time-axis with nan-aware NumPy
# reductions (nanmean, nanmin, nanmax, nanpercentile, ...). For pixels
# where the entire time-axis is NaN (e.g. ocean cells in the ERA5-Land
# land-only grid), NumPy correctly returns NaN but emits a RuntimeWarning
# per call. Over thousands of pixels in a 30-year batch this floods the
# ArcGIS Pro tool output with cosmetic warnings. Suppress at module load.
warnings.filterwarnings(
    "ignore", category=RuntimeWarning,
    message="All-NaN slice encountered",
)
warnings.filterwarnings(
    "ignore", category=RuntimeWarning,
    message="Mean of empty slice",
)
warnings.filterwarnings(
    "ignore", category=RuntimeWarning,
    message="invalid value encountered in",
)


# ===========================================================================
# Configuration
# ===========================================================================

# Candidate names for the time / lat / lon dimensions, in priority order.
# Lower-cased comparison is used. ERA5-Land may use 'valid_time' in some
# downloads.
TIME_DIM_CANDIDATES = ("time", "valid_time", "t", "stdtime")
LAT_DIM_CANDIDATES = ("latitude", "lat", "y")
LON_DIM_CANDIDATES = ("longitude", "lon", "x")

# Coordinate and bound variables to exclude from data-variable auto-detect.
COORD_LIKE_NAMES = {
    "time", "valid_time", "t", "stdtime",
    "latitude", "lat", "y",
    "longitude", "lon", "x",
    "time_bnds", "valid_time_bnds", "lat_bnds", "lon_bnds",
    "latitude_bnds", "longitude_bnds",
    "height", "level", "depth", "expver", "number", "realization",
    "crs", "spatial_ref",
}

# GeoTIFF creation options
GTIFF_OPTIONS = [
    "COMPRESS=DEFLATE",
    "PREDICTOR=2",
    "TILED=YES",
    "BLOCKXSIZE=256",
    "BLOCKYSIZE=256",
    "BIGTIFF=IF_SAFER",
]

NODATA_VAL = -9999.0

# Number of threads used by read_paths_into_stack() to read daily TIFFs
# in parallel. GDAL is thread-safe in read-only mode (each OpenEx returns
# a fresh dataset handle) and Python releases the GIL during the C-level
# I/O, so this scales near-linearly until disk bandwidth saturates.
# Override via env var ERA5LAND_READ_WORKERS=N (set to 1 to force serial).
try:
    _PARALLEL_READ_WORKERS = max(
        1, int(os.environ.get("ERA5LAND_READ_WORKERS", "8"))
    )
except ValueError:
    _PARALLEL_READ_WORKERS = 8

# Variable short-name classification (used by Tool 3 for metric applicability
# and automatic unit conversion). Short names follow ECMWF GRIB conventions
# and are what netCDF4 returns as variable names in ERA5-Land files.
TEMPERATURE_VARIABLES = {
    "t2m", "d2m", "skt",                  # air, dewpoint, skin
    "stl1", "stl2", "stl3", "stl4",       # soil temperature 4 levels
    "tsn",                                # snow temperature
}

PRECIPITATION_VARIABLES = {
    "tp",                                 # total precipitation
}

RADIATION_VARIABLES = {
    "ssrd", "strd", "ssr", "str",         # downward / net SW & LW radiation
}

WIND_VARIABLES = {
    "u10", "v10",                         # 10m wind components
}

# Variables in Kelvin that get auto-converted to Celsius by Tool 3 before
# any metric computation. Identical to TEMPERATURE_VARIABLES today but kept
# separate for future cases where a temperature var might not need conversion.
KELVIN_VARIABLES = set(TEMPERATURE_VARIABLES)

# Reverse map used by Tool 3 for category detection and warnings.
def variable_type_of(short_name):
    sn = str(short_name).lower()
    if sn in TEMPERATURE_VARIABLES:
        return "temperature"
    if sn in PRECIPITATION_VARIABLES:
        return "precipitation"
    if sn in RADIATION_VARIABLES:
        return "radiation"
    if sn in WIND_VARIABLES:
        return "wind"
    return "other"

# Polygon mask cache keyed by (variable shape signature)
_MASK_CACHE = {}
_REGION_LABEL_CACHE = {}


# ===========================================================================
# NetCDF structure inspection
# ===========================================================================

def find_dim(ds, candidates):
    """Find a dimension whose name matches one of candidates (case-insensitive).
    Returns the actual name in the dataset, or None."""
    dim_names = {d.lower(): d for d in ds.dimensions.keys()}
    for c in candidates:
        if c.lower() in dim_names:
            return dim_names[c.lower()]
    return None


def find_coord_var(ds, candidates):
    """Find a coordinate variable matching candidates (case-insensitive)."""
    var_names = {v.lower(): v for v in ds.variables.keys()}
    for c in candidates:
        if c.lower() in var_names:
            return var_names[c.lower()]
    return None


def detect_data_variable(ds, user_choice=None):
    """Find the primary data variable (the one we want to convert).

    Strategy:
      1. If user_choice is given and exists, use it.
      2. Else, find variables with >= 2 dimensions that are NOT coord-like.
      3. Prefer the one with the largest size (most data).
    """
    if user_choice:
        if user_choice not in ds.variables:
            raise ValueError(
                "Variable '{}' not found. Available: {}".format(
                    user_choice, list(ds.variables.keys())
                )
            )
        return user_choice

    candidates = []
    for vname, v in ds.variables.items():
        if vname.lower() in COORD_LIKE_NAMES:
            continue
        if len(v.dimensions) >= 2:
            candidates.append((vname, v.size))

    if not candidates:
        raise ValueError(
            "No data variable found. All variables: {}".format(
                list(ds.variables.keys())
            )
        )

    # Sort by size descending; return name of largest
    candidates.sort(key=lambda x: -x[1])
    return candidates[0][0]


def get_var_metadata(ds, var_name):
    """Return a dict of useful metadata for the variable."""
    v = ds.variables[var_name]
    return {
        "name": var_name,
        "standard_name": getattr(v, "standard_name", ""),
        "long_name": getattr(v, "long_name", ""),
        "units": getattr(v, "units", ""),
        "shape": tuple(v.shape),
        "dimensions": tuple(v.dimensions),
        "fill_value": getattr(v, "_FillValue", getattr(v, "missing_value", None)),
    }


def get_years_from_time(ds, time_dim):
    """Read the time dimension and return a list of (year, month, day) tuples
    plus the year for each timestep.

    Uses netCDF4.num2date for robust calendar handling.
    """
    if time_dim not in ds.variables:
        raise ValueError(
            "Time dimension '{}' has no corresponding variable. "
            "Cannot decode time values.".format(time_dim)
        )

    tvar = ds.variables[time_dim]
    units = getattr(tvar, "units", None)
    calendar = getattr(tvar, "calendar", "standard")

    if not units:
        raise ValueError(
            "Time variable '{}' has no 'units' attribute. Cannot decode.".format(
                time_dim
            )
        )

    times = netCDF4.num2date(
        tvar[:],
        units=units,
        calendar=calendar,
        only_use_cftime_datetimes=False,
        only_use_python_datetimes=False,
    )

    dates = []
    for t in times:
        # Both datetime.datetime and cftime objects have these attributes
        dates.append((int(t.year), int(t.month), int(t.day)))
    return dates


# ===========================================================================
# Filename-based daily statistic detection
# ===========================================================================

# Tokens recognised as the "short stat" segment of source NetCDF filenames
# produced by the era5land-downloads concat step:
#   era5land_<variable>_<short_stat>_<year>.nc           (annual)
#   era5land_<short_var>_<short_stat>_<year>-<MM>.nc     (monthly)
# short_stat = "daily_<token>".replace("daily_", "")
STATISTIC_TOKENS = {
    "mean", "maximum", "minimum", "sum", "stddev",
    "max", "min",  # tolerated short variants
}

# Marker used when the source filename does not match the expected pattern.
# Downstream tools (Tool 3) treat this as "statistic unknown, validation
# warnings may apply".
UNKNOWN_STATISTIC = "unknown_stat"


def detect_statistic_from_filename(nc_path):
    """Infer the daily statistic from a NetCDF source filename.

    Returns the normalised "daily_<token>" form (e.g. "daily_maximum"),
    or None if the pattern is not recognised. Callers typically substitute
    UNKNOWN_STATISTIC when None is returned.
    """
    stem = os.path.splitext(os.path.basename(nc_path))[0]
    parts = stem.split("_")
    if len(parts) < 4:
        return None
    last = parts[-1]
    if not (re.match(r"^\d{4}$", last) or re.match(r"^\d{4}-\d{2}$", last)):
        return None
    candidate = parts[-2].lower()
    if candidate not in STATISTIC_TOKENS:
        return None
    # Normalise short variants
    if candidate == "max":
        candidate = "maximum"
    elif candidate == "min":
        candidate = "minimum"
    return "daily_" + candidate


# ===========================================================================
# Grid utilities
# ===========================================================================

def grid_signature(lat, lon):
    """Hashable key identifying a grid by extent + size."""
    return (
        len(lat), len(lon),
        round(float(lat[0]), 6), round(float(lat[-1]), 6),
        round(float(lon[0]), 6), round(float(lon[-1]), 6),
    )


def grids_match(g1, g2, tolerance=1e-4):
    """Compare two grid signatures. Returns True if effectively identical."""
    if g1[0] != g2[0] or g1[1] != g2[1]:
        return False
    for a, b in zip(g1[2:], g2[2:]):
        if abs(a - b) > tolerance:
            return False
    return True


def compute_geotransform(lon, lat, lat_ascending):
    """GDAL geotransform for north-up output: (origin_x, dx, 0, origin_y, 0, -dy)."""
    cs_x = float(abs(lon[1] - lon[0]))
    cs_y = float(abs(lat[1] - lat[0]))
    origin_x = float(lon[0]) - cs_x / 2.0
    if lat_ascending:
        origin_y = float(lat[-1]) + cs_y / 2.0
    else:
        origin_y = float(lat[0]) + cs_y / 2.0
    return (origin_x, cs_x, 0.0, origin_y, 0.0, -cs_y)


def compute_bbox_slices(lat, lon, clip_extent, buffer_cells=1):
    """Return (lat_slice, lon_slice) that covers the clip extent bbox."""
    cs_x = abs(lon[1] - lon[0])
    cs_y = abs(lat[1] - lat[0])
    buf_x = buffer_cells * cs_x
    buf_y = buffer_cells * cs_y

    lon_mask = (lon >= clip_extent.XMin - buf_x) & (lon <= clip_extent.XMax + buf_x)
    lat_mask = (lat >= clip_extent.YMin - buf_y) & (lat <= clip_extent.YMax + buf_y)

    lon_idx = np.where(lon_mask)[0]
    lat_idx = np.where(lat_mask)[0]

    if len(lon_idx) == 0 or len(lat_idx) == 0:
        raise ValueError(
            "Clip extent does not intersect NetCDF grid. "
            "Clip: [{}, {}, {}, {}], NetCDF lon: [{}, {}], lat: [{}, {}]".format(
                clip_extent.XMin, clip_extent.YMin,
                clip_extent.XMax, clip_extent.YMax,
                lon[0], lon[-1], lat[0], lat[-1]
            )
        )

    return (slice(lat_idx[0], lat_idx[-1] + 1),
            slice(lon_idx[0], lon_idx[-1] + 1))


# ===========================================================================
# Mask rasterisation (same pattern as ecde-arcgis-tools v1.0.1)
# ===========================================================================

def rasterise_clip_mask(clip_fc, lat, lon, lat_ascending):
    """Rasterise the clip polygon to match the (sliced) lat/lon grid.

    Returns boolean numpy array, same orientation as data (lat-native).
    Uses a temporary SHORT field to avoid 64-bit OID issues.
    """
    cs_x = float(abs(lon[1] - lon[0]))
    cs_y = float(abs(lat[1] - lat[0]))

    left = float(lon[0]) - cs_x / 2.0
    right = float(lon[-1]) + cs_x / 2.0
    if lat_ascending:
        bottom = float(lat[0]) - cs_y / 2.0
        top = float(lat[-1]) + cs_y / 2.0
    else:
        bottom = float(lat[-1]) - cs_y / 2.0
        top = float(lat[0]) + cs_y / 2.0

    tmp = "in_memory/_era5_mask_temp"
    if arcpy.Exists(tmp):
        arcpy.management.Delete(tmp)

    tmp_field = "_era5_mask_one"
    source = clip_fc
    added_to_source = False
    try:
        existing = [f.name for f in arcpy.ListFields(clip_fc)]
        if tmp_field not in existing:
            try:
                arcpy.management.AddField(source, tmp_field, "SHORT")
                arcpy.management.CalculateField(source, tmp_field, "1", "PYTHON3")
                added_to_source = True
            except arcpy.ExecuteError:
                source = "in_memory/_era5_mask_src"
                if arcpy.Exists(source):
                    arcpy.management.Delete(source)
                arcpy.management.CopyFeatures(clip_fc, source)
                arcpy.management.AddField(source, tmp_field, "SHORT")
                arcpy.management.CalculateField(source, tmp_field, "1", "PYTHON3")

        with arcpy.EnvManager(
            extent="{} {} {} {}".format(left, bottom, right, top),
            cellSize=cs_x,
            snapRaster=None,
            outputCoordinateSystem=arcpy.SpatialReference(4326),
        ):
            arcpy.conversion.PolygonToRaster(
                in_features=source,
                value_field=tmp_field,
                out_rasterdataset=tmp,
                cell_assignment="CELL_CENTER",
                cellsize=cs_x,
            )
            mask_arr = arcpy.RasterToNumPyArray(tmp, nodata_to_value=-1)
            arcpy.management.Delete(tmp)
    finally:
        if added_to_source:
            try:
                arcpy.management.DeleteField(clip_fc, [tmp_field])
            except Exception:
                pass
        if source != clip_fc and arcpy.Exists(source):
            try:
                arcpy.management.Delete(source)
            except Exception:
                pass

    # arcpy returns north-up (row 0 = highest lat); flip to lat-native if ascending
    if lat_ascending:
        mask_arr = mask_arr[::-1, :]

    nrows_target, ncols_target = len(lat), len(lon)
    if mask_arr.shape != (nrows_target, ncols_target):
        result = np.full((nrows_target, ncols_target), -1, dtype=mask_arr.dtype)
        r = min(mask_arr.shape[0], nrows_target)
        c = min(mask_arr.shape[1], ncols_target)
        result[:r, :c] = mask_arr[:r, :c]
        mask_arr = result

    return mask_arr != -1


def get_or_compute_mask(clip_fc, lat, lon, lat_ascending):
    if clip_fc is None:
        return None
    sig = grid_signature(lat, lon)
    if sig not in _MASK_CACHE:
        _MASK_CACHE[sig] = rasterise_clip_mask(clip_fc, lat, lon, lat_ascending)
    return _MASK_CACHE[sig]


# ===========================================================================
# Regional zonal stats (Tool 4) helpers
# ===========================================================================

def rasterise_regions_to_labels(region_fc, region_id_field, lat, lon,
                                  lat_ascending):
    """Rasterise polygon FC to an int label grid where each pixel holds the
    value of region_id_field for the polygon that contains it. Pixels outside
    any polygon are set to a sentinel (-1).

    Pattern parallels rasterise_clip_mask but with an attribute-driven
    value_field instead of a constant-1 field. The region_id_field MUST be
    integer-typed (SHORT, LONG, or OID); a clear error is raised otherwise.
    OID values are copied to a LONG proxy to dodge the 64-bit OID/GRID
    restriction (same workaround as rasterise_clip_mask, design decision #7).

    Returns: 2D int32 array (lat-native orientation), shape (len(lat), len(lon)).
    """
    fields = {f.name: f.type for f in arcpy.ListFields(region_fc)}
    if region_id_field not in fields:
        raise ValueError(
            "Region ID field '{}' not in {}; available: {}".format(
                region_id_field, region_fc, sorted(fields.keys())
            )
        )
    field_type = fields[region_id_field]
    if field_type not in ("Integer", "SmallInteger", "OID"):
        raise ValueError(
            "Region ID field '{}' has type '{}'; must be integer "
            "(SHORT/LONG/OID). Add an int column to the FC or pick a "
            "different field.".format(region_id_field, field_type)
        )

    cs_x = float(abs(lon[1] - lon[0]))
    cs_y = float(abs(lat[1] - lat[0]))
    left = float(lon[0]) - cs_x / 2.0
    right = float(lon[-1]) + cs_x / 2.0
    if lat_ascending:
        bottom = float(lat[0]) - cs_y / 2.0
        top = float(lat[-1]) + cs_y / 2.0
    else:
        bottom = float(lat[-1]) - cs_y / 2.0
        top = float(lat[0]) + cs_y / 2.0

    tmp = "in_memory/_era5_region_labels"
    if arcpy.Exists(tmp):
        arcpy.management.Delete(tmp)

    # 64-bit OID -> LONG proxy. PolygonToRaster cannot write OBJECTID
    # values directly to a GRID (ERROR 003658). Same fix as Tool 1's
    # rasterise_clip_mask.
    source = region_fc
    use_field = region_id_field
    proxy_source = None
    proxy_field = "_era5_region_int"
    if field_type == "OID":
        proxy_source = "in_memory/_era5_region_src"
        if arcpy.Exists(proxy_source):
            arcpy.management.Delete(proxy_source)
        arcpy.management.CopyFeatures(region_fc, proxy_source)
        arcpy.management.AddField(proxy_source, proxy_field, "LONG")
        arcpy.management.CalculateField(
            proxy_source, proxy_field,
            "!{}!".format(region_id_field), "PYTHON3"
        )
        source = proxy_source
        use_field = proxy_field

    try:
        with arcpy.EnvManager(
            extent="{} {} {} {}".format(left, bottom, right, top),
            cellSize=cs_x,
            snapRaster=None,
            outputCoordinateSystem=arcpy.SpatialReference(4326),
        ):
            arcpy.conversion.PolygonToRaster(
                in_features=source,
                value_field=use_field,
                out_rasterdataset=tmp,
                cell_assignment="CELL_CENTER",
                cellsize=cs_x,
            )
            labels = arcpy.RasterToNumPyArray(tmp, nodata_to_value=-1)
            arcpy.management.Delete(tmp)
    finally:
        if proxy_source is not None and arcpy.Exists(proxy_source):
            try:
                arcpy.management.Delete(proxy_source)
            except Exception:
                pass

    # arcpy returns north-up; flip to lat-native if ascending
    if lat_ascending:
        labels = labels[::-1, :]

    nrows_t, ncols_t = len(lat), len(lon)
    if labels.shape != (nrows_t, ncols_t):
        result = np.full((nrows_t, ncols_t), -1, dtype=np.int32)
        r = min(labels.shape[0], nrows_t)
        c = min(labels.shape[1], ncols_t)
        result[:r, :c] = labels[:r, :c]
        labels = result

    return labels.astype(np.int32)


def get_or_compute_region_labels(region_fc, region_id_field, lat, lon,
                                   lat_ascending):
    """Cached wrapper. Cache key includes the field name so different
    fields on the same FC + grid get distinct entries."""
    if region_fc is None or region_id_field is None:
        return None
    sig = (grid_signature(lat, lon), str(region_fc), str(region_id_field))
    if sig not in _REGION_LABEL_CACHE:
        _REGION_LABEL_CACHE[sig] = rasterise_regions_to_labels(
            region_fc, region_id_field, lat, lon, lat_ascending
        )
    return _REGION_LABEL_CACHE[sig]


def compute_zonal_stats(data, labels, stat_names):
    """Vectorised zonal statistics per region using numpy bincount.

    Args:
        data: 2D float array (one TIFF). NaN = invalid pixel.
        labels: 2D int array same shape as data. -1 = no region.
        stat_names: iterable from {'mean','min','max','std','count'}.

    Returns dict {region_id: {stat_name: value}}. Regions with zero
    valid pixels are omitted.
    """
    if data.shape != labels.shape:
        raise ValueError("data and labels must have the same shape; got {} vs {}"
                         .format(data.shape, labels.shape))
    stat_names = list(stat_names)
    flat_d = data.ravel()
    flat_l = labels.ravel()
    valid = ~np.isnan(flat_d) & (flat_l >= 0)
    if not np.any(valid):
        return {}

    vd = flat_d[valid].astype(np.float64)
    vl = flat_l[valid].astype(np.int64)
    max_lbl = int(vl.max()) + 1

    counts = np.bincount(vl, minlength=max_lbl)

    means = None
    stds = None
    mins = None
    maxs = None

    if "mean" in stat_names or "std" in stat_names:
        sums = np.bincount(vl, weights=vd, minlength=max_lbl)
        means = np.where(counts > 0, sums / np.maximum(counts, 1), np.nan)

    if "std" in stat_names:
        sq_sums = np.bincount(vl, weights=vd * vd, minlength=max_lbl)
        variances = np.where(
            counts > 0,
            (sq_sums / np.maximum(counts, 1)) - means * means,
            np.nan,
        )
        # Floor at 0 to avoid sqrt of small negative from float roundoff
        variances = np.where(np.isnan(variances), np.nan, np.maximum(variances, 0.0))
        stds = np.sqrt(variances)

    if "min" in stat_names:
        mins = np.full(max_lbl, np.inf)
        np.minimum.at(mins, vl, vd)
        mins = np.where(counts > 0, mins, np.nan)

    if "max" in stat_names:
        maxs = np.full(max_lbl, -np.inf)
        np.maximum.at(maxs, vl, vd)
        maxs = np.where(counts > 0, maxs, np.nan)

    out = {}
    for rid in range(max_lbl):
        if counts[rid] == 0:
            continue
        entry = {}
        if "mean" in stat_names:
            entry["mean"] = float(means[rid])
        if "min" in stat_names:
            entry["min"] = float(mins[rid])
        if "max" in stat_names:
            entry["max"] = float(maxs[rid])
        if "std" in stat_names:
            entry["std"] = float(stds[rid])
        if "count" in stat_names:
            entry["count"] = int(counts[rid])
        out[rid] = entry

    return out


# ===========================================================================
# GeoTIFF writer
# ===========================================================================

def write_geotiff_gdal(arr, out_path, geotransform, srs_wkt, nodata=NODATA_VAL):
    """Write a 2D NumPy array as a compressed GeoTIFF via GDAL."""
    driver = gdal.GetDriverByName("GTiff")
    nrows, ncols = arr.shape
    ds_out = driver.Create(
        out_path, ncols, nrows, 1, gdal.GDT_Float32, options=GTIFF_OPTIONS
    )
    ds_out.SetGeoTransform(geotransform)
    ds_out.SetProjection(srs_wkt)
    band = ds_out.GetRasterBand(1)
    band.WriteArray(arr.astype(np.float32))
    band.SetNoDataValue(float(nodata))
    band.FlushCache()
    ds_out = None


# ===========================================================================
# Climate indicators (Tool 3) - I/O helpers
# ===========================================================================

def read_tif_to_array(tif_path):
    """Read a TIFF as a float32 NumPy array with NoData mapped to NaN.

    Uses GDAL directly (consistency with write_geotiff_gdal). The input is
    assumed north-up (the orientation written by Tool 1). NoData -> NaN is
    done in-place to avoid an extra full-array allocation per file.
    """
    ds = gdal.OpenEx(tif_path, gdal.OF_RASTER | gdal.OF_READONLY)
    if ds is None:
        raise RuntimeError("Failed to open TIFF: {}".format(tif_path))
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray().astype(np.float32, copy=False)
    nd = band.GetNoDataValue()
    ds = None
    if not arr.flags.writeable:
        arr = arr.copy()
    if nd is not None:
        ndf = np.float32(nd)
        # In-place NaN substitution avoids a second full-array copy
        np.place(arr, arr == ndf, np.nan)
    np.place(arr, arr == np.float32(NODATA_VAL), np.nan)
    return arr


def read_paths_into_stack(paths):
    """Load N daily TIFFs into a single pre-allocated stack (N, lat, lon).

    Streams each TIFF read directly into its slice of the destination
    buffer (peak memory ~1x stack size, versus ~2x with np.stack of a
    list comprehension) and parallelises the reads via ThreadPoolExecutor
    when more than one worker is configured (see _PARALLEL_READ_WORKERS).

    Slice assignments `stack[i] = ...` write to disjoint memory regions
    and are thread-safe in NumPy. GDAL read-only access via OpenEx is
    thread-safe per dataset.
    """
    n = len(paths)
    if n == 0:
        raise ValueError("read_paths_into_stack: empty paths list")
    sample = read_tif_to_array(paths[0])
    stack = np.empty((n,) + sample.shape, dtype=sample.dtype)
    stack[0] = sample
    del sample

    if n == 1:
        return stack

    remaining_indexed = list(enumerate(paths[1:], start=1))
    n_workers = min(_PARALLEL_READ_WORKERS, len(remaining_indexed))

    if n_workers <= 1:
        for i, p in remaining_indexed:
            stack[i] = read_tif_to_array(p)
        return stack

    def _read_into_slot(idx_path):
        i, p = idx_path
        stack[i] = read_tif_to_array(p)

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        # Consume the iterator so any per-task exception propagates.
        for _ in executor.map(_read_into_slot, remaining_indexed):
            pass
    return stack


def read_grid_metadata(tif_path):
    """Return (geotransform, srs_wkt, nrows, ncols) for a TIFF."""
    ds = gdal.OpenEx(tif_path, gdal.OF_RASTER | gdal.OF_READONLY)
    if ds is None:
        raise RuntimeError("Failed to open TIFF: {}".format(tif_path))
    gt = ds.GetGeoTransform()
    srs_wkt = ds.GetProjection()
    nrows = ds.RasterYSize
    ncols = ds.RasterXSize
    ds = None
    return gt, srs_wkt, nrows, ncols


def derive_lat_lon_from_gt(gt, nrows, ncols):
    """Build (lat, lon) cell-centre arrays from a north-up geotransform.

    For TIFFs written by Tool 1 (always north-up), gt[5] is negative and the
    resulting lat array is descending (row 0 = highest lat). lon is ascending.
    """
    origin_x, dx, _, origin_y, _, dy = gt
    cs_x = abs(dx)
    lon = origin_x + cs_x / 2.0 + np.arange(ncols, dtype=np.float64) * cs_x
    lat = origin_y + dy / 2.0 + np.arange(nrows, dtype=np.float64) * dy
    return lat, lon


def discover_inputs_from_manifest(out_root):
    """Read Tool 1 manifest CSV and return inputs index.

    Returns: ({index, units}, None) on success, (None, error_message) on
    failure. The index maps (variable, statistic) -> {year: [sorted paths]}.
    """
    manifest_path = os.path.join(out_root, "manifest.csv")
    if not os.path.isfile(manifest_path):
        return None, "manifest.csv not found at {}".format(manifest_path)

    index = {}
    units_by_var = {}
    try:
        with open(manifest_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            required = {"raster", "variable", "year"}
            missing = required - set(reader.fieldnames or [])
            if missing:
                return None, "manifest missing required columns: {}".format(missing)
            for row in reader:
                var = row.get("variable", "")
                stat = row.get("statistic") or UNKNOWN_STATISTIC
                try:
                    year = int(row.get("year", "0"))
                except ValueError:
                    continue
                path = row.get("raster", "")
                if not (var and path):
                    continue
                index.setdefault((var, stat), {}).setdefault(year, []).append(path)
                if var not in units_by_var:
                    units_by_var[var] = row.get("units", "")
    except (OSError, csv.Error) as e:
        return None, "manifest read error: {}".format(e)

    for k in index:
        for yr in index[k]:
            index[k][yr].sort()
    return {"index": index, "units": units_by_var}, None


def discover_t4_inputs(out_root):
    """Auto-detect Tool 1 vs Tool 3 manifest and return a uniform list of
    input items for Tool 4 to iterate.

    Detection rule: if manifest.csv has a 'metric' column it is a Tool 3
    output (annual indicators / climatology / IPMA normals); otherwise it
    is a Tool 1 output (per-day raw values).

    Returns (items, source, error) where:
      items: list of dicts with uniform keys
        - raster (str path)
        - variable (str)
        - statistic (str)
        - metric (str or '')
        - period_kind (str): 'daily' for Tool 1; 'annual', 'climatology',
                              'ipma_normals_1981_2010', etc. for Tool 3
        - year (int or None)
        - start_year (int or None)  -- only Tool 3 climatology/normals
        - end_year (int or None)    -- ditto
        - date (str or '')          -- Tool 1 only, YYYY-MM-DD
        - units (str)
      source: 'tool1' or 'tool3'
      error: None on success, message on failure
    """
    manifest_path = os.path.join(out_root, "manifest.csv")
    if not os.path.isfile(manifest_path):
        return [], None, "manifest.csv not found at {}".format(manifest_path)

    try:
        with open(manifest_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            headers = set(reader.fieldnames or [])
            source = "tool3" if "metric" in headers else "tool1"
            items = []
            for row in reader:
                raster = (row.get("raster") or "").strip()
                if not raster or not os.path.isfile(raster):
                    continue
                var = (row.get("variable") or "").strip()
                stat = (row.get("statistic") or "").strip() or UNKNOWN_STATISTIC
                units = (row.get("units") or "").strip()
                if source == "tool3":
                    metric = (row.get("metric") or "").strip()
                    period_kind = (row.get("period_kind") or "annual").strip()
                    year = _safe_int(row.get("year"))
                    start_year = _safe_int(row.get("start_year"))
                    end_year = _safe_int(row.get("end_year"))
                    items.append({
                        "raster": raster, "variable": var, "statistic": stat,
                        "metric": metric, "period_kind": period_kind,
                        "year": year, "start_year": start_year,
                        "end_year": end_year, "date": "", "units": units,
                    })
                else:
                    year = _safe_int(row.get("year"))
                    date = (row.get("date") or "").strip()
                    items.append({
                        "raster": raster, "variable": var, "statistic": stat,
                        "metric": "", "period_kind": "daily",
                        "year": year, "start_year": None, "end_year": None,
                        "date": date, "units": units,
                    })
    except (OSError, csv.Error) as e:
        return [], None, "manifest read error: {}".format(e)

    return items, source, None


def _safe_int(value):
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def discover_inputs_from_walk(out_root):
    """Walk filesystem assuming <var>/<stat>/<year>/<file>.tif layout."""
    if not os.path.isdir(out_root):
        return None, "input folder does not exist: {}".format(out_root)

    index = {}
    for var in sorted(os.listdir(out_root)):
        var_dir = os.path.join(out_root, var)
        if not os.path.isdir(var_dir):
            continue
        for stat in sorted(os.listdir(var_dir)):
            stat_dir = os.path.join(var_dir, stat)
            if not os.path.isdir(stat_dir):
                continue
            for yr in sorted(os.listdir(stat_dir)):
                year_dir = os.path.join(stat_dir, yr)
                if not os.path.isdir(year_dir):
                    continue
                try:
                    year_int = int(yr)
                except ValueError:
                    continue
                paths = sorted(
                    os.path.join(year_dir, f)
                    for f in os.listdir(year_dir)
                    if f.lower().endswith(".tif")
                )
                if paths:
                    index.setdefault((var, stat), {})[year_int] = paths
    return {"index": index, "units": {}}, None


def _audit_index_against_disk(index, log):
    """Remove entries from the inputs index whose .tif files do not exist.

    Defends against the manifest being out of sync with disk reality
    (e.g. user deleted some Tool 1 outputs between runs, or OneDrive
    offloaded files marked as "online only"). Modifies index in place.
    Returns (missing_files_count, removed_keys).
    """
    missing_files = 0
    keys_removed = []
    for key in list(index.keys()):
        empty_years = []
        for year, paths in list(index[key].items()):
            existing = [p for p in paths if os.path.isfile(p)]
            if len(existing) < len(paths):
                missing_files += len(paths) - len(existing)
            if existing:
                index[key][year] = existing
            else:
                empty_years.append(year)
        for y in empty_years:
            del index[key][y]
        if not index[key]:
            keys_removed.append(key)
    for k in keys_removed:
        del index[k]
    return missing_files, keys_removed


def build_inputs_index(out_root, log):
    """Build the inputs index, preferring manifest over filesystem walk.

    After loading, audits each registered .tif path against the filesystem
    and drops missing entries (with warnings). Prevents mid-run crashes
    from manifest/disk drift.
    """
    result, err = discover_inputs_from_manifest(out_root)
    used_manifest = bool(result and result["index"])
    if used_manifest:
        log("  inputs: loaded from manifest.csv ({} (var,stat) groups)".format(
            len(result["index"])
        ))
    else:
        if err:
            log("  manifest unavailable ({}), trying filesystem walk".format(err))
        result, err = discover_inputs_from_walk(out_root)
        if err:
            raise RuntimeError(err)
        if not result["index"]:
            raise RuntimeError(
                "No (variable, statistic, year) inputs found under {}".format(out_root)
            )
        log("  inputs: discovered via filesystem walk ({} (var,stat) groups)".format(
            len(result["index"])
        ))

    # Audit only matters for manifest-loaded index (walk is filesystem-truth
    # by construction). Still cheap to run on the walk result; defensive.
    missing, removed = _audit_index_against_disk(result["index"], log)
    if missing:
        log("  audit: {} TIFF path(s) listed in manifest are missing on disk; "
            "those entries dropped.".format(missing), "warn")
    if removed:
        log("  audit: {} (variable, statistic) group(s) fully dropped: {}".format(
            len(removed), sorted(removed)
        ), "warn")
    if not result["index"]:
        raise RuntimeError(
            "After audit, no usable (variable, statistic, year) inputs remain "
            "under {}. Re-run Tool 1 or restore missing files.".format(out_root)
        )
    return result


# ===========================================================================
# Climate indicators (Tool 3) - Metric registry
# ===========================================================================
#
# Single source of truth for metric definitions. Each entry maps a metric
# name to a spec dict with these keys:
#
#   category:           one of "generic", "temperature", "precipitation",
#                       "radiation", "wind", "dtr".
#   label:              human-readable name shown in the multi-select UI.
#   applicable_to:      set of variable types ("temperature", "precipitation",
#                       "radiation", "wind", "other") for which this metric
#                       is meaningful. None means "any type".
#   expects_statistic:  recommended input statistic ("daily_mean",
#                       "daily_maximum", "daily_minimum", "daily_sum"), or
#                       None. Used only for validation warnings.
#   params:             tuple of param keys the metric reads from the UI
#                       (e.g. ("hdd_base",)). Empty tuple = no params.
#   compute:            callable(stack, **params) -> 2D NaN-aware array.
#                       stack shape (n_days, lat, lon), already in working
#                       units (Celsius for Kelvin vars, raw otherwise).
#   units_out:          str literal OR callable(working_units) -> str.
#   climatology:        aggregator name for combining annual values into
#                       climatology / IPMA normals. Default "mean".
#
# Increment 1 registers only "mean" as a placeholder; further metrics arrive
# in subsequent increments without touching this skeleton.

def _passthrough_units(u):
    return u or "unknown"


def _silent_log(_msg, _level="msg"):
    """No-op logger used to suppress duplicate messages from helper calls
    that run twice within one logical operation (e.g. DTR converting both
    Tmax and Tmin K to degC)."""
    return None


def _nan_if_all_invalid(result, stack):
    """Force result to NaN where the entire time-axis of stack is NaN.
    Avoids spurious zeros (from nansum / count) over fully-masked pixels."""
    all_invalid = ~np.any(~np.isnan(stack), axis=0)
    return np.where(all_invalid, np.nan, result)


def _max_consecutive_run(bool_stack):
    """Longest consecutive run of True along axis 0, per pixel.

    NaN in the source data evaluates to False in the comparison that builds
    bool_stack, so NaN naturally breaks runs. Caller may post-process with
    _nan_if_all_invalid to distinguish "0 days above" from "no data".
    """
    n_t = bool_stack.shape[0]
    if n_t == 0:
        return np.zeros(bool_stack.shape[1:], dtype=np.float32)
    current = np.zeros(bool_stack.shape[1:], dtype=np.int32)
    max_run = np.zeros(bool_stack.shape[1:], dtype=np.int32)
    for t in range(n_t):
        current = np.where(bool_stack[t], current + 1, 0)
        max_run = np.maximum(max_run, current)
    return max_run.astype(np.float32)


# Generic metric implementations (NaN-aware; valid-mask preserved)

def _gen_sum(s):
    return _nan_if_all_invalid(np.nansum(s, axis=0), s)

def _gen_count_above(s, generic_threshold):
    valid = ~np.isnan(s)
    counts = np.sum((s > generic_threshold) & valid, axis=0).astype(np.float32)
    return _nan_if_all_invalid(counts, s)

def _gen_count_below(s, generic_threshold):
    valid = ~np.isnan(s)
    counts = np.sum((s < generic_threshold) & valid, axis=0).astype(np.float32)
    return _nan_if_all_invalid(counts, s)

def _gen_seconds_above(s, generic_threshold):
    return _gen_count_above(s, generic_threshold) * 86400.0

def _gen_seconds_below(s, generic_threshold):
    return _gen_count_below(s, generic_threshold) * 86400.0

def _gen_accum_above(s, generic_threshold):
    above = np.where(s > generic_threshold, s, 0.0)
    return _nan_if_all_invalid(np.nansum(above, axis=0), s)

def _gen_accum_below(s, generic_threshold):
    below = np.where(s < generic_threshold, s, 0.0)
    return _nan_if_all_invalid(np.nansum(below, axis=0), s)

def _gen_max_consec_above(s, generic_threshold):
    return _nan_if_all_invalid(
        _max_consecutive_run(s > generic_threshold), s
    )

def _gen_max_consec_below(s, generic_threshold):
    return _nan_if_all_invalid(
        _max_consecutive_run(s < generic_threshold), s
    )


# Temperature-specific metric implementations.
# Stack is in degC (Kelvin conversion already applied by apply_kelvin_conversion).
# Thresholds (frost/ice/summer/tropical/hot) follow ETCCDI standards and are
# hard-coded; bases for HDD/CDD/GDD are user-configurable.

def _temp_hdd_spinoni(s, hdd_base):
    """Heating Degree Days, Spinoni method.
    HDD = sum over the period of max(hdd_base - T_mean, 0), in degC*day.
    """
    diff = np.where(s < hdd_base, hdd_base - s, 0.0)
    return _nan_if_all_invalid(np.nansum(diff, axis=0), s)


def _temp_cdd_spinoni(s, cdd_base):
    """Cooling Degree Days, Spinoni method.
    CDD = sum over the period of max(T_mean - cdd_base, 0), in degC*day.
    """
    diff = np.where(s > cdd_base, s - cdd_base, 0.0)
    return _nan_if_all_invalid(np.nansum(diff, axis=0), s)


def _temp_gdd(s, gdd_base):
    """Growing Degree Days.
    GDD = sum over the period of max(T_mean - gdd_base, 0), in degC*day.
    """
    diff = np.where(s > gdd_base, s - gdd_base, 0.0)
    return _nan_if_all_invalid(np.nansum(diff, axis=0), s)


def _temp_frost_days(s):
    """Days with Tmin < 0 degC (ETCCDI FD).

    Numerically identical to _temp_ice_days (both count days below 0 degC).
    The difference is semantic: this metric expects daily-min input and
    counts the days where the COLDEST hour fell below freezing. Tool 3
    emits a statistic-mismatch warning if daily_minimum is not supplied.
    """
    return _gen_count_below(s, 0.0)


def _temp_ice_days(s):
    """Days with Tmax < 0 degC (ETCCDI ID).

    Numerically identical to _temp_frost_days (both count days below 0 degC).
    The difference is semantic: this metric expects daily-max input and
    counts the days where even the WARMEST hour stayed below freezing.
    Tool 3 emits a statistic-mismatch warning if daily_maximum is not supplied.
    """
    return _gen_count_below(s, 0.0)


def _temp_summer_days(s):
    """Days with Tmax > 25 degC (ETCCDI SU)."""
    return _gen_count_above(s, 25.0)


def _temp_tropical_nights(s):
    """Nights with Tmin > 20 degC (ETCCDI TR)."""
    return _gen_count_above(s, 20.0)


def _temp_hot_days(s):
    """Days with Tmax > 30 degC."""
    return _gen_count_above(s, 30.0)


def _temp_heating_season_length(s, heating_threshold):
    """Heating Season Length: count of days with Tmean < heating_threshold.

    Distinct from frost_days/ice_days which use Tmin/Tmax against a fixed
    0 degC threshold; HSL uses daily-mean temperature with a user-tunable
    threshold typically aligned with HDD base (Spinoni default 15.5 degC).

    Designed as direct input for shallow geothermal potential mapping
    (G.POT method, Casasso & Sethi 2016): the heating-season length tc
    parameter feeds the empirical correlation for borehole heat exchanger
    sizing.
    """
    return _gen_count_below(s, heating_threshold)


def _temp_cooling_season_length(s, cooling_threshold):
    """Cooling Season Length: count of days with Tmean > cooling_threshold.

    Distinct from summer_days/hot_days/tropical_nights which use Tmax/Tmin
    against fixed thresholds; CSL uses daily-mean temperature with a
    user-tunable threshold typically aligned with CDD base (Spinoni
    default 22 degC).

    Designed as direct input for shallow geothermal potential mapping
    (G.POT method, Casasso & Sethi 2016) in cooling-mode scenarios.
    """
    return _gen_count_above(s, cooling_threshold)


def _temp_txn(s):
    """Annual minimum of daily-max temperature (ETCCDI TXn)."""
    return np.nanmin(s, axis=0)


def _temp_tnx(s):
    """Annual maximum of daily-min temperature (ETCCDI TNx)."""
    return np.nanmax(s, axis=0)


# Radiation-specific metric implementations.
#
# CRITICAL ASSUMPTION: the daily stack is in J/m^2/day where each daily
# value represents the TOTAL accumulated energy per square metre during
# that day. ERA5-Land hourly `ssrd` is "accumulated since 00:00 UTC"; the
# correct CDS statistic that yields the daily total is `daily_sum` (only
# available from the `derived-era5-land-daily-aggregations` dataset, NOT
# from `derived-era5-land-daily-statistics`). If the user instead provides
# daily_mean from `derived-era5-land-daily-statistics`, the output is
# scaled wrongly by a factor that depends on how CDS handles decumulation
# for that dataset (commonly ~24). See docs/metrics.md for guidance.
#
# Tool 3 emits a warning when a user selects a radiation metric with a
# statistic other than "daily_sum" via the standard expects_statistic
# validation path.

def _rad_w_per_m2_mean(s):
    """Annual mean radiation power in W/m^2.
    Input assumption: daily-total J/m^2/day. Annual mean of (J/day) divided
    by seconds-per-day (86400) yields the mean power in W/m^2."""
    return np.nanmean(s, axis=0) / 86400.0


def _rad_kwh_per_m2_year(s):
    """Annual radiation energy in kWh/m^2.
    Input assumption: daily-total J/m^2/day. 1 kWh = 3.6e6 J, so the
    annual sum divided by 3.6e6 yields kWh/m^2/year."""
    return _nan_if_all_invalid(np.nansum(s, axis=0) / 3.6e6, s)


def _rad_pv_potential(s, pv_performance_ratio):
    """Photovoltaic energy potential = kWh/m^2/year * PR.
    Input assumption: daily-total J/m^2/day. PR (performance ratio)
    accounts for inverter, soiling, temperature losses; default 0.75 is
    a typical value for residential PV."""
    energy = _nan_if_all_invalid(np.nansum(s, axis=0) / 3.6e6, s)
    return energy * pv_performance_ratio


def _rad_above_threshold_days(s, radiation_threshold_w_per_m2):
    """Days with daily mean radiation power above threshold (W/m^2).
    Input assumption: daily-total J/m^2/day. The threshold (W/m^2) is
    converted to J/m^2/day via x 86400 before comparing against the stack."""
    threshold_j_per_m2_day = float(radiation_threshold_w_per_m2) * 86400.0
    return _gen_count_above(s, threshold_j_per_m2_day)


# Wind-specific metric implementations.
# Working stack is the daily wind-speed magnitude in m/s, derived from u10
# and v10 components by process_year_wind() before any of these are called.

def _wind_power_density(s, wind_air_density):
    """Wind power density = 0.5 * rho * <v^3>, in W/m^2.
    Simple instantaneous-mean version; does not assume a Weibull fit."""
    return 0.5 * float(wind_air_density) * np.nanmean(s ** 3, axis=0)


def _wind_calm_days(s, wind_calm_threshold):
    """Days with mean wind speed below threshold (default 2 m/s = light air)."""
    return _gen_count_below(s, wind_calm_threshold)


# Precipitation-specific metric implementations (working units: mm/day after
# auto-conversion from m). Threshold default 1 mm follows ETCCDI conventions.

def _precip_total(s):
    """Annual total precipitation (sum over days). Units: mm."""
    return _nan_if_all_invalid(np.nansum(s, axis=0), s)


def _precip_wet_days(s, precip_wet_threshold):
    """Days with precip >= threshold (ETCCDI Rxxmm family)."""
    valid = ~np.isnan(s)
    counts = np.sum((s >= precip_wet_threshold) & valid, axis=0).astype(np.float32)
    return _nan_if_all_invalid(counts, s)


def _precip_consecutive_dry_days(s, precip_wet_threshold):
    """Max consecutive days with precip < threshold (ETCCDI CDD index)."""
    return _nan_if_all_invalid(
        _max_consecutive_run(s < precip_wet_threshold), s
    )


def _precip_consecutive_wet_days(s, precip_wet_threshold):
    """Max consecutive days with precip >= threshold (ETCCDI CWD index)."""
    return _nan_if_all_invalid(
        _max_consecutive_run(s >= precip_wet_threshold), s
    )


def _precip_max_1day(s):
    """Maximum daily precipitation (ETCCDI Rx1day)."""
    return np.nanmax(s, axis=0)


def _temp_gsl(s):
    """Growing Season Length (ETCCDI GSL, Northern Hemisphere variant).

    Start: index of first window of 6 consecutive days with T > 5 degC.
    End:   index (after Julian day 182) of first window of 6 consecutive
           days with T < 5 degC. If no end found, season extends to end of
           year (= n_days).
    GSL = end - start (days). Zero if no start window exists.
    """
    n_t = s.shape[0]
    if n_t < 6:
        return np.full(s.shape[1:], np.nan, dtype=np.float32)

    window = 6
    above = (s > 5.0).astype(np.int32)  # NaN > 5 is False
    below = (s < 5.0).astype(np.int32)
    csa = np.cumsum(above, axis=0)
    csb = np.cumsum(below, axis=0)
    n_w = n_t - window + 1

    win_above = np.empty((n_w,) + s.shape[1:], dtype=np.int32)
    win_above[0] = csa[window - 1]
    win_above[1:] = csa[window:] - csa[:n_w - 1]

    win_below = np.empty((n_w,) + s.shape[1:], dtype=np.int32)
    win_below[0] = csb[window - 1]
    win_below[1:] = csb[window:] - csb[:n_w - 1]

    all_above = (win_above == window)
    all_below = (win_below == window)

    # Start: first window with all > 5 degC anywhere in the year.
    any_start = np.any(all_above, axis=0)
    start_day = np.argmax(all_above, axis=0).astype(np.float32)

    # End: first all-below window with start index >= Julian day 182 (idx 181).
    july1_idx = 181
    end_search = all_below.copy()
    if july1_idx < n_w:
        end_search[:july1_idx] = False
    else:
        end_search[:] = False
    any_end = np.any(end_search, axis=0)
    end_day = np.argmax(end_search, axis=0).astype(np.float32)
    end_day = np.where(any_end, end_day, float(n_t))

    gsl = end_day - start_day
    # If no start window, GSL = 0 (no growing season this year)
    gsl = np.where(any_start, gsl, 0.0)
    return _nan_if_all_invalid(gsl, s)


METRIC_REGISTRY = {
    # ---- Generic: agnostic statistics ----
    "mean": {
        "category": "generic",
        "label": "Mean",
        "applicable_to": None,
        "expects_statistic": None,
        "params": (),
        "compute": lambda s: np.nanmean(s, axis=0),
        "units_out": _passthrough_units,
        "climatology": "mean",
    },
    "min": {
        "category": "generic", "label": "Minimum",
        "applicable_to": None, "expects_statistic": None, "params": (),
        "compute": lambda s: np.nanmin(s, axis=0),
        "units_out": _passthrough_units, "climatology": "mean",
    },
    "max": {
        "category": "generic", "label": "Maximum",
        "applicable_to": None, "expects_statistic": None, "params": (),
        "compute": lambda s: np.nanmax(s, axis=0),
        "units_out": _passthrough_units, "climatology": "mean",
    },
    "range": {
        "category": "generic", "label": "Range (max-min)",
        "applicable_to": None, "expects_statistic": None, "params": (),
        "compute": lambda s: np.nanmax(s, axis=0) - np.nanmin(s, axis=0),
        "units_out": _passthrough_units, "climatology": "mean",
    },
    "std": {
        "category": "generic", "label": "Standard deviation",
        "applicable_to": None, "expects_statistic": None, "params": (),
        "compute": lambda s: np.nanstd(s, axis=0),
        "units_out": _passthrough_units, "climatology": "mean",
    },
    "sum": {
        "category": "generic", "label": "Sum (annual total)",
        "applicable_to": None, "expects_statistic": None, "params": (),
        "compute": _gen_sum,
        "units_out": _passthrough_units, "climatology": "mean",
    },
    "median": {
        "category": "generic", "label": "Median",
        "applicable_to": None, "expects_statistic": None, "params": (),
        "compute": lambda s: np.nanmedian(s, axis=0),
        "units_out": _passthrough_units, "climatology": "mean",
    },
    "p10": {
        "category": "generic", "label": "10th percentile",
        "applicable_to": None, "expects_statistic": None, "params": (),
        "compute": lambda s: np.nanpercentile(s, 10, axis=0),
        "units_out": _passthrough_units, "climatology": "mean",
    },
    "p90": {
        "category": "generic", "label": "90th percentile",
        "applicable_to": None, "expects_statistic": None, "params": (),
        "compute": lambda s: np.nanpercentile(s, 90, axis=0),
        "units_out": _passthrough_units, "climatology": "mean",
    },

    # ---- Generic: threshold-based ----
    "count_above_threshold": {
        "category": "generic", "label": "Count days above threshold",
        "applicable_to": None, "expects_statistic": None,
        "params": ("generic_threshold",),
        "compute": _gen_count_above,
        "units_out": "days", "climatology": "mean",
    },
    "count_below_threshold": {
        "category": "generic", "label": "Count days below threshold",
        "applicable_to": None, "expects_statistic": None,
        "params": ("generic_threshold",),
        "compute": _gen_count_below,
        "units_out": "days", "climatology": "mean",
    },
    "seconds_above_threshold": {
        "category": "generic", "label": "Seconds above threshold (count*86400)",
        "applicable_to": None, "expects_statistic": None,
        "params": ("generic_threshold",),
        "compute": _gen_seconds_above,
        "units_out": "s", "climatology": "mean",
    },
    "seconds_below_threshold": {
        "category": "generic", "label": "Seconds below threshold (count*86400)",
        "applicable_to": None, "expects_statistic": None,
        "params": ("generic_threshold",),
        "compute": _gen_seconds_below,
        "units_out": "s", "climatology": "mean",
    },
    "accumulation_above_threshold": {
        "category": "generic", "label": "Sum of values above threshold",
        "applicable_to": None, "expects_statistic": None,
        "params": ("generic_threshold",),
        "compute": _gen_accum_above,
        "units_out": _passthrough_units, "climatology": "mean",
    },
    "accumulation_below_threshold": {
        "category": "generic", "label": "Sum of values below threshold",
        "applicable_to": None, "expects_statistic": None,
        "params": ("generic_threshold",),
        "compute": _gen_accum_below,
        "units_out": _passthrough_units, "climatology": "mean",
    },
    "max_consecutive_above": {
        "category": "generic", "label": "Max consecutive days above threshold",
        "applicable_to": None, "expects_statistic": None,
        "params": ("generic_threshold",),
        "compute": _gen_max_consec_above,
        "units_out": "days", "climatology": "mean",
    },
    "max_consecutive_below": {
        "category": "generic", "label": "Max consecutive days below threshold",
        "applicable_to": None, "expects_statistic": None,
        "params": ("generic_threshold",),
        "compute": _gen_max_consec_below,
        "units_out": "days", "climatology": "mean",
    },

    # ---- Temperature-specific (working units: degC) ----
    "hdd_spinoni": {
        "category": "temperature",
        "label": "Heating Degree Days (Spinoni)",
        "applicable_to": {"temperature"},
        "expects_statistic": "daily_mean",
        "params": ("hdd_base",),
        "compute": _temp_hdd_spinoni,
        "units_out": lambda u: (u or "degC") + "*day",
        "climatology": "mean",
    },
    "cdd_spinoni": {
        "category": "temperature",
        "label": "Cooling Degree Days (Spinoni)",
        "applicable_to": {"temperature"},
        "expects_statistic": "daily_mean",
        "params": ("cdd_base",),
        "compute": _temp_cdd_spinoni,
        "units_out": lambda u: (u or "degC") + "*day",
        "climatology": "mean",
    },
    "gdd": {
        "category": "temperature",
        "label": "Growing Degree Days",
        "applicable_to": {"temperature"},
        "expects_statistic": "daily_mean",
        "params": ("gdd_base",),
        "compute": _temp_gdd,
        "units_out": lambda u: (u or "degC") + "*day",
        "climatology": "mean",
    },
    "frost_days": {
        "category": "temperature",
        "label": "Frost days (Tmin < 0 degC)",
        "applicable_to": {"temperature"},
        "expects_statistic": "daily_minimum",
        "params": (),
        "compute": _temp_frost_days,
        "units_out": "days", "climatology": "mean",
    },
    "ice_days": {
        "category": "temperature",
        "label": "Ice days (Tmax < 0 degC)",
        "applicable_to": {"temperature"},
        "expects_statistic": "daily_maximum",
        "params": (),
        "compute": _temp_ice_days,
        "units_out": "days", "climatology": "mean",
    },
    "summer_days": {
        "category": "temperature",
        "label": "Summer days (Tmax > 25 degC)",
        "applicable_to": {"temperature"},
        "expects_statistic": "daily_maximum",
        "params": (),
        "compute": _temp_summer_days,
        "units_out": "days", "climatology": "mean",
    },
    "tropical_nights": {
        "category": "temperature",
        "label": "Tropical nights (Tmin > 20 degC)",
        "applicable_to": {"temperature"},
        "expects_statistic": "daily_minimum",
        "params": (),
        "compute": _temp_tropical_nights,
        "units_out": "days", "climatology": "mean",
    },
    "hot_days": {
        "category": "temperature",
        "label": "Hot days (Tmax > 30 degC)",
        "applicable_to": {"temperature"},
        "expects_statistic": "daily_maximum",
        "params": (),
        "compute": _temp_hot_days,
        "units_out": "days", "climatology": "mean",
    },
    "heating_season_length": {
        "category": "temperature",
        "label": "Heating season length (G.POT, Tmean < threshold)",
        "applicable_to": {"temperature"},
        "expects_statistic": "daily_mean",
        "params": ("heating_threshold",),
        "compute": _temp_heating_season_length,
        "units_out": "days", "climatology": "mean",
    },
    "cooling_season_length": {
        "category": "temperature",
        "label": "Cooling season length (G.POT, Tmean > threshold)",
        "applicable_to": {"temperature"},
        "expects_statistic": "daily_mean",
        "params": ("cooling_threshold",),
        "compute": _temp_cooling_season_length,
        "units_out": "days", "climatology": "mean",
    },
    "gsl": {
        "category": "temperature",
        "label": "Growing Season Length (WMO, NH)",
        "applicable_to": {"temperature"},
        "expects_statistic": "daily_mean",
        "params": (),
        "compute": _temp_gsl,
        "units_out": "days", "climatology": "mean",
    },
    "txn": {
        "category": "temperature",
        "label": "Annual min of daily-max (TXn)",
        "applicable_to": {"temperature"},
        "expects_statistic": "daily_maximum",
        "params": (),
        "compute": _temp_txn,
        "units_out": _passthrough_units, "climatology": "mean",
    },
    "tnx": {
        "category": "temperature",
        "label": "Annual max of daily-min (TNx)",
        "applicable_to": {"temperature"},
        "expects_statistic": "daily_minimum",
        "params": (),
        "compute": _temp_tnx,
        "units_out": _passthrough_units, "climatology": "mean",
    },

    # ---- Precipitation-specific (working units: mm) ----
    "total_precip": {
        "category": "precipitation",
        "label": "Total precipitation (annual)",
        "applicable_to": {"precipitation"},
        "expects_statistic": None,
        "params": (),
        "compute": _precip_total,
        "units_out": _passthrough_units, "climatology": "mean",
    },
    "wet_days": {
        "category": "precipitation",
        "label": "Wet days (precip >= threshold mm)",
        "applicable_to": {"precipitation"},
        "expects_statistic": None,
        "params": ("precip_wet_threshold",),
        "compute": _precip_wet_days,
        "units_out": "days", "climatology": "mean",
    },
    "consecutive_dry_days": {
        "category": "precipitation",
        "label": "Consecutive dry days (CDD ETCCDI index)",
        "applicable_to": {"precipitation"},
        "expects_statistic": None,
        "params": ("precip_wet_threshold",),
        "compute": _precip_consecutive_dry_days,
        "units_out": "days", "climatology": "mean",
    },
    "consecutive_wet_days": {
        "category": "precipitation",
        "label": "Consecutive wet days (CWD ETCCDI index)",
        "applicable_to": {"precipitation"},
        "expects_statistic": None,
        "params": ("precip_wet_threshold",),
        "compute": _precip_consecutive_wet_days,
        "units_out": "days", "climatology": "mean",
    },
    "max_1day_precip": {
        "category": "precipitation",
        "label": "Maximum 1-day precipitation (Rx1day)",
        "applicable_to": {"precipitation"},
        "expects_statistic": None,
        "params": (),
        "compute": _precip_max_1day,
        "units_out": _passthrough_units, "climatology": "mean",
    },

    # ---- DTR (diurnal temperature range, multi-input) ----
    # multi_input="dtr" tells the engine to build a synthetic Tmax-Tmin stack
    # from the variable's daily_maximum AND daily_minimum directories before
    # invoking compute(). The metric does NOT see the regular (var, stat)
    # statistic the user selected for the rest of the run.
    "dtr_mean": {
        "category": "dtr",
        "label": "DTR mean (annual mean of Tmax-Tmin)",
        "applicable_to": {"temperature"},
        "expects_statistic": None,
        "multi_input": "dtr",
        "params": (),
        "compute": lambda s: np.nanmean(s, axis=0),
        "units_out": "degC", "climatology": "mean",
    },
    "dtr_max": {
        "category": "dtr",
        "label": "DTR maximum (annual max of Tmax-Tmin)",
        "applicable_to": {"temperature"},
        "expects_statistic": None,
        "multi_input": "dtr",
        "params": (),
        "compute": lambda s: np.nanmax(s, axis=0),
        "units_out": "degC", "climatology": "mean",
    },
    "dtr_p90": {
        "category": "dtr",
        "label": "DTR 90th percentile (annual)",
        "applicable_to": {"temperature"},
        "expects_statistic": None,
        "multi_input": "dtr",
        "params": (),
        "compute": lambda s: np.nanpercentile(s, 90, axis=0),
        "units_out": "degC", "climatology": "mean",
    },
    "dtr_std": {
        "category": "dtr",
        "label": "DTR standard deviation (annual)",
        "applicable_to": {"temperature"},
        "expects_statistic": None,
        "multi_input": "dtr",
        "params": (),
        "compute": lambda s: np.nanstd(s, axis=0),
        "units_out": "degC", "climatology": "mean",
    },

    # ---- Radiation-specific (working units: raw J/m^2/day) ----
    "radiation_w_per_m2_mean": {
        "category": "radiation",
        "label": "Mean radiation power (W/m^2)",
        "applicable_to": {"radiation"},
        "expects_statistic": "daily_sum",
        "params": (),
        "compute": _rad_w_per_m2_mean,
        "units_out": "W/m2", "climatology": "mean",
    },
    "radiation_kwh_per_m2_year": {
        "category": "radiation",
        "label": "Annual radiation energy (kWh/m^2)",
        "applicable_to": {"radiation"},
        "expects_statistic": "daily_sum",
        "params": (),
        "compute": _rad_kwh_per_m2_year,
        "units_out": "kWh/m2", "climatology": "mean",
    },
    "photovoltaic_potential": {
        "category": "radiation",
        "label": "Photovoltaic potential (kWh/m^2/year * PR)",
        "applicable_to": {"radiation"},
        "expects_statistic": "daily_sum",
        "params": ("pv_performance_ratio",),
        "compute": _rad_pv_potential,
        "units_out": "kWh/m2", "climatology": "mean",
    },
    "radiation_above_threshold_days": {
        "category": "radiation",
        "label": "Days with daily mean radiation > threshold (W/m^2)",
        "applicable_to": {"radiation"},
        "expects_statistic": "daily_sum",
        "params": ("radiation_threshold_w_per_m2",),
        "compute": _rad_above_threshold_days,
        "units_out": "days", "climatology": "mean",
    },

    # ---- Wind-specific (multi-input: requires u10 AND v10) ----
    # multi_input="wind" tells the engine to build sqrt(u^2+v^2) magnitude
    # stack from the u10 and v10 directories before invoking compute().
    "wind_speed_mean": {
        "category": "wind",
        "label": "Wind speed mean (annual)",
        "applicable_to": {"wind"},
        "expects_statistic": None,
        "multi_input": "wind",
        "params": (),
        "compute": lambda s: np.nanmean(s, axis=0),
        "units_out": "m/s", "climatology": "mean",
    },
    "wind_speed_max": {
        "category": "wind",
        "label": "Wind speed max (annual)",
        "applicable_to": {"wind"},
        "expects_statistic": None,
        "multi_input": "wind",
        "params": (),
        "compute": lambda s: np.nanmax(s, axis=0),
        "units_out": "m/s", "climatology": "mean",
    },
    "wind_speed_p90": {
        "category": "wind",
        "label": "Wind speed 90th percentile (annual)",
        "applicable_to": {"wind"},
        "expects_statistic": None,
        "multi_input": "wind",
        "params": (),
        "compute": lambda s: np.nanpercentile(s, 90, axis=0),
        "units_out": "m/s", "climatology": "mean",
    },
    "wind_power_density": {
        "category": "wind",
        "label": "Wind power density (0.5*rho*<v^3>)",
        "applicable_to": {"wind"},
        "expects_statistic": None,
        "multi_input": "wind",
        "params": ("wind_air_density",),
        "compute": _wind_power_density,
        "units_out": "W/m2", "climatology": "mean",
    },
    "calm_days": {
        "category": "wind",
        "label": "Calm days (wind < threshold m/s)",
        "applicable_to": {"wind"},
        "expects_statistic": None,
        "multi_input": "wind",
        "params": ("wind_calm_threshold",),
        "compute": _wind_calm_days,
        "units_out": "days", "climatology": "mean",
    },
}

CLIMATOLOGY_AGGREGATORS = {
    "mean":   lambda s: np.nanmean(s, axis=0),
    "median": lambda s: np.nanmedian(s, axis=0),
    "min":    lambda s: np.nanmin(s, axis=0),
    "max":    lambda s: np.nanmax(s, axis=0),
}


def apply_kelvin_conversion(stack, variable_short, log):
    """Convert stack from Kelvin to Celsius in-place if variable is in
    KELVIN_VARIABLES. Returns "degC" if converted, None otherwise."""
    if str(variable_short).lower() in KELVIN_VARIABLES:
        stack -= 273.15
        log("    K to degC conversion applied to stack")
        return "degC"
    return None


def apply_precip_conversion(stack, variable_short, log):
    """Convert precipitation stack from m to mm in-place. ERA5-Land tp is
    delivered in metres of water equivalent per day; the toolbox standardises
    to mm so user-facing thresholds (e.g. wet day >= 1 mm) are intuitive.
    Returns "mm" if converted, None otherwise."""
    if str(variable_short).lower() in PRECIPITATION_VARIABLES:
        stack *= 1000.0
        log("    m to mm conversion applied to stack")
        return "mm"
    return None


def resolve_working_units(variable_short, input_units):
    """Mirror the unit conversion logic of apply_*_conversion without
    requiring the stack. Used by climatology / IPMA-normal writers to
    resolve units_out for metrics that pass through the working units."""
    sn = str(variable_short).lower()
    if sn in KELVIN_VARIABLES:
        return "degC"
    if sn in PRECIPITATION_VARIABLES:
        return "mm"
    return input_units or "unknown"


def _to_disk_array(arr, mask):
    """Replace NaN with NODATA_VAL and apply boolean mask if given.

    Single allocation: builds the invalid-mask once (NaN OR outside polygon)
    and writes NODATA_VAL where invalid. Saves one full-array copy per
    metric vs. the naive two-pass np.where chain.
    """
    out = arr.astype(np.float32, copy=True)
    if mask is None:
        invalid = np.isnan(out)
    else:
        invalid = np.isnan(out) | (~mask)
    out[invalid] = NODATA_VAL
    return out


# ===========================================================================
# Climate indicators (Tool 3) - Annual + climatology processing
# ===========================================================================

def process_year_metrics(variable, statistic, year, daily_paths,
                          metrics_to_run, out_root, mask, gt, srs_wkt,
                          input_units, params, log):
    """Process one (variable, statistic, year): load daily stack, compute
    each selected metric, write annual TIFFs.

    Returns (annual_nan_arrays, records):
      annual_nan_arrays: dict {metric_name: 2D NaN-aware array} used by
                        caller to build the per-metric climatology buffer.
      records:          list of manifest dicts for the TIFFs just written.
    """
    log("  loading {} daily TIFFs into stack...".format(len(daily_paths)))
    stack = read_paths_into_stack(daily_paths)
    log("    stack shape={}, dtype={}".format(stack.shape, stack.dtype))

    # Auto-unit conversion: K -> degC for temperature, m -> mm for precip.
    # At most one applies (variables belong to at most one class).
    working_units = (
        apply_kelvin_conversion(stack, variable, log)
        or apply_precip_conversion(stack, variable, log)
        or (input_units or "unknown")
    )

    annual_nan = {}
    records = []
    for metric_name in metrics_to_run:
        spec = METRIC_REGISTRY.get(metric_name)
        if spec is None:
            log("    unknown metric '{}', skipping".format(metric_name), "warn")
            continue
        try:
            kw = {k: params[k] for k in spec.get("params", ()) if k in params}
            arr = spec["compute"](stack, **kw)
        except Exception as e:
            log("    ERROR computing '{}' on year {}: {}".format(
                metric_name, year, e
            ), "warn")
            continue

        units_out_spec = spec.get("units_out", _passthrough_units)
        units_out = (units_out_spec(working_units) if callable(units_out_spec)
                     else str(units_out_spec))

        annual_dir = os.path.join(
            out_root, "annual", variable, statistic, metric_name
        )
        os.makedirs(annual_dir, exist_ok=True)
        out_name = "{}_{}_{}_{}.tif".format(
            variable, statistic, metric_name, year
        )
        out_path = os.path.join(annual_dir, out_name)
        write_geotiff_gdal(_to_disk_array(arr, mask), out_path, gt, srs_wkt)

        annual_nan[metric_name] = arr
        records.append({
            "raster": out_path,
            "variable": variable,
            "statistic": statistic,
            "metric": metric_name,
            "category": spec.get("category", "generic"),
            "period_kind": "annual",
            "year": year,
            "start_year": year,
            "end_year": year,
            "units": units_out,
        })

    del stack  # release the ~hundreds-of-MB stack before next year
    return annual_nan, records


def process_year_dtr(variable, year, tmax_paths, tmin_paths,
                      dtr_metrics, out_root, mask, gt, srs_wkt, log):
    """DTR (diurnal temperature range) annual pipeline.

    Builds a synthetic daily stack of (Tmax - Tmin) in degC, then applies
    each selected DTR metric. Tmax/Tmin paths must be aligned 1:1 by day
    (caller is responsible for ensuring this; we validate length only).
    Outputs go under <out>/annual/<variable>/dtr/<metric>/.
    """
    if len(tmax_paths) != len(tmin_paths):
        raise ValueError(
            "DTR year {}: Tmax has {} files, Tmin has {} files".format(
                year, len(tmax_paths), len(tmin_paths)
            )
        )
    log("  DTR: loading {} Tmax + {} Tmin daily TIFFs...".format(
        len(tmax_paths), len(tmin_paths)
    ))
    tmax_stack = read_paths_into_stack(tmax_paths)
    tmin_stack = read_paths_into_stack(tmin_paths)

    # Convert once; suppress duplicate log line on second call
    apply_kelvin_conversion(tmax_stack, variable, log)
    apply_kelvin_conversion(tmin_stack, variable, _silent_log)

    # In-place subtraction reuses tmax_stack buffer; saves one full-stack
    # allocation (~38 MB peak for Portugal). dtr_stack and tmax_stack are
    # the same memory; alias for readability.
    np.subtract(tmax_stack, tmin_stack, out=tmax_stack)
    dtr_stack = tmax_stack
    del tmin_stack

    annual_nan = {}
    records = []
    for metric_name in dtr_metrics:
        spec = METRIC_REGISTRY.get(metric_name)
        if spec is None or spec.get("multi_input") != "dtr":
            log("    DTR pipeline rejected metric '{}' (not a DTR spec)".format(
                metric_name
            ), "warn")
            continue
        try:
            arr = spec["compute"](dtr_stack)
        except Exception as e:
            log("    ERROR computing DTR '{}' for year {}: {}".format(
                metric_name, year, e
            ), "warn")
            continue

        units_out_spec = spec.get("units_out", "degC")
        units_out = (units_out_spec("degC") if callable(units_out_spec)
                     else str(units_out_spec))

        annual_dir = os.path.join(
            out_root, "annual", variable, "dtr", metric_name
        )
        os.makedirs(annual_dir, exist_ok=True)
        out_name = "{}_dtr_{}_{}.tif".format(variable, metric_name, year)
        out_path = os.path.join(annual_dir, out_name)
        write_geotiff_gdal(_to_disk_array(arr, mask), out_path, gt, srs_wkt)

        annual_nan[metric_name] = arr
        records.append({
            "raster": out_path,
            "variable": variable,
            "statistic": "dtr",
            "metric": metric_name,
            "category": "dtr",
            "period_kind": "annual",
            "year": year,
            "start_year": year,
            "end_year": year,
            "units": units_out,
        })

    del dtr_stack
    return annual_nan, records


def process_year_wind(year, u_paths, v_paths, wind_metrics,
                       out_root, mask, gt, srs_wkt, params, log):
    """Wind annual pipeline.

    Builds a synthetic daily stack of wind speed magnitude = sqrt(u^2 + v^2)
    in m/s, then applies each selected wind metric. Output path uses a
    synthetic variable label "wind" so wind outputs sit next to the per-
    variable trees in the output root.
    """
    if len(u_paths) != len(v_paths):
        raise ValueError(
            "WIND year {}: u10 has {} files, v10 has {} files".format(
                year, len(u_paths), len(v_paths)
            )
        )
    log("  WIND: loading {} u10 + {} v10 daily TIFFs...".format(
        len(u_paths), len(v_paths)
    ))
    u_stack = read_paths_into_stack(u_paths)
    v_stack = read_paths_into_stack(v_paths)
    # In-place magnitude: u <- sqrt(u^2 + v^2). Reuses u_stack buffer; saves
    # one full-stack allocation (~38 MB peak for Portugal).
    np.square(u_stack, out=u_stack)
    np.square(v_stack, out=v_stack)
    np.add(u_stack, v_stack, out=u_stack)
    np.sqrt(u_stack, out=u_stack)
    speed_stack = u_stack
    del v_stack

    annual_nan = {}
    records = []
    for metric_name in wind_metrics:
        spec = METRIC_REGISTRY.get(metric_name)
        if spec is None or spec.get("multi_input") != "wind":
            log("    WIND pipeline rejected metric '{}' (not a wind spec)".format(
                metric_name
            ), "warn")
            continue
        try:
            kw = {k: params[k] for k in spec.get("params", ()) if k in params}
            arr = spec["compute"](speed_stack, **kw)
        except Exception as e:
            log("    ERROR computing wind '{}' for year {}: {}".format(
                metric_name, year, e
            ), "warn")
            continue

        units_out_spec = spec.get("units_out", "m/s")
        units_out = (units_out_spec("m/s") if callable(units_out_spec)
                     else str(units_out_spec))

        annual_dir = os.path.join(
            out_root, "annual", "wind", metric_name
        )
        os.makedirs(annual_dir, exist_ok=True)
        out_name = "wind_{}_{}.tif".format(metric_name, year)
        out_path = os.path.join(annual_dir, out_name)
        write_geotiff_gdal(_to_disk_array(arr, mask), out_path, gt, srs_wkt)

        annual_nan[metric_name] = arr
        records.append({
            "raster": out_path,
            "variable": "wind",
            "statistic": "derived",
            "metric": metric_name,
            "category": "wind",
            "period_kind": "annual",
            "year": year,
            "start_year": year,
            "end_year": year,
            "units": units_out,
        })

    del speed_stack
    return annual_nan, records


def process_climatology(variable, statistic, metric_name, buffer,
                         out_root, mask, gt, srs_wkt, units_out,
                         period_kind, period_label, start_year, end_year):
    """Aggregate the annual buffer into a climatology / IPMA normal TIFF.

    period_kind: "climatology" -> writes to climatology/ folder
                 "ipma_normal" -> writes to ipma_normals/ folder
    period_label: e.g. "1981-2010"; if empty, derived from start/end.
    """
    if not buffer:
        return None
    spec = METRIC_REGISTRY.get(metric_name, {})
    aggregator = CLIMATOLOGY_AGGREGATORS[spec.get("climatology", "mean")]
    stack = np.stack([a for _, a in buffer], axis=0)
    clim = aggregator(stack)

    folder_name = "climatology" if period_kind == "climatology" else "ipma_normals"
    out_dir = os.path.join(out_root, folder_name, variable, statistic)
    os.makedirs(out_dir, exist_ok=True)
    suffix = period_label if period_label else "{}-{}".format(start_year, end_year)
    kind_token = "climatology" if period_kind == "climatology" else "normal"
    out_name = "{}_{}_{}_{}_{}.tif".format(
        variable, statistic, metric_name, kind_token, suffix
    )
    out_path = os.path.join(out_dir, out_name)
    write_geotiff_gdal(_to_disk_array(clim, mask), out_path, gt, srs_wkt)

    return {
        "raster": out_path,
        "variable": variable,
        "statistic": statistic,
        "metric": metric_name,
        "category": spec.get("category", "generic"),
        "period_kind": period_kind,
        "year": "",
        "start_year": start_year,
        "end_year": end_year,
        "units": units_out,
    }


# ===========================================================================
# Core conversion
# ===========================================================================

def convert_one_nc(nc_path, out_root, user_variable=None, clip_fc=None,
                   messages=None, verbose=False, advance_progressor=None):
    """Convert one ERA5-Land NetCDF into per-day GeoTIFFs.

    Output structure: out_root/<variable>/<year>/<variable>_YYYY-MM-DD.tif

    Returns list of record dicts for manifest.
    """
    def log(msg, level="msg"):
        if messages is None:
            arcpy.AddMessage(msg)
            return
        if level == "warn":
            messages.addWarningMessage(msg)
        elif level == "err":
            messages.addErrorMessage(msg)
        else:
            messages.addMessage(msg)

    t_start = _time.time()

    # Detect daily statistic from filename (era5land-downloads pattern).
    # When not detectable, use a literal marker so the output layout still
    # has a statistic level (forces awareness downstream rather than silent
    # collision when multiple statistics share a variable name).
    statistic = detect_statistic_from_filename(nc_path) or UNKNOWN_STATISTIC
    if statistic == UNKNOWN_STATISTIC:
        log("  statistic: could not infer from filename, using '{}'".format(
            UNKNOWN_STATISTIC
        ), "warn")

    if verbose:
        log("  opening NetCDF...")
    ds = netCDF4.Dataset(nc_path)
    try:
        time_dim = find_dim(ds, TIME_DIM_CANDIDATES)
        lat_dim = find_dim(ds, LAT_DIM_CANDIDATES)
        lon_dim = find_dim(ds, LON_DIM_CANDIDATES)

        missing = [n for n, v in (("time", time_dim), ("lat", lat_dim),
                                  ("lon", lon_dim)) if not v]
        if missing:
            raise ValueError(
                "Missing dimensions {}. Available: {}".format(
                    missing, list(ds.dimensions.keys())
                )
            )

        var_name = detect_data_variable(ds, user_variable)
        var_meta = get_var_metadata(ds, var_name)

        # Read coordinates
        lat_var = find_coord_var(ds, LAT_DIM_CANDIDATES)
        lon_var = find_coord_var(ds, LON_DIM_CANDIDATES)
        if lat_var is None or lon_var is None:
            raise ValueError("Could not find lat/lon coordinate variables.")

        lat = ds.variables[lat_var][:]
        lon = ds.variables[lon_var][:]
        if isinstance(lat, np.ma.MaskedArray):
            lat = lat.filled(np.nan)
        if isinstance(lon, np.ma.MaskedArray):
            lon = lon.filled(np.nan)

        lat_ascending = bool(lat[1] > lat[0]) if len(lat) > 1 else True

        # Time
        dates = get_years_from_time(ds, time_dim)
        n_t = len(dates)

        # Consolidated 1-line header (or verbose multi-line)
        if verbose:
            log("  variable='{}', shape={}, units='{}', long_name='{}'".format(
                var_name, var_meta["shape"], var_meta["units"],
                var_meta["long_name"]
            ))
            log("  time steps={}, range {}-{:02d}-{:02d} to "
                "{}-{:02d}-{:02d}".format(
                    n_t,
                    dates[0][0], dates[0][1], dates[0][2],
                    dates[-1][0], dates[-1][1], dates[-1][2]
                ))
        else:
            log("  {} {} ({}), {} steps {}-{}".format(
                var_name, statistic, var_meta["units"],
                n_t, dates[0][0], dates[-1][0]
            ))

        # Spatial filtering
        if clip_fc:
            clip_extent = arcpy.Describe(clip_fc).extent
            lat_slice, lon_slice = compute_bbox_slices(lat, lon, clip_extent)
            sliced_lat = lat[lat_slice]
            sliced_lon = lon[lon_slice]
            if verbose:
                log("  bbox slice: {} rows x {} cols (from {} x {})".format(
                    len(sliced_lat), len(sliced_lon), len(lat), len(lon)
                ))
            mask = get_or_compute_mask(clip_fc, sliced_lat, sliced_lon, lat_ascending)
        else:
            lat_slice = slice(None)
            lon_slice = slice(None)
            sliced_lat = lat
            sliced_lon = lon
            mask = None

        var = ds.variables[var_name]
        fill_value = var_meta["fill_value"]

        gt = compute_geotransform(sliced_lon, sliced_lat, lat_ascending)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        srs_wkt = srs.ExportToWkt()

        # Variable + statistic subfolder. Statistic level is always present
        # in the output path to prevent silent collision when the user
        # processes multiple statistics of the same variable in one run.
        stat_dir = os.path.join(out_root, var_name, statistic)
        os.makedirs(stat_dir, exist_ok=True)

        records = []
        for i in range(n_t):
            year, month, day = dates[i]
            year_dir = os.path.join(stat_dir, str(year))
            if not os.path.isdir(year_dir):
                os.makedirs(year_dir, exist_ok=True)

            out_name = "{}_{}_{:04d}-{:02d}-{:02d}.tif".format(
                var_name, statistic, year, month, day
            )
            out_path = os.path.join(year_dir, out_name)

            try:
                arr = var[i, lat_slice, lon_slice]
                if isinstance(arr, np.ma.MaskedArray):
                    arr = arr.filled(NODATA_VAL)
                arr = np.asarray(arr, dtype=np.float32)

                if fill_value is not None:
                    try:
                        fv = float(fill_value)
                        arr = np.where(arr == np.float32(fv), NODATA_VAL, arr)
                    except (TypeError, ValueError):
                        pass

                if mask is not None:
                    arr = np.where(mask, arr, NODATA_VAL)

                if lat_ascending:
                    arr_out = arr[::-1, :]
                else:
                    arr_out = arr

                write_geotiff_gdal(arr_out, out_path, gt, srs_wkt)

                records.append({
                    "raster": out_path,
                    "variable": var_name,
                    "statistic": statistic,
                    "standard_name": var_meta["standard_name"],
                    "long_name": var_meta["long_name"],
                    "units": var_meta["units"],
                    "year": year,
                    "month": month,
                    "day": day,
                    "date": "{:04d}-{:02d}-{:02d}".format(year, month, day),
                    "source_nc": os.path.basename(nc_path),
                })

                if verbose and (i + 1) % 30 == 0:
                    log("    {}/{} timesteps written".format(i + 1, n_t))

            except Exception as e:
                log("  ERROR at index {} (date {:04d}-{:02d}-{:02d}): {}".format(
                    i, year, month, day, e
                ), "warn")
            finally:
                if advance_progressor:
                    advance_progressor()

        elapsed = _time.time() - t_start
        per_slice = elapsed / len(records) if records else 0
        log("  wrote {} GeoTIFFs in {:.1f}s ({:.3f}s/slice)".format(
            len(records), elapsed, per_slice
        ))
        return records, {
            "variable": var_name,
            "statistic": statistic,
            "grid": grid_signature(lat, lon),
            "lat_size": len(lat),
            "lon_size": len(lon),
            "lat_min": float(lat.min()),
            "lat_max": float(lat.max()),
            "lon_min": float(lon.min()),
            "lon_max": float(lon.max()),
        }

    finally:
        ds.close()


def count_time_slices(nc_path):
    """Quick count of time steps in a NetCDF."""
    try:
        ds = netCDF4.Dataset(nc_path)
        try:
            time_dim = find_dim(ds, TIME_DIM_CANDIDATES)
            if time_dim:
                return ds.dimensions[time_dim].size
        finally:
            ds.close()
    except Exception:
        pass
    return 0


# ===========================================================================
# Layer resolution (for GPRasterLayer dropdown input)
# ===========================================================================

def get_layer_source_nc(layer_name):
    """Resolve a map layer name to its underlying .nc file path."""
    try:
        aprx = arcpy.mp.ArcGISProject("CURRENT")
    except Exception:
        return None

    candidate = None
    for m in aprx.listMaps():
        for lyr in m.listLayers():
            if lyr.name == layer_name:
                candidate = lyr
                break
            if getattr(lyr, "longName", None) == layer_name:
                candidate = lyr
                break
        if candidate:
            break

    if candidate is None:
        return None

    try:
        cp = candidate.connectionProperties
        if cp:
            ci = cp.get("connection_info", {}) or {}
            ds = cp.get("dataset") or ci.get("dataset")
            db = (cp.get("database") or ci.get("database")
                  or ci.get("server"))
            if ds:
                if os.path.isabs(str(ds)) and str(ds).lower().endswith(".nc"):
                    if os.path.isfile(ds):
                        return ds
                if str(ds).lower().endswith(".nc") and db and os.path.isdir(db):
                    full = os.path.join(db, ds)
                    if os.path.isfile(full):
                        return full
    except Exception:
        pass

    try:
        src = candidate.dataSource
        if src:
            if src.lower().endswith(".nc") and os.path.isfile(src):
                return src
            parts = re.split(r"[\\/]", src)
            for j in range(len(parts), 0, -1):
                cand = os.sep.join(parts[:j])
                if cand.lower().endswith(".nc") and os.path.isfile(cand):
                    return cand
    except Exception:
        pass

    return None


def resolve_nc_inputs(file_param_raw, layer_param_raw, messages):
    """Combine file paths and layer references."""
    nc_files = []

    if file_param_raw:
        for p in file_param_raw.split(";"):
            p = p.strip().strip("'\"")
            if p:
                nc_files.append(p)

    if layer_param_raw:
        for lyr in layer_param_raw.split(";"):
            lyr = lyr.strip().strip("'\"")
            if not lyr:
                continue
            path = get_layer_source_nc(lyr)
            if path:
                nc_files.append(path)
            else:
                messages.addWarningMessage(
                    "Could not resolve layer '{}' to a .nc file.".format(lyr)
                )

    seen = set()
    deduped = []
    for f in nc_files:
        key = os.path.normcase(os.path.abspath(f))
        if key not in seen:
            seen.add(key)
            deduped.append(f)
    return deduped


# ===========================================================================
# Toolbox
# ===========================================================================

class Toolbox(object):
    def __init__(self):
        self.label = "ERA5-Land ArcGIS Tools"
        self.alias = "era5land_tools"
        self.tools = [
            ConvertERA5LandToGeoTIFF,
            InspectERA5LandNetCDF,
            ComputeERA5LandIndicators,
            ExtractRegionalTimeSeries,
        ]


# ---------------------------------------------------------------------------
# Tool 1: Convert
# ---------------------------------------------------------------------------

class ConvertERA5LandToGeoTIFF(object):
    def __init__(self):
        self.label = "1. Convert ERA5-Land NetCDF to GeoTIFF"
        self.description = (
            "Convert ERA5-Land daily statistics NetCDFs to per-day GeoTIFFs. "
            "Variable-agnostic: works with any CF-compliant variable having "
            "(time, lat, lon) dimensions. Output organised as "
            "<output>/<variable>/<statistic>/<year>/"
            "<variable>_<statistic>_YYYY-MM-DD.tif. "
            "Statistic (mean, maximum, minimum, sum, stddev) is auto-detected "
            "from source filename pattern era5land_<variable>_<stat>_<year>.nc; "
            "falls back to 'unknown_stat' if not detectable. Multi-NetCDF "
            "inputs are validated for grid consistency."
        )
        self.canRunInBackground = False

    def getParameterInfo(self):
        p_files = arcpy.Parameter(
            displayName="Input NetCDF files (browse to disk)",
            name="in_nc_files", datatype="DEFile",
            parameterType="Optional", direction="Input", multiValue=True,
        )
        p_files.filter.list = ["nc"]

        p_layers = arcpy.Parameter(
            displayName="OR input multidimensional raster layers (from map)",
            name="in_nc_layers", datatype="GPRasterLayer",
            parameterType="Optional", direction="Input", multiValue=True,
        )

        p_out = arcpy.Parameter(
            displayName="Output folder",
            name="out_folder", datatype="DEFolder",
            parameterType="Required", direction="Input",
        )

        p_var = arcpy.Parameter(
            displayName="Variable name (blank = auto-detect)",
            name="variable", datatype="GPString",
            parameterType="Optional", direction="Input",
        )

        p_clip = arcpy.Parameter(
            displayName="Clip polygon feature class (optional)",
            name="clip_fc", datatype="GPFeatureLayer",
            parameterType="Optional", direction="Input",
        )
        p_clip.filter.list = ["Polygon"]

        p_manifest = arcpy.Parameter(
            displayName="Write manifest CSV",
            name="write_manifest", datatype="GPBoolean",
            parameterType="Optional", direction="Input",
        )
        p_manifest.value = True

        p_strict_grid = arcpy.Parameter(
            displayName="Strict grid validation (error if inputs disagree)",
            name="strict_grid", datatype="GPBoolean",
            parameterType="Optional", direction="Input",
        )
        p_strict_grid.value = True

        p_verbose = arcpy.Parameter(
            displayName="Verbose logging",
            name="verbose", datatype="GPBoolean",
            parameterType="Optional", direction="Input",
        )
        p_verbose.value = False

        return [p_files, p_layers, p_out, p_var, p_clip, p_manifest,
                p_strict_grid, p_verbose]

    def isLicensed(self):
        return True

    def updateMessages(self, parameters):
        p_files = parameters[0]
        p_layers = parameters[1]
        if not p_files.value and not p_layers.value:
            p_files.setErrorMessage(
                "Provide at least one input source."
            )

    def execute(self, parameters, messages):
        file_raw = parameters[0].valueAsText
        layer_raw = parameters[1].valueAsText
        out_folder = parameters[2].valueAsText
        user_variable = parameters[3].valueAsText or None
        clip_fc = parameters[4].valueAsText if parameters[4].value else None
        write_manifest = parameters[5].value
        strict_grid = parameters[6].value
        verbose = bool(parameters[7].value)

        nc_files = resolve_nc_inputs(file_raw, layer_raw, messages)
        if not nc_files:
            messages.addErrorMessage("No valid NetCDF inputs.")
            return

        if not os.path.isdir(out_folder):
            os.makedirs(out_folder)

        slice_counts = [count_time_slices(nc) for nc in nc_files]
        total_slices = sum(slice_counts)

        messages.addMessage("Processing {} NetCDF file(s), {} total time slices".format(
            len(nc_files), total_slices
        ))
        messages.addMessage("Engine: netCDF4 + GDAL (fast path)")
        messages.addMessage("Output: {}".format(out_folder))
        if user_variable:
            messages.addMessage("Variable override: {}".format(user_variable))
        if clip_fc:
            messages.addMessage("Clip: {}".format(clip_fc))
        if strict_grid:
            messages.addMessage("Grid validation: STRICT")
        if verbose:
            messages.addMessage("Verbose logging enabled")

        _MASK_CACHE.clear()

        progressor_state = {"pos": 0}

        def advance():
            progressor_state["pos"] += 1
            arcpy.SetProgressorPosition(progressor_state["pos"])

        arcpy.SetProgressor(
            type="step", message="Processing NetCDF time slices...",
            min_range=0, max_range=max(total_slices, 1), step_value=1,
        )

        all_records = []
        grid_signatures_seen = {}  # var_name -> grid_signature
        t_total_start = _time.time()
        try:
            for idx, nc in enumerate(nc_files):
                fname = os.path.basename(nc)
                arcpy.SetProgressorLabel(
                    "File {}/{}: {}".format(idx + 1, len(nc_files), fname)
                )
                messages.addMessage("--- [{}/{}] {}".format(
                    idx + 1, len(nc_files), fname
                ))
                try:
                    recs, file_info = convert_one_nc(
                        nc, out_folder,
                        user_variable=user_variable,
                        clip_fc=clip_fc, messages=messages,
                        verbose=verbose, advance_progressor=advance,
                    )

                    # Grid consistency check (per variable)
                    var = file_info["variable"]
                    grid = file_info["grid"]
                    if var in grid_signatures_seen:
                        prev_grid = grid_signatures_seen[var]
                        if not grids_match(prev_grid, grid):
                            msg = (
                                "Grid mismatch for variable '{}' in {}: "
                                "previous {} vs this {}.".format(
                                    var, fname, prev_grid, grid
                                )
                            )
                            if strict_grid:
                                messages.addErrorMessage(msg + " Aborting (strict mode).")
                                raise RuntimeError(msg)
                            else:
                                messages.addWarningMessage(msg)
                    else:
                        grid_signatures_seen[var] = grid

                    all_records.extend(recs)
                except Exception as e:
                    messages.addErrorMessage("FAILED on {}: {}".format(fname, e))
                    if strict_grid and "Grid mismatch" in str(e):
                        break
        finally:
            arcpy.ResetProgressor()

        if write_manifest and all_records:
            manifest_path = os.path.join(out_folder, "manifest.csv")
            fields = ["raster", "variable", "statistic", "standard_name",
                      "long_name", "units", "year", "month", "day", "date",
                      "source_nc"]
            with open(manifest_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=fields)
                writer.writeheader()
                writer.writerows(all_records)
            messages.addMessage("Manifest: {}".format(manifest_path))

        total_elapsed = _time.time() - t_total_start
        avg = total_elapsed / len(all_records) if all_records else 0
        messages.addMessage(
            "Done. {} GeoTIFFs from {} NetCDF(s) in {:.1f}s ({:.3f}s/slice avg).".format(
                len(all_records), len(nc_files), total_elapsed, avg
            )
        )

        if grid_signatures_seen:
            messages.addMessage(
                "Variables processed: {}".format(", ".join(grid_signatures_seen.keys()))
            )


# ---------------------------------------------------------------------------
# Tool 2: Inspect
# ---------------------------------------------------------------------------

class InspectERA5LandNetCDF(object):
    def __init__(self):
        self.label = "2. Inspect ERA5-Land NetCDF (debug)"
        self.description = (
            "Print dimensions, variables, attributes and time-axis decoding "
            "for ERA5-Land NetCDFs. Use before Tool 1 to verify file "
            "structure and detect anomalies."
        )
        self.canRunInBackground = False

    def getParameterInfo(self):
        p_files = arcpy.Parameter(
            displayName="Input NetCDF files (browse to disk)",
            name="in_nc_files", datatype="DEFile",
            parameterType="Optional", direction="Input", multiValue=True,
        )
        p_files.filter.list = ["nc"]

        p_layers = arcpy.Parameter(
            displayName="OR input multidimensional raster layers (from map)",
            name="in_nc_layers", datatype="GPRasterLayer",
            parameterType="Optional", direction="Input", multiValue=True,
        )
        return [p_files, p_layers]

    def isLicensed(self):
        return True

    def updateMessages(self, parameters):
        if not parameters[0].value and not parameters[1].value:
            parameters[0].setErrorMessage("Provide at least one input source.")

    def execute(self, parameters, messages):
        file_raw = parameters[0].valueAsText
        layer_raw = parameters[1].valueAsText
        nc_files = resolve_nc_inputs(file_raw, layer_raw, messages)

        for nc in nc_files:
            messages.addMessage("=" * 70)
            messages.addMessage(os.path.basename(nc))
            messages.addMessage("=" * 70)

            try:
                ds = netCDF4.Dataset(nc)
                try:
                    # Global attrs
                    messages.addMessage("Global attributes:")
                    for attr in ds.ncattrs()[:5]:  # top 5 to avoid spam
                        val = str(ds.getncattr(attr))[:100]
                        messages.addMessage("  {} = {}".format(attr, val))

                    # Dimensions
                    messages.addMessage("Dimensions:")
                    for dname, d in ds.dimensions.items():
                        messages.addMessage("  {}: size={}".format(dname, d.size))

                    # Variables summary
                    messages.addMessage("Variables:")
                    for vname, v in ds.variables.items():
                        units = getattr(v, "units", "")
                        std_name = getattr(v, "standard_name", "")
                        long_name = getattr(v, "long_name", "")
                        messages.addMessage(
                            "  {}: dims={}, shape={}, units='{}', std_name='{}'".format(
                                vname, v.dimensions, v.shape, units, std_name
                            )
                        )
                        if long_name and long_name != vname:
                            messages.addMessage("      long_name: '{}'".format(long_name))

                    # Time decoding check
                    time_dim = find_dim(ds, TIME_DIM_CANDIDATES)
                    if time_dim:
                        messages.addMessage("Time decoding ({}):".format(time_dim))
                        try:
                            dates = get_years_from_time(ds, time_dim)
                            messages.addMessage(
                                "  {} timesteps: {:04d}-{:02d}-{:02d} to {:04d}-{:02d}-{:02d}".format(
                                    len(dates),
                                    dates[0][0], dates[0][1], dates[0][2],
                                    dates[-1][0], dates[-1][1], dates[-1][2]
                                )
                            )
                        except Exception as e:
                            messages.addWarningMessage(
                                "  Time decoding failed: {}".format(e)
                            )

                    # Detected data variable
                    try:
                        v = detect_data_variable(ds)
                        meta = get_var_metadata(ds, v)
                        messages.addMessage("Detected data variable: '{}'".format(v))
                        messages.addMessage("  standard_name: '{}'".format(meta["standard_name"]))
                        messages.addMessage("  long_name: '{}'".format(meta["long_name"]))
                        messages.addMessage("  units: '{}'".format(meta["units"]))
                        messages.addMessage("  fill_value: {}".format(meta["fill_value"]))
                    except Exception as e:
                        messages.addWarningMessage(
                            "  Auto-detect failed: {}".format(e)
                        )

                finally:
                    ds.close()
            except Exception as e:
                messages.addErrorMessage("Inspection failed: {}".format(e))


# ---------------------------------------------------------------------------
# Tool 3: Compute climate indicators
# ---------------------------------------------------------------------------

class ComputeERA5LandIndicators(object):
    def __init__(self):
        self.label = "3. Compute climate indicators"
        self.description = (
            "Compute annual climate indicators and climatology from per-day "
            "GeoTIFFs produced by Tool 1. Input is the Tool 1 output root; "
            "Tool 3 reads manifest.csv to discover available "
            "(variable, statistic, year) combinations. Always computes "
            "annual outputs; climatology (mean of annual values across the "
            "full year range) is written when enabled. Kelvin temperature "
            "variables are auto-converted to Celsius before metric "
            "computation. Output layout: "
            "annual/<variable>/<statistic>/<metric>/, "
            "climatology/<variable>/<statistic>/."
        )
        self.canRunInBackground = False

    def getParameterInfo(self):
        p_input = arcpy.Parameter(
            displayName="Input folder (Tool 1 output root with manifest.csv)",
            name="in_folder", datatype="DEFolder",
            parameterType="Required", direction="Input",
        )

        p_output = arcpy.Parameter(
            displayName="Output folder",
            name="out_folder", datatype="DEFolder",
            parameterType="Required", direction="Input",
        )

        p_variable = arcpy.Parameter(
            displayName="Variable",
            name="variable", datatype="GPString",
            parameterType="Required", direction="Input",
        )
        p_variable.filter.type = "ValueList"
        p_variable.filter.list = []

        p_statistic = arcpy.Parameter(
            displayName="Statistic",
            name="statistic", datatype="GPString",
            parameterType="Required", direction="Input",
        )
        p_statistic.filter.type = "ValueList"
        p_statistic.filter.list = []

        p_metrics = arcpy.Parameter(
            displayName="Metrics to compute",
            name="metrics", datatype="GPString",
            parameterType="Required", direction="Input",
            multiValue=True,
        )
        p_metrics.filter.type = "ValueList"
        p_metrics.filter.list = sorted(METRIC_REGISTRY.keys())

        p_generic_threshold = arcpy.Parameter(
            displayName="Generic threshold (for count/seconds/accumulation/consecutive metrics)",
            name="generic_threshold", datatype="GPDouble",
            parameterType="Optional", direction="Input",
        )
        p_generic_threshold.value = 15.0

        p_hdd_base = arcpy.Parameter(
            displayName="HDD base temperature (degC, Spinoni default 15.5)",
            name="hdd_base", datatype="GPDouble",
            parameterType="Optional", direction="Input",
        )
        p_hdd_base.value = 15.5

        p_cdd_base = arcpy.Parameter(
            displayName="CDD base temperature (degC, Spinoni default 22.0)",
            name="cdd_base", datatype="GPDouble",
            parameterType="Optional", direction="Input",
        )
        p_cdd_base.value = 22.0

        p_gdd_base = arcpy.Parameter(
            displayName="GDD base temperature (degC, default 10.0)",
            name="gdd_base", datatype="GPDouble",
            parameterType="Optional", direction="Input",
        )
        p_gdd_base.value = 10.0

        p_precip_wet = arcpy.Parameter(
            displayName="Precipitation wet-day threshold (mm, default 1.0)",
            name="precip_wet_threshold", datatype="GPDouble",
            parameterType="Optional", direction="Input",
        )
        p_precip_wet.value = 1.0

        p_pv_pr = arcpy.Parameter(
            displayName="PV performance ratio (default 0.75)",
            name="pv_performance_ratio", datatype="GPDouble",
            parameterType="Optional", direction="Input",
        )
        p_pv_pr.value = 0.75

        p_rad_threshold = arcpy.Parameter(
            displayName="Radiation threshold for above-threshold count (W/m^2, default 100.0)",
            name="radiation_threshold_w_per_m2", datatype="GPDouble",
            parameterType="Optional", direction="Input",
        )
        p_rad_threshold.value = 100.0

        p_wind_rho = arcpy.Parameter(
            displayName="Wind air density rho (kg/m^3, default 1.225)",
            name="wind_air_density", datatype="GPDouble",
            parameterType="Optional", direction="Input",
        )
        p_wind_rho.value = 1.225

        p_wind_calm = arcpy.Parameter(
            displayName="Wind calm threshold (m/s, default 2.0)",
            name="wind_calm_threshold", datatype="GPDouble",
            parameterType="Optional", direction="Input",
        )
        p_wind_calm.value = 2.0

        p_year_min = arcpy.Parameter(
            displayName="Year minimum (inclusive, optional)",
            name="year_min", datatype="GPLong",
            parameterType="Optional", direction="Input",
        )

        p_year_max = arcpy.Parameter(
            displayName="Year maximum (inclusive, optional)",
            name="year_max", datatype="GPLong",
            parameterType="Optional", direction="Input",
        )

        p_clip = arcpy.Parameter(
            displayName="Clip polygon feature class (optional)",
            name="clip_fc", datatype="GPFeatureLayer",
            parameterType="Optional", direction="Input",
        )
        p_clip.filter.list = ["Polygon"]

        p_climatology = arcpy.Parameter(
            displayName="Compute climatology (mean of annual values)",
            name="compute_climatology", datatype="GPBoolean",
            parameterType="Optional", direction="Input",
        )
        p_climatology.value = True

        p_ipma_normals = arcpy.Parameter(
            displayName="Emulate IPMA normals (1981-2010 and 1991-2020)",
            name="compute_ipma_normals", datatype="GPBoolean",
            parameterType="Optional", direction="Input",
        )
        p_ipma_normals.value = False

        p_manifest = arcpy.Parameter(
            displayName="Write manifest CSV",
            name="write_manifest", datatype="GPBoolean",
            parameterType="Optional", direction="Input",
        )
        p_manifest.value = True

        p_verbose = arcpy.Parameter(
            displayName="Verbose logging",
            name="verbose", datatype="GPBoolean",
            parameterType="Optional", direction="Input",
        )
        p_verbose.value = False

        # G.POT-oriented heating/cooling season thresholds. Appended at the
        # END of the parameter list (not next to HDD/CDD bases) so the
        # positional indices of all pre-existing parameters in execute()
        # do not shift. UI ordering is cosmetic; correctness is positional.
        p_heating_threshold = arcpy.Parameter(
            displayName="Heating season threshold (degC, default 15.5)",
            name="heating_threshold", datatype="GPDouble",
            parameterType="Optional", direction="Input",
        )
        p_heating_threshold.value = 15.5

        p_cooling_threshold = arcpy.Parameter(
            displayName="Cooling season threshold (degC, default 22.0)",
            name="cooling_threshold", datatype="GPDouble",
            parameterType="Optional", direction="Input",
        )
        p_cooling_threshold.value = 22.0

        # WMO standard normals: 1961-1990, 1981-2010, 1991-2020.
        # Appended at the end of the return list for the same index-stability
        # reason as the heating/cooling thresholds. Disjoint from IPMA
        # (1981-2010 + 1991-2020); enabling both produces duplicate writes
        # of the overlapping periods into ipma_normals/ and wmo_normals/.
        p_wmo_normals = arcpy.Parameter(
            displayName="Compute WMO standard normals (1961-1990, 1981-2010, 1991-2020)",
            name="compute_wmo_normals", datatype="GPBoolean",
            parameterType="Optional", direction="Input",
        )
        p_wmo_normals.value = False

        # Parameter return order is what drives the UI rendering. Grouping:
        #  - input/output (0-1)
        #  - data selection (2-4): variable, statistic, metrics
        #  - thresholds and bases (5-15), all together
        #  - year range (16-17)
        #  - spatial (18): clip polygon
        #  - climatology checkboxes (19-21), all together
        #  - output flags (22-23): manifest, verbose
        # `execute()` and `updateParameters()` read by positional index;
        # any reorder requires synchronised updates there.
        return [
            p_input, p_output,                                          # 0-1
            p_variable, p_statistic, p_metrics,                         # 2-4
            p_generic_threshold,                                        # 5
            p_hdd_base, p_cdd_base, p_gdd_base,                         # 6-8
            p_heating_threshold, p_cooling_threshold,                   # 9-10
            p_precip_wet,                                               # 11
            p_pv_pr, p_rad_threshold,                                   # 12-13
            p_wind_rho, p_wind_calm,                                    # 14-15
            p_year_min, p_year_max,                                     # 16-17
            p_clip,                                                     # 18
            p_climatology, p_ipma_normals, p_wmo_normals,               # 19-21
            p_manifest, p_verbose,                                      # 22-23
        ]

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        """Auto-populate variable and statistic dropdowns from the input
        folder. Lookups are silent on error (do not block the UI)."""
        p_input = parameters[0]
        p_variable = parameters[2]
        p_statistic = parameters[3]

        if not p_input.value:
            return

        in_path = p_input.valueAsText
        index = None
        try:
            result, _ = discover_inputs_from_manifest(in_path)
            if result and result["index"]:
                index = result["index"]
            else:
                result, _ = discover_inputs_from_walk(in_path)
                if result and result["index"]:
                    index = result["index"]
        except Exception:
            return

        if index is None:
            return

        variables = sorted({v for v, _ in index.keys()})
        if p_variable.filter.list != variables:
            p_variable.filter.list = variables

        if p_variable.value:
            current_var = p_variable.valueAsText
            stats = sorted({s for v, s in index.keys() if v == current_var})
            if p_statistic.filter.list != stats:
                p_statistic.filter.list = stats

            # Filter metrics list to those applicable to the selected
            # variable's type. Metrics with applicable_to=None are generic
            # (apply to any variable). Reduces 48 entries to ~25 (typical
            # for temperature-only runs); avoids the user accidentally
            # selecting `wind_power_density` for a t2m run.
            var_type = variable_type_of(current_var)
            filtered_metrics = [
                name for name in sorted(METRIC_REGISTRY.keys())
                if (METRIC_REGISTRY[name].get("applicable_to") is None
                    or var_type in METRIC_REGISTRY[name].get(
                        "applicable_to", set()
                    ))
            ]
            if parameters[4].filter.list != filtered_metrics:
                parameters[4].filter.list = filtered_metrics

        # Dynamic parameter enable/disable: only metrics actually selected
        # drive which advanced threshold/base params are interactable.
        # Reduces UI noise for single-domain runs (e.g. only temperature).
        # The check is string-based on the multi-select raw value because
        # parameters[4].values returns None until the user opens the picker
        # at least once on some ArcGIS Pro versions.
        selected_raw = (parameters[4].valueAsText or "").lower()
        needs = {
            "generic":   any(k in selected_raw for k in (
                "count_above_threshold", "count_below_threshold",
                "seconds_above_threshold", "seconds_below_threshold",
                "accumulation_above_threshold", "accumulation_below_threshold",
                "consecutive_above_threshold", "consecutive_below_threshold",
            )),
            "hdd":       "hdd_spinoni" in selected_raw,
            "cdd":       "cdd_spinoni" in selected_raw,
            "gdd":       "gdd" in selected_raw,
            "precip":    any(k in selected_raw for k in (
                "wet_days", "r10mm", "r20mm", "rx1day", "r95p",
            )),
            "pv":        "solar_pv" in selected_raw,
            "rad":       any(k in selected_raw for k in (
                "rad_", "solar_", "thermal_",
            )),
            "wind":      "wind_" in selected_raw,
            "heating":   "heating_season_length" in selected_raw,
            "cooling":   "cooling_season_length" in selected_raw,
        }
        # Positional indices defined in getParameterInfo's return list
        # (grouped: thresholds 5-15, climatology checkboxes 19-21).
        try:
            parameters[5].enabled  = needs["generic"]   # generic_threshold
            parameters[6].enabled  = needs["hdd"]       # hdd_base
            parameters[7].enabled  = needs["cdd"]       # cdd_base
            parameters[8].enabled  = needs["gdd"]       # gdd_base
            parameters[9].enabled  = needs["heating"]   # heating_threshold
            parameters[10].enabled = needs["cooling"]   # cooling_threshold
            parameters[11].enabled = needs["precip"]    # precip_wet_threshold
            parameters[12].enabled = needs["pv"]        # pv_performance_ratio
            parameters[13].enabled = needs["rad"]       # radiation_threshold_w_per_m2
            parameters[14].enabled = needs["wind"]      # wind_air_density
            parameters[15].enabled = needs["wind"]      # wind_calm_threshold
        except IndexError:
            # Defensive: if the parameter list shape changes in the future,
            # fail open (all enabled) rather than crash the UI.
            pass

    def updateMessages(self, parameters):
        pass

    @staticmethod
    def _eta_label(idx_year, total_years, year, pipeline_label,
                    pipeline_start, years_done):
        """Build a progressor label including a rolling ETA after the first
        year has finished (so we have empirical rate)."""
        base = "{} year {}/{}: {}".format(
            pipeline_label, idx_year, total_years, year
        )
        if years_done <= 0:
            return base
        elapsed = _time.time() - pipeline_start
        avg = elapsed / years_done
        remaining = avg * (total_years - years_done)
        if remaining < 60:
            return "{} (ETA {:.0f}s)".format(base, remaining)
        return "{} (ETA {:.0f}m{:02.0f}s)".format(
            base, remaining // 60, remaining % 60
        )

    def _emit_climatology(self, buffer, variable, statistic_label,
                           sorted_years, working_units,
                           out_folder, mask, gt, srs_wkt,
                           do_climatology, do_ipma, log,
                           do_wmo=False):
        """Run climatology and IPMA / WMO normal aggregation on an annual buffer.

        Single point of truth used by the regular, DTR and wind pipelines.
        statistic_label is the string written into output paths and manifest
        ('daily_mean', 'daily_maximum', ... for regular; 'dtr' for DTR;
        'derived' for wind).

        Aggregation outputs:
          - do_climatology: one TIFF per metric covering year_min..year_max
            of the input (folder: climatology/).
          - do_ipma: 1981-2010 + 1991-2020 (folder: ipma_normals/).
          - do_wmo:  1961-1990 + 1981-2010 + 1991-2020 (folder: wmo_normals/).
            Periods overlapping IPMA produce duplicate writes if both flags
            are on; this is intentional (per-checkbox separation, the user
            chose to enable both).

        Returns the list of records for written climatology / normal TIFFs.
        """
        if not buffer:
            return []
        records = []

        if do_climatology:
            start_y, end_y = sorted_years[0], sorted_years[-1]
            log("--- climatology {}-{} ({}) ---".format(
                start_y, end_y, statistic_label
            ))
            for metric_name, buf in buffer.items():
                spec = METRIC_REGISTRY.get(metric_name, {})
                units_out_spec = spec.get("units_out", _passthrough_units)
                units_out = (units_out_spec(working_units)
                             if callable(units_out_spec)
                             else str(units_out_spec))
                rec = process_climatology(
                    variable, statistic_label, metric_name, buf,
                    out_folder, mask, gt, srs_wkt, units_out,
                    period_kind="climatology",
                    period_label="",
                    start_year=start_y, end_year=end_y,
                )
                if rec:
                    records.append(rec)
                    log("  wrote climatology for '{}'".format(metric_name))

        def _emit_period_set(periods, period_kind, label_for_log):
            """Helper: iterate (start_year, end_year) tuples and write one
            normal per metric per period to <period_kind>s/ subfolder."""
            for ys, ye in periods:
                log("--- {} {}-{} ({}) ---".format(
                    label_for_log, ys, ye, statistic_label
                ))
                period_years = [y for y in sorted_years if ys <= y <= ye]
                expected = ye - ys + 1
                if not period_years:
                    log("  no years available in {}-{}, {} skipped".format(
                        ys, ye, label_for_log
                    ), "warn")
                    continue
                if len(period_years) < expected:
                    log("  partial coverage: only {} of {} years available "
                        "for {}-{}. Normal computed on available years.".format(
                            len(period_years), expected, ys, ye
                        ), "warn")
                for metric_name, buf in buffer.items():
                    buf_period = [(y, a) for (y, a) in buf if ys <= y <= ye]
                    if not buf_period:
                        continue
                    spec = METRIC_REGISTRY.get(metric_name, {})
                    units_out_spec = spec.get("units_out", _passthrough_units)
                    units_out = (units_out_spec(working_units)
                                 if callable(units_out_spec)
                                 else str(units_out_spec))
                    rec = process_climatology(
                        variable, statistic_label, metric_name, buf_period,
                        out_folder, mask, gt, srs_wkt, units_out,
                        period_kind=period_kind,
                        period_label="{}-{}".format(ys, ye),
                        start_year=ys, end_year=ye,
                    )
                    if rec:
                        records.append(rec)
                        log("  wrote {} {}-{} for '{}'".format(
                            label_for_log, ys, ye, metric_name
                        ))

        if do_ipma:
            _emit_period_set(
                [(1981, 2010), (1991, 2020)],
                period_kind="ipma_normal",
                label_for_log="IPMA normal",
            )
        if do_wmo:
            _emit_period_set(
                [(1961, 1990), (1981, 2010), (1991, 2020)],
                period_kind="wmo_normal",
                label_for_log="WMO normal",
            )
        return records

    def execute(self, parameters, messages):
        # Index map (mirrors the getParameterInfo return list grouping):
        #  0-1   input/output folders
        #  2-4   variable, statistic, metrics
        #  5     generic_threshold
        #  6-8   hdd_base, cdd_base, gdd_base
        #  9-10  heating_threshold, cooling_threshold
        #  11    precip_wet_threshold
        #  12-13 pv_performance_ratio, radiation_threshold_w_per_m2
        #  14-15 wind_air_density, wind_calm_threshold
        #  16-17 year_min, year_max
        #  18    clip_fc
        #  19-21 compute_climatology, compute_ipma_normals, compute_wmo_normals
        #  22-23 write_manifest, verbose
        in_folder = parameters[0].valueAsText
        out_folder = parameters[1].valueAsText
        variable = parameters[2].valueAsText
        statistic = parameters[3].valueAsText
        metrics_raw = parameters[4].valueAsText
        generic_threshold = (parameters[5].value
                             if parameters[5].value is not None else 15.0)
        hdd_base = (parameters[6].value
                    if parameters[6].value is not None else 15.5)
        cdd_base = (parameters[7].value
                    if parameters[7].value is not None else 22.0)
        gdd_base = (parameters[8].value
                    if parameters[8].value is not None else 10.0)
        heating_threshold = (parameters[9].value
                             if parameters[9].value is not None else 15.5)
        cooling_threshold = (parameters[10].value
                             if parameters[10].value is not None else 22.0)
        precip_wet_threshold = (parameters[11].value
                                if parameters[11].value is not None else 1.0)
        pv_performance_ratio = (parameters[12].value
                                if parameters[12].value is not None else 0.75)
        radiation_threshold_w_per_m2 = (parameters[13].value
                                        if parameters[13].value is not None
                                        else 100.0)
        wind_air_density = (parameters[14].value
                            if parameters[14].value is not None else 1.225)
        wind_calm_threshold = (parameters[15].value
                               if parameters[15].value is not None else 2.0)
        year_min = parameters[16].value
        year_max = parameters[17].value
        clip_fc = parameters[18].valueAsText if parameters[18].value else None
        compute_climatology_flag = bool(parameters[19].value)
        compute_ipma_normals = bool(parameters[20].value)
        compute_wmo_normals = bool(parameters[21].value)
        write_manifest = bool(parameters[22].value)
        verbose = bool(parameters[23].value)

        metrics = [m.strip().strip("'\"") for m in (metrics_raw or "").split(";")
                   if m.strip()]
        if not metrics:
            messages.addErrorMessage("Select at least one metric.")
            return

        # Separate metrics requiring multi-input handling (DTR, wind) from
        # regular metrics that run on the single (variable, statistic) stack
        # the user selected.
        dtr_metrics = [
            m for m in metrics
            if METRIC_REGISTRY.get(m, {}).get("multi_input") == "dtr"
        ]
        wind_metrics_sel = [
            m for m in metrics
            if METRIC_REGISTRY.get(m, {}).get("multi_input") == "wind"
        ]
        regular_metrics = [
            m for m in metrics
            if METRIC_REGISTRY.get(m, {}).get("multi_input") is None
        ]

        if not os.path.isdir(out_folder):
            os.makedirs(out_folder)

        _MASK_CACHE.clear()

        def log(msg, level="msg"):
            if level == "warn":
                messages.addWarningMessage(msg)
            elif level == "err":
                messages.addErrorMessage(msg)
            else:
                messages.addMessage(msg)

        log("Tool 3: Compute climate indicators")
        log("Input:  {}".format(in_folder))
        log("Output: {}".format(out_folder))

        try:
            idx = build_inputs_index(in_folder, log)
        except RuntimeError as e:
            messages.addErrorMessage(str(e))
            return

        # ---------------- Resolve regular pipeline inputs ----------------
        # Failures here demote regular metrics to skipped (rather than abort
        # the whole tool) so DTR / wind pipelines can still run on the same
        # invocation.
        regular_year_files = {}
        if regular_metrics:
            key = (variable, statistic)
            if key not in idx["index"]:
                log("Regular pipeline: no inputs for variable='{}', "
                    "statistic='{}'. {} regular metric(s) will be skipped. "
                    "Available (var,stat) keys: {}".format(
                        variable, statistic, len(regular_metrics),
                        sorted(idx["index"].keys())
                    ), "warn")
                regular_metrics = []
            else:
                yf = dict(idx["index"][key])
                if year_min is not None:
                    yf = {y: f for y, f in yf.items() if y >= year_min}
                if year_max is not None:
                    yf = {y: f for y, f in yf.items() if y <= year_max}
                if not yf:
                    log("Regular pipeline: no years for ({}, {}) in range "
                        "[{}, {}]. {} regular metric(s) will be skipped.".format(
                            variable, statistic, year_min, year_max,
                            len(regular_metrics)
                        ), "warn")
                    regular_metrics = []
                else:
                    regular_year_files = yf

        # ---------------- Resolve WIND pipeline inputs ----------------
        wind_year_files = {}  # {year: (u_paths, v_paths)}
        if wind_metrics_sel:
            # Wind always pairs u10 + v10. Statistic comes from p_statistic
            # (typically daily_mean for ERA5-Land wind components).
            key_u = ("u10", statistic)
            key_v = ("v10", statistic)
            missing = [k for k in (key_u, key_v) if k not in idx["index"]]
            if missing:
                log("WIND requested but missing required inputs (u10, v10) "
                    "at statistic '{}': {}. Wind metrics will be skipped.".format(
                        statistic, [k[0] for k in missing]
                    ), "warn")
                wind_metrics_sel = []
            else:
                yu = set(idx["index"][key_u])
                yv = set(idx["index"][key_v])
                common = yu & yv
                if year_min is not None:
                    common = {y for y in common if y >= year_min}
                if year_max is not None:
                    common = {y for y in common if y <= year_max}
                if not common:
                    log("WIND: no years where both u10 and v10 are available "
                        "at statistic '{}'.".format(statistic), "warn")
                    wind_metrics_sel = []
                else:
                    for y in common:
                        wind_year_files[y] = (
                            sorted(idx["index"][key_u][y]),
                            sorted(idx["index"][key_v][y]),
                        )

        # ---------------- Resolve DTR pipeline inputs ----------------
        dtr_year_files = {}  # {year: (tmax_paths, tmin_paths)}
        if dtr_metrics:
            key_max = (variable, "daily_maximum")
            key_min = (variable, "daily_minimum")
            missing = [k for k in (key_max, key_min) if k not in idx["index"]]
            if missing:
                log("DTR requested but missing statistics for '{}': {}. "
                    "DTR metrics will be skipped.".format(
                        variable, [k[1] for k in missing]
                    ), "warn")
                dtr_metrics = []
            else:
                ymax = set(idx["index"][key_max])
                ymin = set(idx["index"][key_min])
                common = ymax & ymin
                if year_min is not None:
                    common = {y for y in common if y >= year_min}
                if year_max is not None:
                    common = {y for y in common if y <= year_max}
                if not common:
                    log("DTR: no years where both daily_maximum and "
                        "daily_minimum are available for variable '{}'.".format(
                            variable
                        ), "warn")
                    dtr_metrics = []
                else:
                    for y in common:
                        dtr_year_files[y] = (
                            sorted(idx["index"][key_max][y]),
                            sorted(idx["index"][key_min][y]),
                        )

        if (not regular_year_files
                and not dtr_year_files
                and not wind_year_files):
            messages.addErrorMessage(
                "No work to do after input resolution. Check variable, "
                "statistic, year range, and DTR/Wind pre-requisites."
            )
            return

        # ---------------- Reference TIFF + grid metadata ----------------
        if regular_year_files:
            sample_tif = regular_year_files[sorted(regular_year_files)[0]][0]
        elif dtr_year_files:
            sample_tif = dtr_year_files[sorted(dtr_year_files)[0]][0][0]
        else:
            sample_tif = wind_year_files[sorted(wind_year_files)[0]][0][0]

        gt, srs_wkt, nrows, ncols = read_grid_metadata(sample_tif)
        log("Grid: {} rows x {} cols".format(nrows, ncols))
        if verbose:
            log("  geotransform: {}".format(gt))

        # Cross-pipeline grid consistency: DTR and wind pipelines may read
        # TIFFs from a different (variable, statistic) tree than the primary
        # sample. If those have a different grid (rare but possible when the
        # user mixed Tool 1 runs with different clip polygons), skip the
        # affected pipeline rather than write outputs with wrong geotransform.
        def _grid_consistent(other_tif):
            try:
                ogt, _, onrows, oncols = read_grid_metadata(other_tif)
            except RuntimeError:
                return False
            if (onrows, oncols) != (nrows, ncols):
                return False
            return all(abs(a - b) < 1e-6 for a, b in zip(ogt, gt))

        if dtr_year_files:
            dtr_sample = dtr_year_files[sorted(dtr_year_files)[0]][0][0]
            if dtr_sample != sample_tif and not _grid_consistent(dtr_sample):
                log("DTR pipeline grid mismatches primary sample; DTR "
                    "metrics will be skipped to avoid wrong georeferencing.",
                    "warn")
                dtr_year_files = {}
                dtr_metrics = []

        if wind_year_files:
            wind_sample = wind_year_files[sorted(wind_year_files)[0]][0][0]
            if wind_sample != sample_tif and not _grid_consistent(wind_sample):
                log("WIND pipeline grid mismatches primary sample; wind "
                    "metrics will be skipped to avoid wrong georeferencing.",
                    "warn")
                wind_year_files = {}
                wind_metrics_sel = []

        if (not regular_year_files
                and not dtr_year_files
                and not wind_year_files):
            messages.addErrorMessage(
                "All pipelines disabled after grid consistency / input checks."
            )
            return

        mask = None
        if clip_fc:
            lat, lon = derive_lat_lon_from_gt(gt, nrows, ncols)
            log("Clip polygon: {}".format(clip_fc))
            mask = get_or_compute_mask(clip_fc, lat, lon, lat_ascending=False)

        input_units = idx["units"].get(variable, "")
        params = {
            "generic_threshold": float(generic_threshold),
            "hdd_base": float(hdd_base),
            "cdd_base": float(cdd_base),
            "gdd_base": float(gdd_base),
            "precip_wet_threshold": float(precip_wet_threshold),
            "pv_performance_ratio": float(pv_performance_ratio),
            "radiation_threshold_w_per_m2": float(radiation_threshold_w_per_m2),
            "wind_air_density": float(wind_air_density),
            "wind_calm_threshold": float(wind_calm_threshold),
            "heating_threshold": float(heating_threshold),
            "cooling_threshold": float(cooling_threshold),
        }
        var_type = variable_type_of(variable)
        log("Variable: {} (type: {})".format(variable, var_type))
        log("Statistic: {}".format(statistic))
        log("Metrics regular: {}".format(", ".join(regular_metrics) or "(none)"))
        log("Metrics DTR:     {}".format(", ".join(dtr_metrics) or "(none)"))
        log("Metrics WIND:    {}".format(", ".join(wind_metrics_sel) or "(none)"))
        log("Params: generic_threshold={}, hdd_base={}, cdd_base={}, "
            "gdd_base={}, precip_wet_threshold={}, "
            "pv_performance_ratio={}, radiation_threshold_w_per_m2={}, "
            "wind_air_density={}, wind_calm_threshold={}, "
            "heating_threshold={}, cooling_threshold={}".format(
                generic_threshold, hdd_base, cdd_base, gdd_base,
                precip_wet_threshold,
                pv_performance_ratio, radiation_threshold_w_per_m2,
                wind_air_density, wind_calm_threshold,
                heating_threshold, cooling_threshold
            ))

        # Validation warnings (variable applicability, expected statistic)
        for m in metrics:
            spec = METRIC_REGISTRY.get(m)
            if spec is None:
                continue
            applicable = spec.get("applicable_to")
            if applicable is not None and var_type not in applicable:
                log("  WARNING: metric '{}' expects variable type in {} but "
                    "variable '{}' is type '{}'. Output may not be meaningful.".format(
                        m, sorted(applicable), variable, var_type
                    ), "warn")
            expected_stat = spec.get("expects_statistic")
            if (expected_stat and statistic
                    and spec.get("multi_input") is None
                    and statistic != expected_stat):
                log("  WARNING: metric '{}' expects statistic '{}' but input "
                    "is '{}'. Output may be incorrect.".format(
                        m, expected_stat, statistic
                    ), "warn")

        working_units = resolve_working_units(variable, input_units)
        all_records = []
        t_start = _time.time()

        # ---------------- Regular pipeline ----------------
        if regular_year_files and regular_metrics:
            sorted_years = sorted(regular_year_files.keys())
            log("Regular pipeline: years {} to {} ({} years)".format(
                sorted_years[0], sorted_years[-1], len(sorted_years)
            ))
            arcpy.SetProgressor(
                type="step", message="Regular: processing years...",
                min_range=0, max_range=max(len(sorted_years), 1), step_value=1,
            )
            regular_buffer = {}
            pipeline_start = _time.time()
            years_done = 0
            try:
                for idx_year, year in enumerate(sorted_years, 1):
                    arcpy.SetProgressorLabel(self._eta_label(
                        idx_year, len(sorted_years), year,
                        "Regular", pipeline_start, years_done,
                    ))
                    log("--- regular year {} ({}/{}) ---".format(
                        year, idx_year, len(sorted_years)
                    ))
                    try:
                        annual_nan, recs = process_year_metrics(
                            variable, statistic, year,
                            regular_year_files[year],
                            regular_metrics, out_folder, mask, gt, srs_wkt,
                            input_units, params, log,
                        )
                        all_records.extend(recs)
                        for m, arr in annual_nan.items():
                            regular_buffer.setdefault(m, []).append((year, arr))
                        years_done += 1
                    except Exception as e:
                        log("FAILED regular year {}: {}".format(year, e), "warn")
                    finally:
                        arcpy.SetProgressorPosition(idx_year)
            finally:
                arcpy.ResetProgressor()

            all_records.extend(self._emit_climatology(
                regular_buffer, variable, statistic, sorted_years,
                working_units, out_folder, mask, gt, srs_wkt,
                compute_climatology_flag, compute_ipma_normals, log,
                do_wmo=compute_wmo_normals,
            ))

        # ---------------- DTR pipeline ----------------
        if dtr_year_files and dtr_metrics:
            sorted_dtr = sorted(dtr_year_files.keys())
            log("DTR pipeline: years {} to {} ({} years)".format(
                sorted_dtr[0], sorted_dtr[-1], len(sorted_dtr)
            ))
            arcpy.SetProgressor(
                type="step", message="DTR: processing years...",
                min_range=0, max_range=max(len(sorted_dtr), 1), step_value=1,
            )
            dtr_buffer = {}
            pipeline_start = _time.time()
            years_done = 0
            try:
                for idx_year, year in enumerate(sorted_dtr, 1):
                    arcpy.SetProgressorLabel(self._eta_label(
                        idx_year, len(sorted_dtr), year,
                        "DTR", pipeline_start, years_done,
                    ))
                    log("--- DTR year {} ({}/{}) ---".format(
                        year, idx_year, len(sorted_dtr)
                    ))
                    tmax_paths, tmin_paths = dtr_year_files[year]
                    try:
                        annual_nan, recs = process_year_dtr(
                            variable, year, tmax_paths, tmin_paths,
                            dtr_metrics, out_folder, mask, gt, srs_wkt, log,
                        )
                        all_records.extend(recs)
                        for m, arr in annual_nan.items():
                            dtr_buffer.setdefault(m, []).append((year, arr))
                        years_done += 1
                    except Exception as e:
                        log("FAILED DTR year {}: {}".format(year, e), "warn")
                    finally:
                        arcpy.SetProgressorPosition(idx_year)
            finally:
                arcpy.ResetProgressor()

            all_records.extend(self._emit_climatology(
                dtr_buffer, variable, "dtr", sorted_dtr,
                "degC", out_folder, mask, gt, srs_wkt,
                compute_climatology_flag, compute_ipma_normals, log,
                do_wmo=compute_wmo_normals,
            ))

        # ---------------- WIND pipeline ----------------
        if wind_year_files and wind_metrics_sel:
            sorted_wind = sorted(wind_year_files.keys())
            log("WIND pipeline: years {} to {} ({} years)".format(
                sorted_wind[0], sorted_wind[-1], len(sorted_wind)
            ))
            arcpy.SetProgressor(
                type="step", message="Wind: processing years...",
                min_range=0, max_range=max(len(sorted_wind), 1), step_value=1,
            )
            wind_buffer = {}
            pipeline_start = _time.time()
            years_done = 0
            try:
                for idx_year, year in enumerate(sorted_wind, 1):
                    arcpy.SetProgressorLabel(self._eta_label(
                        idx_year, len(sorted_wind), year,
                        "Wind", pipeline_start, years_done,
                    ))
                    log("--- WIND year {} ({}/{}) ---".format(
                        year, idx_year, len(sorted_wind)
                    ))
                    u_paths, v_paths = wind_year_files[year]
                    try:
                        annual_nan, recs = process_year_wind(
                            year, u_paths, v_paths, wind_metrics_sel,
                            out_folder, mask, gt, srs_wkt, params, log,
                        )
                        all_records.extend(recs)
                        for m, arr in annual_nan.items():
                            wind_buffer.setdefault(m, []).append((year, arr))
                        years_done += 1
                    except Exception as e:
                        log("FAILED wind year {}: {}".format(year, e), "warn")
                    finally:
                        arcpy.SetProgressorPosition(idx_year)
            finally:
                arcpy.ResetProgressor()

            all_records.extend(self._emit_climatology(
                wind_buffer, "wind", "derived", sorted_wind,
                "m/s", out_folder, mask, gt, srs_wkt,
                compute_climatology_flag, compute_ipma_normals, log,
                do_wmo=compute_wmo_normals,
            ))

        if write_manifest and all_records:
            manifest_path = os.path.join(out_folder, "manifest.csv")
            fields = ["raster", "variable", "statistic", "metric", "category",
                      "period_kind", "year", "start_year", "end_year", "units"]
            with open(manifest_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=fields)
                writer.writeheader()
                writer.writerows(all_records)
            log("Manifest: {}".format(manifest_path))

        elapsed = _time.time() - t_start
        log("Done. {} TIFFs written in {:.1f}s.".format(
            len(all_records), elapsed
        ))


# ---------------------------------------------------------------------------
# Tool 4: Extract regional time series
# ---------------------------------------------------------------------------

class ExtractRegionalTimeSeries(object):
    """Zonal statistics over polygon regions for each TIFF produced by
    Tool 1 (daily values) or Tool 3 (annual indicators / climatology /
    IPMA normals). Emits a long-format CSV ready for downstream analysis
    in R / Pandas / SQL. Auto-detects input source from manifest.csv.
    """

    def __init__(self):
        self.label = "4. Extract regional time series"
        self.description = (
            "Compute zonal statistics (mean and optionally min, max, std, "
            "count) over polygon regions for each TIFF in a Tool 1 or "
            "Tool 3 output folder. Input source is auto-detected from "
            "manifest.csv (presence of 'metric' column => Tool 3). "
            "Output is a long-format CSV with one row per (region, "
            "variable, statistic, [metric,] period_kind, time_key, "
            "zonal_stat)."
        )
        self.canRunInBackground = False

    def getParameterInfo(self):
        p_input = arcpy.Parameter(
            displayName="Input folder (Tool 1 or Tool 3 output, with manifest.csv)",
            name="in_folder", datatype="DEFolder",
            parameterType="Required", direction="Input",
        )

        p_regions = arcpy.Parameter(
            displayName="Region polygon feature class",
            name="region_fc", datatype="GPFeatureLayer",
            parameterType="Required", direction="Input",
        )
        p_regions.filter.list = ["Polygon"]

        p_region_id = arcpy.Parameter(
            displayName="Region ID field (integer; SHORT/LONG/OBJECTID)",
            name="region_id_field", datatype="Field",
            parameterType="Required", direction="Input",
        )
        p_region_id.parameterDependencies = [p_regions.name]
        p_region_id.filter.list = ["Short", "Long", "OID"]

        p_region_label = arcpy.Parameter(
            displayName="Region label field (optional, text for human reading)",
            name="region_label_field", datatype="Field",
            parameterType="Optional", direction="Input",
        )
        p_region_label.parameterDependencies = [p_regions.name]
        p_region_label.filter.list = ["Text"]

        p_out_csv = arcpy.Parameter(
            displayName="Output CSV file",
            name="out_csv", datatype="DEFile",
            parameterType="Required", direction="Output",
        )
        p_out_csv.filter.list = ["csv"]

        p_variable = arcpy.Parameter(
            displayName="Variable filter (optional; default all)",
            name="variable", datatype="GPString",
            parameterType="Optional", direction="Input",
            multiValue=True,
        )
        p_variable.filter.type = "ValueList"
        p_variable.filter.list = []

        p_statistic = arcpy.Parameter(
            displayName="Statistic filter (optional; default all)",
            name="statistic", datatype="GPString",
            parameterType="Optional", direction="Input",
            multiValue=True,
        )
        p_statistic.filter.type = "ValueList"
        p_statistic.filter.list = []

        p_metric = arcpy.Parameter(
            displayName="Metric filter (Tool 3 only; optional, default all)",
            name="metric", datatype="GPString",
            parameterType="Optional", direction="Input",
            multiValue=True,
        )
        p_metric.filter.type = "ValueList"
        p_metric.filter.list = []

        p_year_min = arcpy.Parameter(
            displayName="Year minimum (inclusive, optional)",
            name="year_min", datatype="GPLong",
            parameterType="Optional", direction="Input",
        )

        p_year_max = arcpy.Parameter(
            displayName="Year maximum (inclusive, optional)",
            name="year_max", datatype="GPLong",
            parameterType="Optional", direction="Input",
        )

        p_stats = arcpy.Parameter(
            displayName="Zonal statistics to compute",
            name="zonal_stats", datatype="GPString",
            parameterType="Required", direction="Input",
            multiValue=True,
        )
        p_stats.filter.type = "ValueList"
        p_stats.filter.list = ["mean", "min", "max", "std", "count"]
        p_stats.values = ["mean"]

        p_verbose = arcpy.Parameter(
            displayName="Verbose logging",
            name="verbose", datatype="GPBoolean",
            parameterType="Optional", direction="Input",
        )
        p_verbose.value = False

        return [p_input, p_regions, p_region_id, p_region_label, p_out_csv,
                p_variable, p_statistic, p_metric,
                p_year_min, p_year_max,
                p_stats, p_verbose]

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        """Auto-populate variable / statistic / metric dropdowns from the
        input folder's manifest. Silent on error (do not block UI)."""
        p_input = parameters[0]
        p_variable = parameters[5]
        p_statistic = parameters[6]
        p_metric = parameters[7]

        if not p_input.value:
            return
        try:
            items, source, err = discover_t4_inputs(p_input.valueAsText)
            if err:
                return
        except Exception:
            return
        if not items:
            return

        variables = sorted({it["variable"] for it in items if it["variable"]})
        if p_variable.filter.list != variables:
            p_variable.filter.list = variables

        statistics = sorted({it["statistic"] for it in items if it["statistic"]})
        if p_statistic.filter.list != statistics:
            p_statistic.filter.list = statistics

        if source == "tool3":
            metrics = sorted({it["metric"] for it in items if it["metric"]})
            if p_metric.filter.list != metrics:
                p_metric.filter.list = metrics
        else:
            if p_metric.filter.list != []:
                p_metric.filter.list = []

    def updateMessages(self, parameters):
        pass

    def execute(self, parameters, messages):
        in_folder = parameters[0].valueAsText
        region_fc = parameters[1].valueAsText
        region_id_field = parameters[2].valueAsText
        region_label_field = (parameters[3].valueAsText
                              if parameters[3].value else None)
        out_csv = parameters[4].valueAsText
        variable_filter = parameters[5].valueAsText
        statistic_filter = parameters[6].valueAsText
        metric_filter = parameters[7].valueAsText
        year_min = parameters[8].value
        year_max = parameters[9].value
        zonal_stats_raw = parameters[10].valueAsText
        verbose = bool(parameters[11].value)

        def log(msg, level="info"):
            if level == "warn":
                messages.addWarningMessage(msg)
            elif level == "error":
                messages.addErrorMessage(msg)
            else:
                messages.addMessage(msg)

        # Reset caches at start of each run (Tool 1 / Tool 3 convention).
        _MASK_CACHE.clear()
        _REGION_LABEL_CACHE.clear()

        t_start = _time.time()

        items, source, err = discover_t4_inputs(in_folder)
        if err:
            log(err, "error")
            return
        if not items:
            log("No items in manifest.", "error")
            return
        log("Detected source: {}".format(source))
        log("Manifest items: {}".format(len(items)))

        def _split_multi(raw):
            if not raw:
                return None
            parts = [p.strip("'\" ") for p in raw.split(";") if p.strip()]
            return set(parts) if parts else None

        var_set = _split_multi(variable_filter)
        stat_set = _split_multi(statistic_filter)
        metric_set = _split_multi(metric_filter)
        zonal_stats = [s.strip("'\" ").lower()
                       for s in (zonal_stats_raw or "mean").split(";")
                       if s.strip()]
        if not zonal_stats:
            zonal_stats = ["mean"]
        valid_stats = {"mean", "min", "max", "std", "count"}
        invalid = set(zonal_stats) - valid_stats
        if invalid:
            log("Unknown zonal stats ignored: {}".format(sorted(invalid)), "warn")
            zonal_stats = [s for s in zonal_stats if s in valid_stats]
        if not zonal_stats:
            log("No valid zonal stats requested.", "error")
            return

        # Filter items
        filtered = []
        for it in items:
            if var_set is not None and it["variable"] not in var_set:
                continue
            if stat_set is not None and it["statistic"] not in stat_set:
                continue
            if metric_set is not None and source == "tool3" and it["metric"] not in metric_set:
                continue
            if year_min is not None and it["year"] is not None and it["year"] < int(year_min):
                continue
            if year_max is not None and it["year"] is not None and it["year"] > int(year_max):
                continue
            filtered.append(it)

        if not filtered:
            log("No items match the filters.", "error")
            return
        log("Items after filters: {}".format(len(filtered)))
        log("Zonal stats: {}".format(", ".join(zonal_stats)))

        # Read first TIFF for grid metadata; assume all share same grid
        # (Tool 1 strict mode + Tool 3 cross-pipeline check enforce this).
        sample_path = filtered[0]["raster"]
        gt, srs_wkt, nrows, ncols = read_grid_metadata(sample_path)
        lat, lon = derive_lat_lon_from_gt(gt, nrows, ncols)
        log("Grid: {} rows x {} cols".format(nrows, ncols))

        # Build region label raster (cached).
        labels = get_or_compute_region_labels(
            region_fc, region_id_field, lat, lon, lat_ascending=False
        )
        if labels is None:
            log("Failed to rasterise regions.", "error")
            return

        unique_ids = sorted(set(int(v) for v in np.unique(labels) if v >= 0))
        log("Regions resolved: {} (ids {} to {})".format(
            len(unique_ids),
            min(unique_ids) if unique_ids else "n/a",
            max(unique_ids) if unique_ids else "n/a"
        ))
        if not unique_ids:
            log("No polygon coverage on the data grid. Check that the "
                "region FC overlaps the raster extent.", "error")
            return

        # Build region_id -> label map (optional)
        id_to_label = {}
        if region_label_field:
            try:
                with arcpy.da.SearchCursor(
                    region_fc, [region_id_field, region_label_field]
                ) as cur:
                    for rid, lbl in cur:
                        if rid is None:
                            continue
                        id_to_label[int(rid)] = (lbl if lbl is not None else "")
            except Exception as e:
                log("Could not read region labels: {}".format(e), "warn")

        # Iterate items and accumulate rows.
        rows = []
        arcpy.SetProgressor("step", "Computing zonal stats",
                            0, len(filtered), 1)
        try:
            for idx, it in enumerate(filtered, 1):
                arcpy.SetProgressorLabel(
                    "Zonal {}/{}: {}".format(
                        idx, len(filtered), os.path.basename(it["raster"])
                    )
                )
                try:
                    data = read_tif_to_array(it["raster"])
                except Exception as e:
                    log("Read FAILED ({}): {}".format(it["raster"], e), "warn")
                    arcpy.SetProgressorPosition(idx)
                    continue
                if data.shape != labels.shape:
                    log("Grid mismatch on {} ({} vs {}). Skipped.".format(
                        it["raster"], data.shape, labels.shape), "warn")
                    arcpy.SetProgressorPosition(idx)
                    continue

                per_region = compute_zonal_stats(data, labels, zonal_stats)
                for rid, stats_dict in per_region.items():
                    row = {
                        "region_id": rid,
                        "region_label": id_to_label.get(rid, ""),
                        "variable": it["variable"],
                        "statistic": it["statistic"],
                        "metric": it["metric"],
                        "period_kind": it["period_kind"],
                        "year": it["year"] if it["year"] is not None else "",
                        "start_year": (it["start_year"]
                                       if it["start_year"] is not None else ""),
                        "end_year": (it["end_year"]
                                     if it["end_year"] is not None else ""),
                        "date": it["date"],
                        "units": it["units"],
                    }
                    for sk in zonal_stats:
                        row["zonal_" + sk] = stats_dict.get(sk, "")
                    rows.append(row)

                if verbose and idx % 50 == 0:
                    log("  processed {}/{}".format(idx, len(filtered)))
                arcpy.SetProgressorPosition(idx)
        finally:
            arcpy.ResetProgressor()

        if not rows:
            log("No rows generated (all reads failed or all regions empty).",
                "error")
            return

        # Write CSV
        out_dir = os.path.dirname(out_csv)
        if out_dir and not os.path.isdir(out_dir):
            os.makedirs(out_dir, exist_ok=True)
        base_cols = ["region_id", "region_label", "variable", "statistic",
                     "metric", "period_kind", "year", "start_year",
                     "end_year", "date", "units"]
        stat_cols = ["zonal_" + s for s in zonal_stats]
        fieldnames = base_cols + stat_cols
        try:
            with open(out_csv, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames,
                                        extrasaction="ignore")
                writer.writeheader()
                writer.writerows(rows)
        except OSError as e:
            log("Failed to write {}: {}".format(out_csv, e), "error")
            return

        elapsed = _time.time() - t_start
        log("Wrote {} rows to {} in {:.1f}s.".format(
            len(rows), out_csv, elapsed
        ))
