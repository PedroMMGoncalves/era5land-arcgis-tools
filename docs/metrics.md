# Climate indicators catalogue (Tool 3)

This document defines every metric implemented by Tool 3 of
`era5land-arcgis-tools`, with formulas, expected input, units, and the
peer-reviewed reference (where applicable). Use it as the authoritative
specification when interpreting outputs or writing methods sections.

Conventions used below:
- `s[t, y, x]` denotes the daily stack indexed by day `t`, latitude row
  `y`, longitude col `x`. `n` is the number of days in the year.
- Aggregations are NaN-aware (`nanmean`, `nansum`, etc.). Pixels with
  zero valid days return NaN to distinguish from "zero by computation".
- Auto-unit conversions are applied ONCE per stack before any metric:
  - Kelvin variables (`t2m`, `d2m`, `skt`, `stl1..4`, `tsn`) are
    converted from K to degC by subtracting 273.15.
  - Precipitation (`tp`) is converted from m to mm by multiplying by 1000.
  - Radiation variables (`ssrd`, `strd`, `ssr`, `str`) remain in raw
    J/m^2/day; each radiation metric applies its own per-formula factor.
  - Wind multi-input metrics build a daily magnitude stack
    `sqrt(u10^2 + v10^2)` in m/s before any metric runs.

Statistic-input expectations are warnings, not errors. The user can
proceed past a mismatch; results may be uninformative.

---

## Generic metrics (17)

Applicable to any variable. Output units inherit from working units unless
otherwise stated.

### mean
- Formula: `nanmean(s, axis=0)`.
- Units: working units of input.
- Reference: standard descriptive statistic.

### min
- Formula: `nanmin(s, axis=0)`.

### max
- Formula: `nanmax(s, axis=0)`.

### range
- Formula: `nanmax(s, axis=0) - nanmin(s, axis=0)`.

### std
- Formula: `nanstd(s, axis=0)` (population standard deviation, default
  NumPy ddof=0).

### sum
- Formula: `nansum(s, axis=0)`, NaN where all days invalid.

### median
- Formula: `nanmedian(s, axis=0)`.

### p10, p90
- Formula: `nanpercentile(s, q, axis=0)` for q = 10 and 90.
- Reference: percentile method = linear interpolation (NumPy default).

### count_above_threshold, count_below_threshold
- Formula: `sum((s > threshold) & valid_mask, axis=0)`, where
  `valid_mask = ~isnan(s)`. NaN-only pixels return NaN.
- Units: days.
- Parameter: `generic_threshold` (default 15). For Kelvin variables the
  threshold is interpreted in degC after auto-conversion; for other
  variables it is in raw input units.

### seconds_above_threshold, seconds_below_threshold
- Formula: `count_above_threshold * 86400`. Trivially redundant for
  daily input; kept for parity with hourly variants and unit-consumers
  that expect seconds.
- Units: s.

### accumulation_above_threshold, accumulation_below_threshold
- Formula: `nansum(where(s > threshold, s, 0), axis=0)` and analogous
  below. NaN-only pixels return NaN.
- Units: working units of input (interpret as units * day-equivalent
  depending on variable semantics; for precipitation in mm/day, sum
  reads as total mm above threshold).

### max_consecutive_above, max_consecutive_below
- Formula: longest run of consecutive days where `s > threshold` (or
  `s < threshold`). Implementation: running counter that resets when
  the condition fails; NaN evaluates to False and breaks runs.
- Units: days.

---

## Temperature-specific metrics (11)

Applicable to variable type "temperature" (`t2m`, `d2m`, `skt`, soil
temperatures, snow temperature). Stack is in degC by the time these
run.

### hdd_spinoni (Heating Degree Days, Spinoni method)
- Formula: `sum(max(hdd_base - T_mean, 0))` over days.
- Default base: 15.5 degC.
- Expected input statistic: `daily_mean`.
- Units: degC*day.
- Reference: Spinoni et al. (2018), "Changes of heating and cooling
  degree-days in Europe from 1981 to 2100", International Journal of
  Climatology 38: e191-e208. doi:10.1002/joc.5362.

### cdd_spinoni (Cooling Degree Days, Spinoni method)
- Formula: `sum(max(T_mean - cdd_base, 0))` over days.
- Default base: 22.0 degC.
- Expected input statistic: `daily_mean`.
- Units: degC*day.
- Reference: Spinoni et al. (2018) as above.

### gdd (Growing Degree Days)
- Formula: `sum(max(T_mean - gdd_base, 0))` over days.
- Default base: 10.0 degC.
- Expected input statistic: `daily_mean`.
- Units: degC*day.
- Reference: McMaster & Wilhelm (1997), "Growing degree-days: one
  equation, two interpretations", Agricultural and Forest Meteorology
  87: 291-300. doi:10.1016/S0168-1923(97)00027-0.

### frost_days
- Formula: count of days with `T_min < 0 degC` (ETCCDI FD).
- Expected input statistic: `daily_minimum`.
- Reference: Karl et al. (1999), "CLIVAR/GCOS/WMO workshop on indices
  and indicators for climate extremes", Climatic Change 42: 3-7;
  Zhang et al. (2011), "Indices for monitoring changes in extremes
  based on daily temperature and precipitation data", WIREs Climate
  Change 2: 851-870. doi:10.1002/wcc.147.

### ice_days
- Formula: count of days with `T_max < 0 degC` (ETCCDI ID).
- Expected input statistic: `daily_maximum`.
- Reference: Zhang et al. (2011), as above.

### summer_days
- Formula: count of days with `T_max > 25 degC` (ETCCDI SU).
- Expected input statistic: `daily_maximum`.
- Reference: Zhang et al. (2011), as above.

### tropical_nights
- Formula: count of days with `T_min > 20 degC` (ETCCDI TR).
- Expected input statistic: `daily_minimum`.
- Reference: Zhang et al. (2011), as above.

### hot_days
- Formula: count of days with `T_max > 30 degC`.
- Expected input statistic: `daily_maximum`.
- Reference: non-ETCCDI but commonly used in Iberian heat-wave
  literature; e.g. Pereira et al. (2017), "Heat wave and cold spell
  changes in Iberia for a future climate scenario", International
  Journal of Climatology 37: 5192-5205. doi:10.1002/joc.5158.

### heating_season_length (G.POT input)

- Formula: count of days with `T_mean < heating_threshold`.
- Expected input statistic: `daily_mean`.
- Parameter: `heating_threshold` (degC, default 15.5, matches Spinoni HDD
  base). User-tunable per regional convention.
- Units: days.
- Use case: direct input for the heating-season length `tc` parameter
  in shallow geothermal potential mapping. The G.POT method (eq. 14 in
  the reference) uses `tc` together with ground thermal properties and
  borehole geometry to compute the average thermal load that a Borehole
  Heat Exchanger can sustainably exchange.
- Distinct from `frost_days` / `ice_days` (fixed 0 degC threshold,
  Tmin/Tmax based) and from `count_below_threshold` Generic (any
  variable, not temperature-specific).
- Reference: Casasso, A. & Sethi, R. (2016), "G.POT: A quantitative
  method for the assessment and mapping of the shallow geothermal
  potential", Energy 106: 765-773. doi:10.1016/j.energy.2016.03.091.

### cooling_season_length (G.POT input)

- Formula: count of days with `T_mean > cooling_threshold`.
- Expected input statistic: `daily_mean`.
- Parameter: `cooling_threshold` (degC, default 22.0, matches Spinoni
  CDD base). User-tunable per regional convention.
- Units: days.
- Use case: direct input for the cooling-season length `tc` parameter
  in G.POT (cooling-mode operation).
- Distinct from `summer_days` (Tmax > 25 degC fixed), `hot_days`
  (Tmax > 30 degC fixed), `tropical_nights` (Tmin > 20 degC fixed),
  and `count_above_threshold` Generic (variable-agnostic).
- Reference: Casasso, A. & Sethi, R. (2016), as above.

### gsl (Growing Season Length, WMO Northern Hemisphere variant)
- Definition: number of days between the FIRST 6-day run of `T_mean
  > 5 degC` (season start) and the FIRST 6-day run of `T_mean < 5 degC`
  starting after Julian day 182 (season end). If no season start is
  found, GSL = 0. If a start exists but no end, season extends to
  day-of-year `n`.
- Implementation: cumsum-based sliding-window detection (window = 6).
- Expected input statistic: `daily_mean`.
- Units: days.
- Reference: Frich et al. (2002), "Observed coherent changes in
  climatic extremes during the second half of the twentieth century",
  Climate Research 19: 193-212. doi:10.3354/cr019193. ETCCDI GSL
  definition aligned with WMO (2017), "WMO Guidelines on the
  Calculation of Climate Normals" (WMO-No. 1203).

### txn (Annual minimum of daily-max temperature, ETCCDI TXn)
- Formula: `nanmin(s, axis=0)`.
- Expected input statistic: `daily_maximum`.
- Units: degC.
- Reference: Zhang et al. (2011), as above.

### tnx (Annual maximum of daily-min temperature, ETCCDI TNx)
- Formula: `nanmax(s, axis=0)`.
- Expected input statistic: `daily_minimum`.
- Units: degC.
- Reference: Zhang et al. (2011), as above.

---

## Precipitation-specific metrics (5)

Applicable to variable type "precipitation" (`tp`). Stack is in mm by
the time these run (auto-converted from m).

### total_precip
- Formula: `nansum(s, axis=0)`, NaN where all days invalid.
- Units: mm.
- Reference: standard climatological accumulation; ETCCDI PRCPTOT.

### wet_days
- Formula: count of days with precip >= threshold.
- Default threshold: 1.0 mm. ETCCDI R1mm convention; higher thresholds
  (10 mm, 20 mm) reproduce ETCCDI R10mm / R20mm.
- Units: days.
- Reference: Zhang et al. (2011), as above.

### consecutive_dry_days (ETCCDI CDD index)
- Formula: longest run of consecutive days with precip < threshold.
- Units: days.
- Reference: Frich et al. (2002), as above; Karl et al. (1999).

### consecutive_wet_days (ETCCDI CWD index)
- Formula: longest run of consecutive days with precip >= threshold.
- Units: days.
- Reference: Frich et al. (2002), as above.

### max_1day_precip (ETCCDI Rx1day)
- Formula: `nanmax(s, axis=0)`.
- Units: mm.
- Reference: Zhang et al. (2011), as above.

---

## Radiation-specific metrics (4)

Applicable to variable type "radiation" (`ssrd`, `strd`, `ssr`, `str`).
Stack remains in raw J/m^2/day; each metric applies its own conversion.

### CRITICAL: CDS daily statistic gotcha for accumulated variables

ERA5-Land radiation variables are delivered as joules per square metre
accumulated since 00:00 UTC of each day in the hourly product. To use
the formulas below correctly, the daily stack fed to Tool 3 MUST be
the daily total in J/m^2/day (i.e. the energy received during that day).

The CDS dataset `derived-era5-land-daily-statistics` (the one targeted
by the companion `era5land-downloads` repo for temperature) offers
`daily_mean`, `daily_maximum`, `daily_minimum` for these variables.
For accumulated variables, none of these is the daily total:
- `daily_mean` is the mean of the 24 hourly values after CDS
  decumulation, which yields a per-hour mean (factor ~24 smaller than
  the daily total).
- `daily_maximum` is the largest hourly increment (typically near solar
  noon), which is a power, not an energy.

The correct dataset for radiation totals is
`derived-era5-land-daily-aggregations` with statistic `daily_sum`. The
`expects_statistic` for radiation metrics is therefore declared as
`daily_sum`, and Tool 3 will emit a warning if the user provides a
different statistic.

If you must use `derived-era5-land-daily-statistics` daily_mean for
radiation:
- For `radiation_kwh_per_m2_year`: multiply the output by 24 manually
  (Raster Calculator), or rescale the inputs before Tool 1.
- For the other radiation metrics, factors differ; safer to switch
  dataset.

**Recommended workflow**: download radiation via
`derived-era5-land-daily-aggregations` (statistic `daily_sum`), name
the source NetCDF `era5land_<variable>_sum_<year>.nc` so Tool 1's
auto-detection sets `statistic="daily_sum"` in the manifest, then
Tool 3 runs cleanly.

### radiation_w_per_m2_mean
- Formula: `nanmean(s, axis=0) / 86400`.
- Rationale: ERA5-Land radiation variables are accumulated joules per
  square metre per day; dividing by 86400 s/day yields the daily-mean
  power in W/m^2. Annual mean of daily means.
- Expected input statistic: `daily_mean`.
- Units: W/m^2.
- Reference: ERA5-Land variable definition, Munoz-Sabater et al. (2021),
  "ERA5-Land: a state-of-the-art global reanalysis dataset for land
  applications", Earth System Science Data 13: 4349-4383.
  doi:10.5194/essd-13-4349-2021.

### radiation_kwh_per_m2_year
- Formula: `nansum(s, axis=0) / 3.6e6`.
- Rationale: 1 kWh = 3.6 * 10^6 J. Sum of daily J/m^2 over a year,
  divided by 3.6e6, yields annual energy in kWh/m^2.
- Expected input statistic: `daily_mean` (the daily total derived from
  ERA5-Land hourly accumulation; for CDS daily aggregations this is
  equivalent to `daily_sum` after unit clarification).
- Units: kWh/m^2.
- Reference: standard energy conversion; see IEC 61724-1 (2017),
  "Photovoltaic system performance, Part 1: Monitoring".

### photovoltaic_potential
- Formula: `(nansum(s, axis=0) / 3.6e6) * pv_performance_ratio`.
- Parameter: `pv_performance_ratio` (default 0.75). The PR captures
  inverter losses, module temperature, soiling, mismatch, and DC/AC
  conversion. Default 0.75 is a typical value for well-maintained
  residential PV (IEA-PVPS Task 13 reports 0.7-0.85 for European
  sites).
- Units: kWh/m^2 (per year of input).
- Reference: IEC 61724-1 (2017); Reich et al. (2012), "Performance
  Ratio Revisited: Is PR > 90% Realistic?", Progress in Photovoltaics
  20: 717-726. doi:10.1002/pip.1219.

### radiation_above_threshold_days
- Formula: count of days where `s > threshold_W_m2 * 86400` (threshold
  is provided in W/m^2 and converted internally to J/m^2/day).
- Default threshold: 100.0 W/m^2.
- Expected input statistic: `daily_mean`.
- Units: days.
- Reference: no standard ETCCDI equivalent; useful for solar resource
  availability windows.

---

## Wind-specific metrics (5, multi-input)

Require BOTH `u10` and `v10` at the chosen statistic. Tool 3
auto-detects availability and builds a daily wind-speed magnitude
stack `v = sqrt(u^2 + v^2)` in m/s before applying any metric.

### wind_speed_mean, wind_speed_max, wind_speed_p90
- Formula: `nanmean(v, axis=0)`, `nanmax(v, axis=0)`,
  `nanpercentile(v, 90, axis=0)` respectively.
- Units: m/s.
- Reference: Stull (1988), "An Introduction to Boundary Layer
  Meteorology", Kluwer.

### wind_power_density
- Formula: `0.5 * rho * nanmean(v^3, axis=0)`.
- Parameter: `wind_air_density` rho (default 1.225 kg/m^3, ICAO
  standard atmosphere at sea level, 15 degC).
- Units: W/m^2.
- Caveats:
  - This is the simple instantaneous-mean approximation. For bankable
    wind resource assessment, fit a Weibull distribution to the speed
    sample first. See Manwell et al. (2010).
  - Default rho assumes sea level. For high-altitude sites (Serra da
    Estrela approx 1900 m, Madeira interior, Pico in the Azores)
    actual rho is approx 10-15 percent lower, which translates linearly
    to a same-order overestimate of wind power density. Adjust the
    `wind_air_density` parameter accordingly: rho approx 1.225 *
    exp(-h / 8400) where h is elevation in metres (barometric formula,
    isothermal approximation).
- Reference: Manwell, McGowan & Rogers (2010), "Wind Energy Explained:
  Theory, Design and Application", 2nd ed., Wiley. IEC 61400-12-1
  (2017) for measurement standards.

### calm_days
- Formula: count of days with wind speed below threshold.
- Default threshold: 2.0 m/s (Beaufort 1, "light air").
- Units: days.
- Reference: Beaufort wind force scale; WMO (2017) Manual on Codes
  Volume I.1.

---

## DTR (Diurnal Temperature Range, 4 metrics, multi-input)

Requires BOTH `daily_maximum` and `daily_minimum` for the same
temperature variable. Tool 3 auto-detects availability and builds a
daily DTR stack `dtr = T_max - T_min` in degC before applying any
metric.

### dtr_mean, dtr_max, dtr_p90, dtr_std
- Formula: `nanmean(dtr, axis=0)`, `nanmax(dtr, axis=0)`,
  `nanpercentile(dtr, 90, axis=0)`, `nanstd(dtr, axis=0)`.
- Units: degC.
- Reference: ETCCDI DTR metric definition; Karl, Knight & Plummer
  (1995), "Trends in high-frequency climate variability in the
  twentieth century", Nature 377: 217-220. doi:10.1038/377217a0.

---

## Climatology and IPMA normals

- **Climatology**: when enabled, each annual metric is aggregated
  across the full year range available in the input. The default
  aggregator is mean (i.e. climatology of HDD = mean of annual HDD
  over the period). For metrics where the spec declares a different
  aggregator (currently none in v0.2.0), that is used instead.
  Output path: `climatology/<variable>/<statistic>/...climatology_<start>-<end>.tif`.
- **IPMA normals**: when enabled, the same aggregation is applied to
  the subset of years falling within 1981-2010 and 1991-2020. Warnings
  are emitted if the input does not cover the full period; the normal
  is still computed on the available subset (partial coverage is
  logged). Output path:
  `ipma_normals/<variable>/<statistic>/...normal_1981-2010.tif`
  and `..._normal_1991-2020.tif`.
- Reference: WMO (2017) defines climatological standard normals over
  30-year periods (1961-1990, 1991-2020). IPMA (Instituto Portugues
  do Mar e da Atmosfera) publishes Portuguese normals for 1971-2000
  and 1981-2010; the 1981-2010 period is the most-cited reference in
  recent Portuguese climate studies.

---

## References (canonical list)

1. Frich, P. et al. (2002). Observed coherent changes in climatic
   extremes during the second half of the twentieth century. Climate
   Research 19: 193-212. doi:10.3354/cr019193.
2. IEC 61400-12-1 (2017). Wind energy generation systems, Part 12-1:
   Power performance measurements of electricity producing wind
   turbines.
3. IEC 61724-1 (2017). Photovoltaic system performance, Part 1:
   Monitoring.
4. Karl, T.R., Knight, R.W. & Plummer, N. (1995). Trends in
   high-frequency climate variability in the twentieth century.
   Nature 377: 217-220. doi:10.1038/377217a0.
5. Karl, T.R., Nicholls, N. & Ghazi, A. (1999). CLIVAR/GCOS/WMO
   workshop on indices and indicators for climate extremes:
   workshop summary. Climatic Change 42: 3-7.
6. Manwell, J.F., McGowan, J.G. & Rogers, A.L. (2010). Wind Energy
   Explained: Theory, Design and Application, 2nd ed. Wiley.
7. McMaster, G.S. & Wilhelm, W.W. (1997). Growing degree-days: one
   equation, two interpretations. Agricultural and Forest Meteorology
   87: 291-300. doi:10.1016/S0168-1923(97)00027-0.
8. Munoz-Sabater, J. et al. (2021). ERA5-Land: a state-of-the-art
   global reanalysis dataset for land applications. Earth System
   Science Data 13: 4349-4383. doi:10.5194/essd-13-4349-2021.
9. Pereira, S.C. et al. (2017). Heat wave and cold spell changes in
   Iberia for a future climate scenario. International Journal of
   Climatology 37: 5192-5205. doi:10.1002/joc.5158.
10. Reich, N.H. et al. (2012). Performance Ratio Revisited: Is
    PR > 90% Realistic? Progress in Photovoltaics 20: 717-726.
    doi:10.1002/pip.1219.
11. Spinoni, J. et al. (2018). Changes of heating and cooling
    degree-days in Europe from 1981 to 2100. International Journal
    of Climatology 38: e191-e208. doi:10.1002/joc.5362.
12. Stull, R.B. (1988). An Introduction to Boundary Layer
    Meteorology. Kluwer Academic Publishers.
13. WMO (2017). WMO Guidelines on the Calculation of Climate Normals.
    WMO-No. 1203, World Meteorological Organization, Geneva.
14. Zhang, X. et al. (2011). Indices for monitoring changes in
    extremes based on daily temperature and precipitation data.
    WIREs Climate Change 2: 851-870. doi:10.1002/wcc.147.
