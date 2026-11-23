# bird_summary

Disk-first summaries of BAM bird predictions, following the single/multi pattern
of [Biomass_summary](https://github.com/PredictiveEcology/Biomass_summary/blob/modsForFireSense/Biomass_summary.R).
Exactly one BAM bootstrap is paired with each landscape replicate.

## Single mode: within a replicate run

```r
params <- list(
  bird_summary = list(
    mode = "single",
    years = c(2020L, 2120L),
    species = c("GCKI", "RCKI"),
    predictionDensityUnits = "males/ha",
    scenario = "NRV_ssp370",
    climateScenario = "CanESM5_SSP370",
    .studyAreaName = "ELF_6.2.2"
  )
)
```

`repID` and `bootstrapID` default to the predictor's `repID` and sole
`bootstrapIDs` value. Explicit values must agree with that pairing. With a
provider other than the predictor, specify both. `repID` uses `rep1`, not numeric
`1` or padded `rep01`. Different replicates may use the same bootstrap.

`years` is an exact selection, not a range. If omitted, single mode uses the
predictor calendar. If `species` is omitted, it uses the loader's model index
and species selection. This selection is resolved when saving, not during
`.inputObjects` or summary initialization: the loader creates `modelIndex` in
its `init` event. A missing index still stops the save rather than silently
omitting species. Saving occurs at each selected year at last event
priority, after predictions; no other module's events are scheduled.

Single mode reads the existing GeoTIFF referenced by
`sim$summaries[[species]][[year]]$meanPath` for every requested species. It
validates the single-bootstrap filename and layer, calculates metrics from the
saved raster, and writes/registers only a pair/year CSV with
`SpaDES.core::registerOutputs()`. It does not write or copy TIFFs. The predictor
already registers its own raster files. No permanent summary table is added to
the sim. An alternative provider must supply the same saved-file contract.

Outputs are under:

```text
outputPath/bird_summary/SCENARIO/CLIMATE/STUDY/rep1/year_2020/
  birdPairSummary_year2020_rep1_boot-31.csv
```

Each CSV is also a manifest: species, year, replicate, bootstrap, scenario,
climate, reporting-area label, units, metrics, and absolute `rasterPath`.
`sourceDensityUnits` records the raw raster units, independently of the reported
`densityUnits`, which is always males/ha. Abundance is reported as males.
`rasterPath` points to the predictor's existing file, which must remain
accessible and unchanged for multi mode. Predictor paths do not include the
summary's scenario/climate/area labels: use separate predictor output roots
for runs that would otherwise overwrite the same replicate/year/species/bootstrap.
Missing or mismatched files stop the save; there is no duplicate-write fallback.
Rerunning the same identity overwrites its CSV; different scenarios,
climates, areas, replicates and years have separate directories.
Older CSVs referencing summary-owned TIFF snapshots remain supported. Existing
snapshots are not deleted or rewritten automatically; retain any TIFFs still
referenced by CSVs.

## Multi mode: standalone postprocessing

```r
params <- list(
  bird_summary = list(
    mode = "multi",
    simOutputPath = "/path/to/simulation/outputs",
    reps = 1L:5L,
    years = c(2020L, 2120L),
    species = c("GCKI", "RCKI"),
    predictionDensityUnits = "males/ha",
    currentCondition = "both",
    bamRegion = "Canada",
    scenario = "NRV_ssp370",
    climateScenario = "CanESM5_SSP370",
    .studyAreaName = "ELF_6.2.2"
  )
)
```

Multi mode needs no dummy `summaries`, models, or saved simLists. BAM and both
download missing official products. Pass `outputsDF` as a combined `data.table::rbindlist()` of the
completed runs' `SpaDES.core::outputs(sim)` records if available. Registration
does **not** combine simulations automatically, so providing `objects =
list(outputsDF = outputsDF)` is optional. If omitted (or empty), the module
recursively scans `simOutputPath`. With a nonempty registry, it uses those
paths only: missing registered files cause an error, not a fallback to disk.
Scanning predictor rasters alone is insufficient; single-mode pair CSVs must
exist. Recorded paths need to remain accessible if outputs are moved.

Expected `reps`, `years`, and `species` must be explicit. The module selects one
scenario/climate/area context and rejects missing or duplicate species/year/rep
records, inconsistent bootstrap mappings, missing rasters, mismatched identities,
different grids or NA masks, and metrics inconsistent with their rasters.
It does not resample simulation rasters or silently omit incomplete runs.

For legacy pair CSVs without `sourceDensityUnits`, multi mode first verifies
their stored metrics under their recorded response-unit interpretation, then
regenerates metrics from their existing rasters using the explicitly configured
`predictionDensityUnits`, with a message. Source rasters and per-year CSVs are
not rewritten. New combined CSVs record the confirmed source units, preventing
double conversion. Previous correctly normalized males/km2 records are also
verified and converted to reported males/ha without changing abundance.
Conflicting recorded/configured raw units stop processing.

Under the context's `multi/` directory it writes/registers:

- `birdPairSummary_all.csv`: complete selected pair/year records and source paths.
- `birdAcrossPairSummary.csv`: annual across-pair mean, median, and pair counts
  for area-weighted mean density, pixel median density, and abundance, with
  an explicit units column.
- `birdAcrossReplicateDensity_SPECIES.png`: per-replicate trajectories of
  area-weighted mean density (males per hectare).
- `birdAcrossReplicateTotalAbundance_SPECIES.png`: total predicted males.
- `birdAcrossReplicateMedianDensity_SPECIES.png`: unweighted median pixel
  density (males per hectare).

The multi-mode `plotMetrics` parameter defaults to
`c("meanDensity", "totalAbundance", "medianDensity")`. Supply a subset to select
plots, or `character(0)` to disable them. All CSV metrics are still calculated.
Each plot title contains the species code and metric. Grey trajectories show
individual replicates; the blue trajectory shows their annual median for each
year across the selected replicates. Trajectories are shown without point
markers, and vertical dotted guides mark each plotted year.
With SIM or both selected and 2020 present, a red dashed horizontal line marks the median of the
selected replicates' 2020 values for that species and metric. This is the
model-predicted current condition, not observed bird counts. The individual
2020 values remain in `birdPairSummary_all.csv`; their median is recorded in
`birdAcrossPairSummary.csv` for year 2020. When 2020 is not selected, the module
reports that the baseline line is omitted; it does not substitute another year.
Y-axis ticks and horizontal gridlines use matching, equally spaced intervals.
Each plot retains its own scale, rather than sharing limits across species or metrics.
Legends occupy a separate panel below the axes, outside the curve area.
PNG dimensions are 1200 by 950 pixels.
All three plot types shade the visible portion of simulation years 2720-3020
in pale beige, labelled "FireSense NRV reference years" in the legend. Plots
outside that period have no shading or reference legend entry. This identifies
the FireSense reference period visually only: it does not filter records,
change summary statistics, or calculate NRV bounds.
Only plots from the current selection are registered; old plots from earlier
selections are not deleted.

The annual summaries are descriptive response metrics, not NRV bounds. Mean density is useful
for density comparisons, and abundance for reporting total predicted males
within the reporting area. With identical reporting area/support, total abundance
and mean density are proportional, so they are not independent indices.
Median pixel density describes the typical pixel but does not integrate to total
abundance.

## Optional FireSense-style NRV sample bounds

Enable this only in standalone multi mode, for a scenario you intend to use
as a natural reference:

```r
bird_summary = list(
  mode = "multi",
  # Retain your paths, species, scenario, climateScenario and .studyAreaName.
  reps = 1L:5L,
  years = seq(2020L, 3020L, 100L),
  calculateNRV = TRUE
)
```

The default `calculateNRV = FALSE` retains descriptive summaries. When enabled,
all four exact years **2720, 2820, 2920, 3020** must be included in `years` and
have one record for every selected species and replicate. Other years remain
in the trajectories and annual CSVs but do not contribute to these bounds.
Missing snapshots stop the run instead of yielding a partial reference sample.

For each species and each of mean density, median pixel density and total
abundance, calculate the minimum and maximum across the individual
replicate/year values, **without averaging replicates first**. Five selected
replicates give 20 observations per species/metric. Fewer or more replicates
are supported: the sample always contains four observations per selected
replicate.

The module writes/registers `birdNRVBounds.csv` even when plots are disabled.
It contains `lower`, `upper`, `nObservations`, `nPairs`, `nYears`,
`referenceYears`, `method = "sample_min_max"`, units, species and the selected
scenario/climate/area labels. The existing combined pair CSV retains the
underlying records and their raster/source paths. Selected plots add dashed
bound lines across the trajectory and report replicate/observation counts;
the pale reference-year band remains. Annual statistics and rasters are unchanged.

This is the user-selected FireSense-style **sample envelope**, not a verified
reimplementation of the FireSense report, a confidence interval, or guaranteed
ecological limits. Extreme observations and sample size affect min-max bounds.
Snapshots may be temporally correlated, and paired BAM bootstraps combine
landscape variation with model uncertainty. Climate trends may shift the
reference distribution. Enabling this option does not establish that the selected
scenario represents natural conditions. No current-status classification or
spatial NRV map is calculated. Old bounds files are not deleted when the option
is later disabled; only products from the current run are registered.

## Scientific meaning and limits

Unmodified BAM V5 models return **males/ha**, even on a 1-km predictor grid.
`predictionDensityUnits` therefore defaults to `"males/ha"`; use `"males/km2"`
only for raw rasters already converted to that scale. All reported densities
use **males/ha**, matching BAM. Source rasters remain unchanged. Grid resolution
does not define the density denominator.

Versions before 0.1.12 incorrectly interpreted raw V5 values as males/km2.
Regeneration fixes those density labels and increases their abundance by 100.
Version 0.1.12 instead reported correctly converted males/km2 densities:
those density numbers divide by 100 when converted to males/ha, but their
already-correct abundance is unchanged. Neither migration changes predictions.

This follows the response-scale description and prediction call in
[V5's prediction script](https://github.com/borealbirds/LandbirdModelsV5/blob/8190ec97c6622bf5295f08aa6d0517acfd9fc729/analysis/07.Predict/07.Predict.R).
Cell areas are computed with
`terra::cellSize(unit = "ha")`, not assumed to be 100 ha per cell. The implementation
equivalently uses `cellSize(unit = "km") * 100`, retaining `validAreaKm2` as
an area field. Mean density is area
weighted; median density is the unweighted median of valid pixel values.
Total abundance is the sum of density times cell area over nonmissing support.
Valid area and cell counts are recorded. No-data areas are not counted as zero.
All-NA, infinite or negative predictions fail. Common raster support is required
within each species across years and pairs.

Abundance is `sum(densityMalesPerHa * cellAreaHa)` on the reporting mask.

Across-pair statistics describe the paired landscape/model predictions, not
separate landscape and model uncertainties. Years remain separate in annual
summaries; the optional sample envelope pools only the four specified reference
snapshots. Metadata alone does not establish scientific validity.
Temporal changes, spatial summary maps, and numerical current-status-versus-NRV
classifications are outside this version.
The year 2020 is not treated as an NRV reference by default.

## Current condition and official BAM 2020 benchmark

In multi mode, set `currentCondition` to `"BAM"` (the default), `"SIM"`, or
`"both"`. BAM compares the official 2020 estimate with simulated NRV; SIM
uses the simulated 2020 replicate median; both shows both. SIM and both need
2020 in `years` for a current-condition comparison; both requires it.
All modes still need the simulated NRV records when `calculateNRV = TRUE`.

BAM or both automatically downloads missing V5 2020 species products on first
use, using BAMexploreR, and reuses local products on later runs. Products are
stored in `inputPath/bird_summary/BAM_v5_2020/REGION/`. The project setup
installs the declared BAMexploreR dependency if needed; downloads need network
access. Whole-region products can be large. Set `bamRegion` to `"Canada"`
(default), `"Alaska"`, or `"Lower48"` to select the product containing the
study area. Products are not combined for cross-region areas. Missing species,
invalid files and download failures are errors; there is no silent SIM fallback.
No downloads occur in single mode or SIM mode.

Alternatively, supply `objects = list(bamBenchmarkFiles = manifest)` in multi
mode. A supplied manifest takes precedence and must cover every selected species.
It has one selected row per species:

```r
manifest <- data.table::data.table(
  species = "GCKI",
  year = 2020L,
  rasterPath = "/path/to/GCKI_Canada_2020.tif",
  densityUnits = "males/ha",
  modelVersion = "v5",
  sourceURL = "https://doi.org/10.5281/zenodo.21632211",
  nBootstraps = 32L,
  ensembleStatistic = "mean"
)
```

Only the named `mean` band is used, not `standard_deviation`, `cv`, or
`detection_distance`. Automatic acquisition uses the whole-region product,
then `reproducible::prepInputs` loads the persistent local file. Existing products
are reused without downloading or overwriting; invalid files fail explicitly.
The persistent input file is the download cache: caching a file-writing download
as a function result would not restore a deleted file.
SIM mode uses the full simulation reporting area, without a BAM mask.

The official mean density is bilinearly resampled/projected to each species'
simulation grid if necessary, with an explicit message. Separately, its native
valid-data footprint (`!is.na(mean)`) is aligned by nearest neighbour. This
prevents bilinear interpolation from filling gaps in that footprint. Valid
zero predictions remain valid cells; NA cells are excluded, never changed to zero.

For each species, the reporting mask is the intersection of this BAM-2020
footprint and the existing simulation reporting area. **The same fixed mask**
is applied to official BAM 2020, every simulated 2020 replicate, and every
selected NRV/trajectory year and replicate. Metrics are recalculated from
the existing simulation rasters, not scaled from earlier CSV summaries.
Missing simulation values inside the fixed mask stop comparison rather than
shrinking it differently for individual years. An empty overlap also stops.

This is a reporting-only operation. Neither source rasters nor the underlying
LandR/fireSense simulations are modified. The resulting means and abundance
describe the **BAM-supported reporting area**, not the whole ELF or necessarily
all land. BAM also applies data-availability masks. Masks can differ among
species, so use the recorded areas when comparing results across species.
Bilinear interpolation treats density as continuous; it does not guarantee
preservation of the native raster's integrated total. Reported benchmark
abundance is integrated on the aligned simulation reporting grid.

Benchmark-enabled products are isolated under `multi/bam_supported/`, leaving
full-support products in `multi/` untouched. In addition to the usual pair,
annual and optional NRV/plot products, the module saves/registers:

- `birdBAMReportingMask_SPECIES.tif`: 1 on retained cells, NA elsewhere.
- `birdBAMBenchmark2020.csv`, including per-hectare
  metrics, raw units, raster path, source URL, model version, bootstrap count,
  alignment method, original/excluded cell counts and reporting area.

Pair, benchmark, annual and NRV tables record
`reportingSupport = "bam_2020_valid"` and `reportingMaskPath`; retained cell counts
and areas accompany the derived statistics. The pair table's `rasterPath` still points to the original
full raster; apply its saved mask to reproduce the reported metrics. Single-mode
records instead use `reportingSupport = "full_prediction"`.

Selected trajectories show the requested current-condition line or lines.
The benchmark is a metric of the **bootstrap-mean
surface**, not a median of bootstrap-level landscape metrics.

With `calculateNRV = TRUE`, the module also plots the selected current condition
against NRV. SIM requires 2020 in `years` to produce the comparison plot.
The figure shows simulated 2020 as a replicate median and range, BAM 2020 as
one estimate without an uncertainty whisker, and NRV as a median and
2.5th-97.5th percentile interval across the reference replicate/year values.
These are descriptive variability intervals, not confidence intervals. The
separate NRV bounds CSV and trajectory lines remain sample min-max. Pixel-level
SD/CV cannot provide a landscape interval without spatial covariance or
bootstrap-level surfaces. The module does not classify ecological status.

When `calculateNRV = TRUE` and `plotMetrics` includes `meanDensity`, the module
also registers multi-species `birdCurrentNRVWhiskers_MODE_pageNN.png` figures.
They show the NRV median and observed sample min-max, together with the
selected current condition. Select `currentCondition = "both"` to show both
BAM and the simulated 2020 median and replicate range. They display landscape
mean density in males/km2; existing density plots and CSVs remain in males/ha.
Pages contain up to seven species. Sample ranges are not confidence intervals.
SIM requires 2020 in `years` to produce the comparison figures.

The workspace runner exposes `currentCondition` and `bamRegion` in `birdRun`.

Official V5 products undergo truncation, range limitation and data/water masking
before packaging; see
[truncation](https://github.com/borealbirds/LandbirdModelsV5/blob/8190ec97c6622bf5295f08aa6d0517acfd9fc729/analysis/10.Truncate.R)
and [packaging](https://github.com/borealbirds/LandbirdModelsV5/blob/8190ec97c6622bf5295f08aa6d0517acfd9fc729/analysis/11.Package.R).
These steps and different baseline covariates can explain differences from raw
simulation predictions even with matching model versions. Matching units and
coverage alone does not make the processing scientifically identical.

Run local fixture tests from this module directory:
`testthat::test_dir("tests/testthat")`.

The scaffold's license placeholder has not yet been replaced; choose a license
before publishing the module.
