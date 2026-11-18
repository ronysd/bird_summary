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
    simOutputPath = "/home/sdas/SpaDES_book/LandRDemo_coreVeg/outputs",
    reps = 1L:5L,
    years = c(2020L, 2120L),
    species = c("GCKI", "RCKI"),
    scenario = "NRV_ssp370",
    climateScenario = "CanESM5_SSP370",
    .studyAreaName = "ELF_6.2.2"
  )
)
```

Multi mode needs no dummy `summaries`, models, downloads, or saved
simLists. Pass `outputsDF` as a combined `data.table::rbindlist()` of the
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
It does not resample data or silently omit incomplete runs.

Under the context's `multi/` directory it writes/registers:

- `birdPairSummary_all.csv`: complete selected pair/year records and source paths.
- `birdAcrossPairSummary.csv`: annual across-pair mean, median, and pair counts
  for area-weighted mean density, pixel median density, and abundance, with
  an explicit units column.
- `birdAcrossReplicateDensity_SPECIES.png`: per-replicate trajectories of
  area-weighted mean density (males per square kilometre).
- `birdAcrossReplicateTotalAbundance_SPECIES.png`: total predicted males.
- `birdAcrossReplicateMedianDensity_SPECIES.png`: unweighted median pixel
  density (males per square kilometre).

The multi-mode `plotMetrics` parameter defaults to
`c("meanDensity", "totalAbundance", "medianDensity")`. Supply a subset to select
plots, or `character(0)` to disable them. All CSV metrics are still calculated.
Each plot title contains the species code and metric. Grey trajectories show
individual replicates; the blue trajectory shows their annual median for each
year across the selected replicates. Trajectories are shown without point
markers, and vertical dotted guides mark each plotted year.
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

The confirmed model response is **males/km2**. Cell areas are computed with
`terra::cellSize(unit = "km")`, not assumed to be 1 km2. Mean density is area
weighted; median density is the unweighted median of valid pixel values.
Total abundance is the sum of density times cell area over nonmissing support.
Valid area and cell counts are recorded. No-data areas are not counted as zero.
All-NA, infinite or negative predictions fail. Common raster support is required
within each species across years and pairs.

A future males/ha conversion divides density by 100; integrating in hectares
must preserve total predicted males. No unit-conversion parameter is included.

Across-pair statistics describe the paired landscape/model predictions, not
separate landscape and model uncertainties. Years remain separate in annual
summaries; the optional sample envelope pools only the four specified reference
snapshots. Metadata alone does not establish scientific validity.
SD, quantiles, temporal changes, spatial
summary maps, and current-status-versus-NRV comparisons are outside this version.
The year 2020 is not treated as an NRV reference by default.

Run local fixture tests from this module directory:
`testthat::test_dir("tests/testthat")`.

The scaffold's license placeholder has not yet been replaced; choose a license
before publishing the module.
