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
and species selection. Saving occurs at each selected year at last event
priority, after predictions; no other module's events are scheduled.

Single mode validates one-layer `bBOOT_YEAR` mosaics for every requested species,
writes a durable copy and a pair/year CSV, and registers both with
`SpaDES.core::registerOutputs()`. No permanent summary table is added to the sim.

Outputs are under:

```text
outputPath/bird_summary/SCENARIO/CLIMATE/STUDY/rep1/year_2020/
  2020_GCKI_mosaic_rep1_boot-31.tif
  birdPairSummary_year2020_rep1_boot-31.csv
```

Each CSV is also a manifest: species, year, replicate, bootstrap, scenario,
climate, reporting-area label, units, metrics, and absolute `rasterPath`.
Rerunning the same identity overwrites its products; different scenarios,
climates, areas, replicates and years have separate directories.

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

Multi mode needs no dummy `mosaicedPredictions`, models, downloads, or saved
simLists. Pass `outputsDF` as a combined `data.table::rbindlist()` of the
completed runs' `SpaDES.core::outputs(sim)` records if available. Registration
does **not** combine simulations automatically. Without a nonempty registry,
the module recursively scans `simOutputPath`. With a registry, it uses those
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
separate landscape and model uncertainties. Years remain separate. **NRV is not
calculated**: reference scenario, period, burn-in, temporal weighting and bounds
must be explicitly agreed before adding reference distributions. Metadata alone
does not establish an NRV policy. SD, quantiles, temporal changes, and spatial
summary maps are outside this initial version.

Run local fixture tests from this module directory:
`testthat::test_dir("tests/testthat")`.

The scaffold's license placeholder has not yet been replaced; choose a license
before publishing the module.
