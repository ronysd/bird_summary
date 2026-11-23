# bird_summary 0.1.15

- Add multi-species median and sample min-max whisker plots in males/km2
  alongside existing NRV comparison and trajectory plots.
- Preserve NRV bounds when simulated 2020 is not selected.

# bird_summary 0.1.14

- Multi-mode current condition can use BAM (default), SIM, or both. Missing
  official BAM products are downloaded once to the configured input directory;
  supplied manifests remain supported.
- Current-versus-NRV plots show replicate ranges and descriptive NRV percentile
  intervals rather than labelling them as confidence intervals.

# bird_summary 0.1.13 (19 November 2026)

- Report all densities in males/ha, matching BAM. Integrate with cell areas in
  hectares; abundance remains predicted males. Verify and migrate both legacy
  records and previously normalized males/km2 records without double scaling.
- Only when an official BAM benchmark is supplied, use one fixed per-species
  BAM-2020 valid-data mask for the benchmark and all simulated years/replicates.
  Keep valid zeros. Align validity by nearest neighbour separately from
  bilinear density interpolation; fail on missing values within the fixed mask.
- Recalculate masked metrics from existing rasters. Save mask rasters and
  provenance under multi/bam_supported, separate from full-support summaries.
  Source rasters and simulations are unchanged. Add reporting regression
  coverage and documentation (MP-M4/MP-L1 reporting follow-up).

# bird_summary 0.1.12 (19 November 2026)

- Correct the unmodified V5 raw response interpretation to males/ha; normalize
  reported densities to males/km2 before area integration. Add explicit
  predictionDensityUnits and record sourceDensityUnits. Regenerate verified
  legacy metrics from existing rasters without changing rasters or original
  pair CSVs; V5 density/abundance values increase by 100.
- Add optional local official-BAM-2020 benchmark manifests, normalized metrics
  and a distinct trajectory baseline. Log bilinear alignment and require full
  existing simulation reporting coverage; no automatic acquisition or masking
  changes to simulations (MP-M4/MP-L1 reporting follow-up).
- Add simulated-2020/BAM-2020/NRV comparison plots when complete reference years
  and NRV calculation are enabled. Label sample ranges and the official
  bootstrap-mean surface distinctly; do not infer landscape confidence limits
  from pixel uncertainty. Record public-product postprocessing in provenance.

# bird_summary 0.1.11 (18 November 2026)

- Add a red dashed current-condition line at each species/metric's median
  across selected replicates at 2020, when that year is selected. Existing CSVs
  retain individual baseline values and their median; predictions are unchanged.
- Use equally spaced Y-axis ticks with matching horizontal gridlines and
  species/metric-specific scales (MP-M4/MP-L1 reporting follow-up).

# bird_summary 0.1.10 (18 November 2026)

- Remove point markers from replicate and median trajectories and add vertical
  dotted guides for plotted years.

# bird_summary 0.1.9 (18 November 2026)

- Plot the across-replicate median trajectory instead of the mean. Individual
  replicate lines and saved CSV statistics are unchanged.

# bird_summary 0.1.8 (9 November 2026)

- Move metric-plot legends into a separate panel below the axes so text cannot
  cover replicate curves or NRV bound lines. Increase PNG height to retain
  space for the data panel; calculations and filenames are unchanged.
- Test nonoverlapping data/legend panels with and without NRV bounds
  (MP-M4/MP-L1 reporting follow-up).

# bird_summary 0.1.7 (6 November 2026)

- Add opt-in `calculateNRV` in multi mode: sample minimum and maximum for each
  species and metric across individual snapshots at 2720, 2820, 2920 and 3020.
  Use all selected replicates, require complete snapshots, and report counts.
- Save/register `birdNRVBounds.csv` and show bound lines and sample counts on
  selected plots. Default descriptive mode and annual CSV metrics are unchanged.
- Follow the user-selected FireSense-style sample envelope, not confidence
  limits or percentile bounds; add focused regression coverage and documentation
  (MP-M4/MP-L1 reporting follow-up, no predictor changes).

# bird_summary 0.1.6 (6 November 2026)

- Shade 2720-3020 in all multi-mode metric plots and label the visible band
  "FireSense NRV reference years". CSVs and statistics are unchanged; this is
  a visual reference-period annotation, not an NRV range calculation.
- Add plot regression coverage for the annotation (MP-M4/MP-L1 reporting
  follow-up; no predictor changes).

# bird_summary 0.1.5 (5 November 2026)

- Single mode now references the predictor's existing `summaries` meanPath TIFFs
  instead of writing duplicate rasters. Metrics are calculated from those saved
  files; only pair/year CSVs are written and registered by single mode.
- Preserve multi-mode support for older snapshot-backed CSVs. Predictor files
  must remain accessible and unchanged; missing or mismatched files stop clearly.
- Document shared-file lifetime and predictor output collisions (MP-H3, MP-M2,
  MP-M4, MP-L1 shared-interface follow-up; no predictor code changes).

# bird_summary 0.1.4 (4 November 2026)

- Add `plotMetrics` to select mean density, total abundance, and/or median pixel
  density plots; default to all three, or disable with `character(0)`.
- Display density units with a superscript square and retain species codes in
  every plot title. CSV units, calculations, and existing density filenames
  remain unchanged. NRV reference distributions remain deferred.

# bird_summary 0.1.3 (4 November 2026)

- Fix species filtering in density plots: each figure now contains only its
  species rather than both species combined. Saved summary CSVs were unaffected.
- Check actual plotted replicate and mean values with distinct species fixtures.

# bird_summary 0.1.2 (4 November 2026)

- Fix automatic species selection failing during `.inputObjects`, before the
  model loader has created its index. Resolve inferred species at save time
  while retaining early parameter checks and strict save-time validation
  (MP-H3 shared-interface follow-up).

# bird_summary 0.1.1 (4 November 2026)

- Add per-species plots of replicate density trajectories and their annual mean.
- Clarify that multi mode can discover pair tables by scanning `simOutputPath`
  when no combined `outputsDF` is supplied.
- Keep NRV reference years and spatial comparison deferred.

# bird_summary 0.1.0 (3 November 2026)

- Add single-mode pair/year mosaics, area-weighted density, pixel median density
  and abundance CSVs with explicit identity, units, raster paths and registration.
- Add standalone multi-mode discovery from combined outputsDF or disk, complete
  selection validation, and annual across-pair mean/median products.
- Follow Biomass_summary naming and event conventions without loading saved sims
  or requiring dummy single-mode objects in multi mode.
- Defer reference NRV calculations until their sampling policy is explicitly
  defined; no SD, quantiles or NRV spatial summary maps are produced.
- Add deterministic unit and disk/registry integration fixtures (MP-H3, MP-M4,
  MP-L1 shared-interface follow-up).
