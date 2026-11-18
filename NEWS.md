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
