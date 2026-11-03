# bird_summary 0.1.0 (3 November 2026)

- Add single-mode pair/year mosaics, area-weighted density, pixel median density
  and abundance CSVs with explicit identity, units, raster paths and registration.
- Add standalone multi-mode discovery from combined outputsDF or disk, complete
  selection validation, and annual across-pair mean/median products.
- Follow Biomass_summary naming and event conventions without loading saved sims
  or requiring dummy single-mode objects in multi mode.
- Defer reference NRV calculations until their sampling policy is explicitly
  defined; no SD, quantiles or spatial summary maps are produced.
- Add deterministic unit and disk/registry integration fixtures (MP-H3, MP-M4,
  MP-L1 shared-interface follow-up).
