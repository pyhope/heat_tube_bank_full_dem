# Changelog

## 2026-06-23

This update refines the preliminary DEM postprocessing workflow.

- Changed the second panel of `flow_metrics` to use the outlet deletion counter
  from `global_metrics.dat`.  The processed mass flow is now calculated from
  cumulative deleted grains multiplied by the mean mobile-grain mass.
  The previous flow estimate is not reliable for this simulation.
- Updated contact-force handling in `contact_metrics`.  The force panel now
  uses the LAMMPS `compute pair/local` normal force directly.
- Removed the area-coverage panel from `grain_distributions`.  The current
  sampled area-coverage proxy depends strongly on contact-output frequency and
  does not track unique surface locations on rotating grains, so it should not
  be plotted as a final physical metric. Please find more details in `ANALYSIS_OUTPUTS.md`.
- Added `ANALYSIS_OUTPUTS.md`, which explains the CSV files, all figure panels,
  the main formulas used by `analyze_dem.py`, and the limitations of the
  current area-coverage proxy.
- Kept provisional coverage columns in `grain_summary.csv` for diagnostics, but
  documented that accurate area coverage will require high-frequency contact
  output, contact-episode tracking, and a surface-mixing model.
- Replaced the reference `results/` directory with a new postprocessed run using
  `nsteps = 1000000` and `dump_stop = 1000000`.  The main input file was also
  updated to match this longer 1M-step run.
- The current longer run still uses `dump_start = 200000`. Please increase it to
  at least 700000 as I mentioned in the email.

Note: the raw LAMMPS dumps are not stored in this GitHub repository because
they are large.  To regenerate the reference figures from scratch, run
`analyze_dem.py` in a complete LAMMPS run directory containing the raw dump
files.
