# Postprocessing Summary

## Run Status

- Slurm state: `COMPLETED`, exit code `0:0`
- Resources: 16 MPI tasks, 1 OpenMP thread per task
- Wall time: 46 min 52 s
- LAMMPS production segment: 400,000 steps from 40 to 120 ms
- Final reported mobile grains: 121,275
- Final grains in tube-bank region: 97,731
- Final deleted outlet grains: 45,608

## Main Postprocessed Outputs

- `summary.csv`: one-line aggregate summary
- `global_metrics_clean.csv`: cleaned global time series from `global_metrics.dat`
- `bed_time_series.csv`: particle number, mass, energy, and bed-position time series from the synchronized dump
- `contact_time_series.csv`: contact-count, force, and inferred disturbed-depth time series
- `grain_summary.csv`: per-grain residence time, contact counts, coverage proxies, and collision-energy proxies
- `figures/flow_metrics.pdf` and `.png`
- `figures/contact_metrics.pdf` and `.png`
- `figures/grain_distributions.pdf` and `.png`
- `figures/snapshot_dactual.pdf` and `.png`

## Key Numbers

- Mobile grains seen in synchronized dumps: 166,874
- Contact records analyzed: 7,515
- Velocity-resolved contact records: 7,108
- Grain-grain contacts: 7,101
- Grain-tube contacts: 414
- Median force-based `d_actual` proxy at `chi = 0.2`: 0.0249 um
- 90th percentile force-based `d_actual` proxy at `chi = 0.2`: 0.0365 um
- Median velocity-based `d_actual` proxy at `chi = 0.2`: 0.0481 um
- 90th percentile velocity-based `d_actual` proxy at `chi = 0.2`: 0.0763 um
- Median normal contact force: 5.24e-5 N
- 90th percentile normal contact force: 1.65e-4 N
- Median normal relative velocity: 0.0560 m/s
- 90th percentile normal relative velocity: 0.179 m/s
- Median tube-bank residence time in sampled frames: 0.044 s
- 90th percentile tube-bank residence time in sampled frames: 0.076 s
- Tangential work proxy per mobile mass seen: 0.300 J/kg
- Normal collision-energy proxy per mobile mass seen: 1.23e-4 J/kg

## Steady-Flow Check

The tube-bank population is much flatter than the total mobile-particle population near the end of the sampled window:

- Mean `nbank` over the last 10 global samples: 98,035 +/- 293 grains
- Mean `nmobile` over the last 10 global samples: 118,194 +/- 3,160 grains
- Linear trend over the second half of the run:
  - `nmobile`: +1.14e4 grains per 100,000 steps
  - `nbank`: +8.33e3 grains per 100,000 steps
  - outlet deletions: +1.50e4 grains per 100,000 steps

This test is therefore useful as a framework and order-of-magnitude check, but it should not be treated as a converged steady-state production result. Matt should test a later statistics window, longer run time, and sensitivity to feed rate, timestep, particle number, and box size before using time-averaged parameters quantitatively.
