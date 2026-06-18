# Continuous-Feed Tube-Bank DEM Model for Helium-3 Agitation

This directory contains a LAMMPS/DEM starting point for Matt's agitation branch
of the lunar helium-3 extraction project.  The model is inspired by the Helium
Extraction and Acquisition Testbed (HEAT) discussed by Olson, but it is not a
full thermal or engineering reproduction of that apparatus.  It is a mechanical
surrogate for regolith-like grains flowing through a heat-pipe/tube-bank region.

The current input is a baseline model, not a final production run.  Matt should
use it to learn the workflow, inspect the output, and perform convergence tests.
The final production parameters should be chosen only after those tests.

The model is designed to estimate DEM-derived quantities needed by the
reduced-order helium release model:

- tube-bank residence time,
- processed mass rate,
- grain-grain and grain-tube contact statistics,
- contact-force and impact-velocity distributions,
- surface-contact coverage,
- mechanically affected depth proxies, denoted `d_actual`,
- preliminary agitation-energy proxies per processed mass.

The current model uses continuous feeding at the top and outlet deletion at the
bottom. The goal is a statistically steady moving bed: grains enter, interact
with the tube bank, leave through the outlet, and are removed from the
simulation before they become lost atoms.

## Files

- `make_tube_particles.py`: generates fixed tube-surface beads.
- `tube_particles.in`: generated LAMMPS commands for fixed type-2 tube beads.
- `in.heat_tube_bank`: top-level LAMMPS input file.
- `in.heat_tube_bank_body`: LAMMPS body included by
  `in.heat_tube_bank`.
- `run_adroit.sbatch`: Adroit Slurm script for the baseline run.
- `analyze_dem.py`: postprocesses LAMMPS outputs and generates CSV files and
  figures.
- `requirements.txt`: Python packages required by `analyze_dem.py`.
- `results/`: baseline run results completed on Adroit. Some large output files (e.g., *.dump) are not included.

The production values for `nsteps`, `npart`, dump intervals, timestep, geometry, and
material/contact parameters should be decided after convergence tests.

## Geometry

The simulation domain is a narrow non-periodic channel:

| Quantity | Value |
| --- | ---: |
| `lx` | `1.2e-2 m` |
| `ly` | `2.4e-3 m` |
| `lz` | `3.2e-2 m` |
| Tube radius in `make_tube_particles.py` | `4.5e-4 m` |
| Tube-surface bead spacing | about `1.05e-4 m` |
| Number of tube centers | `21` |
| Generated fixed tube beads | `12474` |

The tube centers are defined in `make_tube_particles.py`:

| Row | `z` center (m) | `x` centers (m) |
| ---: | ---: | --- |
| 1 | `6.0e-3` | `2.0e-3`, `5.0e-3`, `8.0e-3`, `1.1e-2` |
| 2 | `1.0e-2` | `3.5e-3`, `6.5e-3`, `9.5e-3` |
| 3 | `1.4e-2` | `2.0e-3`, `5.0e-3`, `8.0e-3`, `1.1e-2` |
| 4 | `1.8e-2` | `3.5e-3`, `6.5e-3`, `9.5e-3` |
| 5 | `2.2e-2` | `2.0e-3`, `5.0e-3`, `8.0e-3`, `1.1e-2` |
| 6 | `2.6e-2` | `3.5e-3`, `6.5e-3`, `9.5e-3` |

Mobile regolith grains are type-1 spheres with diameters sampled from
`8.0e-5` to `1.2e-4 m`.  Fixed tube-surface beads are type-2 spheres.  The tube
bank analysis region is `bank_zlo = 5.2e-3 m` to `bank_zhi = 2.72e-2 m`.

## Boundary Conditions

The LAMMPS boundary is

```lammps
boundary f f f
```

The box is non-periodic in all directions.  Physical side walls are imposed
explicitly by

```lammps
fix wallx mobile wall/gran ... xplane 0.0 ${lx}
fix wally mobile wall/gran ... yplane 0.0 ${ly}
```

Particles that touch the `yz` planes at `x = 0` or `x = lx` interact with
`wallx`.  Particles that touch the `xz` planes at `y = 0` or `y = ly` interact
with `wally`.  These side-wall contacts use the same granular Hooke/history
parameters as other contacts.

The flow direction is `z`.  Grains are inserted near the top and removed in a
bottom outlet region:

```lammps
region outlet block 0 ${lx} 0 ${ly} 0 ${outlet_zhi} units box
fix outlet mobile evaporate ${outlet_every} ${outlet_max} outlet ${outlet_seed}
```

The outlet deletion is a controlled numerical sink.  It records the cumulative
number of deleted grains through `f_outlet`, which the postprocessor labels as
`outlet_deleted`.  The outlet is below the active tube-bank region, so it should
not be interpreted as a contact boundary inside the analysis zone.

## Baseline Parameters

`in.heat_tube_bank` defines the baseline parameter set and then includes
`in.heat_tube_bank_body`.

| Parameter | Value | Meaning |
| --- | ---: | --- |
| `run_label` | `baseline` | Label printed in `log.lammps`. |
| `lx` | `1.2e-2` | Domain length in `x`. |
| `ly` | `2.4e-3` | Domain width in `y`. |
| `lz` | `3.2e-2` | Domain height in `z`. |
| `g_accel` | `9.81` | Gravity in `m/s^2`; use `1.62` for lunar gravity tests. |
| `dlo`, `dhi` | `8.0e-5`, `1.2e-4` | Inserted grain diameter range. |
| `rho` | `3300` | Mobile grain density in `kg/m^3`. |
| `npart` | `500000` | Feed budget for `fix pour`, not target in-box count. |
| `tube_diameter` | `1.4e-4` | Diameter assigned to fixed tube beads. |
| `tube_density` | `7800` | Density assigned to fixed tube beads. |
| `bank_zlo`, `bank_zhi` | `5.2e-3`, `2.72e-2` | Tube-bank analysis region. |
| `outlet_zhi` | `1.2e-3` | Upper boundary of outlet deletion region. |
| `outlet_every` | `1000` | Outlet deletion interval in timesteps. |
| `outlet_max` | `100000` | Maximum grains deleted per outlet event. |
| `outlet_seed` | `928371` | Random seed for `fix evaporate`. |
| `kn` | `1.0e3` | Normal spring stiffness for the soft DEM contact model. |
| `gamma` | `2.0` | Normal damping argument for the granular contact model. |
| `mu` | `0.5` | Friction coefficient. |
| `dt` | `2.0e-7 s` | DEM timestep. |
| `nsteps` | `600000` | Total timesteps. |
| `thermo_every` | `20000` | Log output interval. |
| `global_every` | `2000` | Scalar diagnostic interval. |
| `sync_every` | `20000` | Rich mobile-particle dump interval. |
| `contact_every` | `20000` | Local contact dump interval. |
| `traj_every` | `100000` | Visualization dump interval. |
| `dump_start`, `dump_stop` | `200000`, `600000` | Detailed dump window. |

The feed region is near the top:

| Parameter | Value |
| --- | ---: |
| `feed_xlo`, `feed_xhi` | `5.0e-4`, `1.15e-2` |
| `feed_ylo`, `feed_yhi` | `2.0e-4`, `2.2e-3` |
| `feed_zlo`, `feed_zhi` | `2.85e-2`, `3.1e-2` |

The current insertion velocity is set in `in.heat_tube_bank_body`:

```lammps
vel -0.03 0.03 -0.01 0.01 -0.35
```

The first four numbers are lower and upper bounds for `vx` and `vy`; the last
number sets the downward `vz` component. This is a modeling choice and should
be tested.

## LAMMPS Design

The central LAMMPS commands are:

```lammps
units           si
atom_style      sphere
pair_style      gran/hooke/history ${kn} NULL ${gamma} NULL ${mu} 0
fix             int mobile nve/sphere
fix             grav mobile gravity ${g_accel} vector 0 0 -1
fix             feed mobile pour ...
fix             outlet mobile evaporate ...
```

The model uses a deliberately soft DEM stiffness, `kn = 1.0e3 N/m`, for
computational speed.  LAMMPS contact forces should therefore be treated as
DEM-scale collision diagnostics, not calibrated mineral elastic forces.

Only the `mobile` group is integrated by `fix nve/sphere`; the tube beads are
fixed.  Tube-tube neighbor interactions are excluded:

```lammps
neigh_modify delay 0 every 1 check yes exclude group tubes tubes
```

The remaining shared-body settings are summarized below.  These are part of the
baseline numerical model and should be changed only after Matt understands what
each setting controls.

| Command or setting | Purpose and notes |
| --- | --- |
| `newton off` | Conservative setting for granular contacts with tangential history.  It avoids relying on Newton pair-force communication for history-dependent contacts, at the cost of some extra pair work. |
| `comm_modify mode single vel yes` | Communicates velocities for ghost particles, which granular contact models need for damping and relative-velocity calculations. |
| `processors * 1 *` | Keeps only one MPI processor division in the thin `y` direction.  This avoids inefficient decomposition across a very narrow channel. |
| `create_box 2 box` | Creates two atom types: type 1 mobile grains and type 2 fixed tube beads. |
| `include tube_particles.in` | Reads the generated fixed tube-surface beads.  Regenerate this file after changing the tube geometry. |
| `set type 2 diameter ... density ...` | Assigns diameter and density to fixed tube beads.  Their mass is mostly for contact bookkeeping because they are not integrated. |
| `group tubes type 2`, `group mobile type 1` | Separates fixed tube beads from mobile grains.  Only `mobile` is integrated. |
| `pair_coeff * *` | Required coefficient line for the selected granular pair style. |
| `neighbor 3.0e-5 bin` | Neighbor-list skin.  Increase it if LAMMPS reports dangerous builds; do not reduce it without checking stability. |
| `neigh_modify delay 0 every 1 check yes ...` | Updates neighbor lists frequently, checks whether rebuilding is needed, and excludes tube-tube contacts. |
| `group mobile_bank dynamic ...` | Updates the group of grains inside the tube-bank analysis region at the scalar-output interval. |
| `group mobile_outlet dynamic ...` | Updates the group of grains currently inside the outlet region. |
| `compute contact_count mobile contact/atom` | Counts current contacts for each mobile grain. |
| `compute ke_atom`, `compute erot_atom` | Computes translational and rotational kinetic energies per mobile grain. |
| `compute ... reduce ...` | Converts per-grain quantities into scalar sums or maxima for monitoring. |
| `variable pre_steps`, `production_steps`, `post_steps` | Splits the run into pre-dump, detailed-dump, and post-dump stages. |
| `fix global_out all ave/time ...` | Writes scalar diagnostics to `global_metrics.dat` every `global_every` steps. |
| `thermo_style custom ...` | Chooses the quantities printed to `log.lammps`. |
| `thermo_modify lost warn` | Prints warnings instead of stopping if atoms are lost.  The outlet deletion should remove processed grains first, but this keeps the open-flow calculation from failing immediately if an atom leaves the box. |
| `compute property/local ... cutoff radius` | Outputs contact partner ids and types for finite-size overlap contacts. |
| `compute pair/local ... cutoff radius` | Outputs local contact distances, vectors, forces, and pair-style-specific granular quantities. |
| `dump_modify ... sort id` | Sorts dump entries by atom id so trajectories and postprocessing are more reproducible. |
| `undump` and `uncompute` after the detailed window | Stops expensive output and local-contact computes after `dump_stop`. |

## Contact Parameters

The current contact parameters are baseline DEM values, not calibrated lunar
regolith or ilmenite material constants.  They are suitable for testing whether
the moving-bed workflow runs and whether the needed contact statistics can be
extracted, but they must be checked before quantitative interpretation.

The relevant command is:

```lammps
pair_style gran/hooke/history ${kn} NULL ${gamma} NULL ${mu} 0
```

For this Hookean history model, LAMMPS uses a normal spring, a normal damping
term, a tangential history spring, and a Coulomb friction limit.  Because the
second and fourth arguments are `NULL`, LAMMPS uses $K_t = 2K_n/7$ and
$\gamma_t = \gamma_n/2$.  The final `0` means tangential damping is not
included.  The same parameter meanings are used by the side-wall
`fix wall/gran` commands.

`kn = 1.0e3 N/m` is deliberately soft.  It reduces the collision frequency and
allows a practical timestep, but it should be treated as a numerical stiffness,
not as the true elastic stiffness of mineral grains.  If `kn` is too small,
overlaps become unrealistically large and the contact-force and `d_actual`
proxies will change.  If `kn` is increased, the timestep usually needs to be
reduced.

`gamma = 2.0` is a normal damping parameter in the Hookean contact model.  This
value should be treated as a placeholder.  A more physical choice can be
estimated from a target coefficient of restitution, $e$, using the linear
spring-dashpot relation:

$$
\zeta =
\frac{-\ln e}{\sqrt{\pi^2 + \left(\ln e\right)^2}},
\qquad
\gamma_n = 2\zeta\sqrt{\frac{k_n}{m_\mathrm{eff}}}.
$$

where $m_\mathrm{eff}$ is the effective contact mass.  For equal mobile grains,
$m_\mathrm{eff} = m/2$; for a grain-wall collision, $m_\mathrm{eff} = m$.
This estimate is only a calibration aid, so Matt should still verify the
resulting dynamics in LAMMPS.

`mu = 0.5` is a moderate friction baseline.  The correct value depends on grain
composition, angularity, surface roughness, dust coatings, and whether the
contact is grain-grain or grain-wall.  Since this input uses one `mu` for all
contacts, it should be interpreted as an effective friction coefficient.

The contact parameter tests should compare not just the visual flow, but also
the quantities that enter the helium-release model: residence-time
distributions, contact rates, grain-tube versus grain-grain contact fractions,
normal relative-velocity distributions, force distributions, `d_actual`
percentiles, and agitation-energy proxies per processed mass.

## Output Files

### `global_metrics.dat`

This is the cheapest and most important live-monitoring file.  It is written
throughout the simulation.  The columns are:

| Column | Meaning |
| --- | --- |
| `step` | Timestep. |
| `nmobile` | Number of mobile grains in the simulation box. |
| `nbank` | Number of mobile grains in the tube-bank analysis region. |
| `noutlet` | Number of mobile grains currently in the outlet region. |
| `outlet_deleted` | Cumulative number of mobile grains deleted by `fix evaporate`. |
| `sum_contact` | Sum of per-grain contact counts over all mobile grains. |
| `max_contact` | Maximum contact count on any mobile grain. |
| `sum_ke_J` | Total mobile translational kinetic energy. |
| `sum_erot_J` | Total mobile rotational kinetic energy. |
| `bank_contact` | Contact-count sum restricted to grains in the tube bank. |
| `bank_ke_J` | Translational kinetic energy restricted to the tube bank. |
| `bank_erot_J` | Rotational kinetic energy restricted to the tube bank. |

A useful run should show a startup period followed by a statistical plateau in
`nmobile`, `nbank`, contact counts, kinetic energies, and the slope of
`outlet_deleted`.

### `particles_sync.dump`

This rich mobile-particle dump is synchronized with `contact_events.local`.  It
contains positions, velocities, angular velocities, forces, torques, per-grain
contact counts, translational kinetic energy, and rotational kinetic energy.

### `contact_events.local`

This local dump contains grain-grain and grain-tube contact information from
`compute property/local` and `compute pair/local` with `cutoff radius`.  The
`cutoff radius` option restricts records to finite-size particles whose centers
are within the sum of their radii.

### `trajectory_viz.dump`

This lower-frequency dump contains both mobile grains and fixed tube beads.  It
is intended for visualization.

## Postprocessing

Install Python dependencies:

```bash
python3 -m pip install -r requirements.txt
```

Run:

```bash
python3 analyze_dem.py /path/to/run_directory
```

The script writes:

- `summary.csv`
- `bed_time_series.csv`
- `contact_time_series.csv`
- `grain_summary.csv`
- `global_metrics_clean.csv`
- `figures/flow_metrics.pdf` and `.png`
- `figures/contact_metrics.pdf` and `.png`
- `figures/grain_distributions.pdf` and `.png`
- `figures/snapshot_dactual.pdf` and `.png`

The postprocessor estimates Hertz-style contact radii and mechanically affected
depth proxies.  It uses preliminary elastic constants:

$$
E_\mathrm{grain} = 120~\mathrm{GPa}, \qquad
\nu_\mathrm{grain} = 0.25,
$$

$$
E_\mathrm{tube} = 200~\mathrm{GPa}, \qquad
\nu_\mathrm{tube} = 0.29.
$$

For a contact between particles 1 and 2,

$$
R^* =
\left(\frac{1}{R_1} + \frac{1}{R_2}\right)^{-1},
\qquad
\frac{1}{E^*} =
\frac{1-\nu_1^2}{E_1} + \frac{1-\nu_2^2}{E_2}.
$$

The force-based contact radius is

$$
a_F =
\left(\frac{3F_{n}R^{*}}{4E^{*}}\right)^{1/3}.
$$

The velocity-based contact radius is estimated from

$$
\delta_\mathrm{max} =
\left(
\frac{15m^{*}v_n^{2}}{16E^{*}\sqrt{R^{*}}}
\right)^{2/5},
\qquad
a_{v} = \sqrt{R^{*}\delta_{\mathrm{max}}}.
$$

The affected-depth proxy is

$$
d_\mathrm{actual} = \frac{\chi a}{2},
$$

where $\chi$ is an uncertain calibration factor.  The script reports
$\chi = 0.05$, $0.20$, and $1.00$; summary columns use $\chi = 0.20$.

## Running on Adroit

Please modify the path to LAMMPS binary file and the required modules for LAMMPS before submitting the SLURM script. The
script requests Intel Ice Lake nodes with:

```text
#SBATCH --constraint=ice
```
This line can be removed if the number of available Ice Lake cores are limited (check using `shownodes` command).
The number of requested CPU cores (`--ntasks`) can also be adjusted based on the available resources.

Example run:

```bash
cd heat_tube_bank_full_dem
python3 make_tube_particles.py
sbatch run_adroit.sbatch"
```

Monitor a job:

```bash
ssh adroit "squeue -u ${NETID} --format='%.18i %.9P %.12j %.2t %.10M %.10l %.6D %Z'"
```

During a run, inspect `global_metrics.dat` first.  It is smaller and more
reliable for progress monitoring than waiting for `log.lammps` to flush.

## What Matt Should Check First

1. Run the baseline input and confirm that LAMMPS writes `global_metrics.dat`,
   `particles_sync.dump`, `contact_events.local`, `trajectory_viz.dump`, and
   `log.lammps`.

2. Confirm that outlet deletion is active.  In `global_metrics.dat`,
   `outlet_deleted` should eventually increase.  If it does not, the grains
   have not yet reached the outlet or the run is too short.

3. Identify whether the simulation reaches a statistical moving-bed regime.
   Look for a period where `nmobile`, `nbank`, contact counts, kinetic energies, and the
   slope of `outlet_deleted` fluctuate around a stable time average. If not, the total simulation
   time should be increased and the dump/statistics window should be delayed.

4. Run `analyze_dem.py` on the run directory and inspect the generated CSV files
   and figures.

5. Design convergence tests before choosing production settings.

## Convergence Tests Before Production

Matt can decide production settings after checking whether the main outputs are converged with respect to
numerical and finite-size choices. Important tests include:

- Larger/smaller domains or tube banks with different geometry.
- Total simulation time.
- Number of feed particles `npart`.
- Contact stiffness tests, for example `kn = 1e3`, `3e3`, and `1e4`.
  For each `kn`, reduce `dt` if needed and check that overlaps are small,
  neighbor-list warnings are absent, and the main outputs do not change
  strongly with further stiffening.
- Damping tests based on coefficient of restitution. Choose several target
  values such as $e = 0.5$, $0.8$, and $0.95$, convert them to $\gamma_n$ using
  the linear spring-dashpot estimate in the Contact Parameters section, and
  compare residence times, collision velocities, contact forces, and processed
  mass rate.
- Friction tests, for example `mu = 0.3`, `0.5`, and `0.8`.  This is important
  because friction controls tangential work, rotational motion, clogging, and
  the grain-tube contact history.
- Optional tangential-damping model test.  The current final argument in
  `gran/hooke/history` is `0`; changing it to `1` includes tangential damping
  and should be treated as a separate model choice.
- Grain-size ranges `dlo/dhi`.
- Feed velocity.
- `contact_every`, `sync_every`, and `traj_every`.
- `dump_start` and `dump_stop` values after checking when the
  moving-bed statistics become stable.
- Different random seeds for `fix pour` and `fix evaporate`.
- Earth gravity, `g_accel = 9.81 m/s^2`, versus lunar gravity,
  `g_accel = 1.62 m/s^2`.

The quantities to compare are:

- tube-bank residence-time distributions,
- processed mass rate inferred from `outlet_deleted`,
- contact counts and contact-type ratios,
- force and relative-velocity distributions,
- surface coverage,
- `d_actual` percentiles,
- agitation-energy proxies per processed mass.

## Practical Notes

- Keep `contact_every` and `sync_every` equal when velocity-based contact
  postprocessing is needed for every sampled contact.
- Increase `sync_every` and `contact_every` if dump files become too large.
- If LAMMPS reports many dangerous neighbor builds, reduce `dt` or increase the
  neighbor skin.
- If `fix pour` reports many failed insertions, enlarge the feed region, reduce
  insertion density, reduce feed rate, or revisit the grain-size distribution.
- Do not use periodic `x` or `y` boundaries unless the side-wall model is also
  redesigned.  The current `wallx` and `wally` fixes require non-periodic
  directions.
- The model is mechanical only.  It does not solve heat transfer, helium
  diffusion, or helium release kinetics.

## References and Documentation

- Olson, A. D. S. "Lunar Helium-3: Mining Concepts, Extraction Research, and
  Potential ISRU Synergies." AIAA ASCEND, 2021. NASA NTRS:
  <https://ntrs.nasa.gov/api/citations/20210022801/downloads/AIAA%20ASCEND%202021%20Paper_211018.pdf>
- LAMMPS `fix pour` documentation:
  <https://docs.lammps.org/fix_pour.html>
- LAMMPS `fix evaporate` documentation:
  <https://docs.lammps.org/fix_evaporate.html>
- LAMMPS `fix wall/gran` documentation:
  <https://docs.lammps.org/fix_wall_gran.html>
- LAMMPS granular pair-style documentation:
  <https://docs.lammps.org/pair_gran.html>
- LAMMPS local contact output documentation:
  <https://docs.lammps.org/compute_pair_local.html>
