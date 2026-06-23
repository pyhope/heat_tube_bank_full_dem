# Analysis Outputs

This note explains what `analyze_dem.py` reads, what it writes, and how the
main plotted quantities are calculated.  The current analysis is meant for
preliminary diagnostics of the tube-bank DEM model, not as a final calibrated
helium-release model.

## Input Files

`analyze_dem.py` expects a LAMMPS run directory containing:

- `particles_sync.dump`: synchronized mobile-particle frames with positions,
  velocities, angular velocities, radii, masses, contact counts, and kinetic
  energies.
- `contact_events.local`: local pair-contact data from `compute pair/local`.
- `global_metrics.dat`: inexpensive scalar diagnostics from `fix ave/time`.
- `trajectory_viz.dump`: optional visualization trajectory for the final
  spatial snapshot.
- `in.heat_tube_bank` and `in.heat_tube_bank_body`: used to read model
  parameters such as `dt`, `bank_zlo`, `bank_zhi`, `lx`, and `lz`.

The raw dump files may be too large for GitHub.  The small reference files in
`results/` are included for inspection, but rerunning the full analysis from
scratch requires the raw LAMMPS dumps.

## CSV Outputs

The script writes:

- `bed_time_series.csv`: time series from `particles_sync.dump`, including
  particle counts, bank-region counts, total mobile mass, kinetic energies, and
  vertical-position percentiles.
- `contact_time_series.csv`: per-contact-frame counts, contact type counts,
  normal-force percentiles, and force-based `d_actual` percentiles.
- `grain_summary.csv`: per-grain residence time, contact counts, provisional
  sampled area-coverage proxies, maximum force-based and velocity-based
  `d_actual` values, and collision-energy proxies.
- `summary.csv`: one-row aggregate summary.
- `global_metrics_clean.csv`: cleaned version of `global_metrics.dat`.

## Common Definitions

The tube-bank region is defined by

$$
z_\mathrm{bank,lo} \le z_i \le z_\mathrm{bank,hi}.
$$

The script uses preliminary elastic constants

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

The force-based Hertz contact radius is

$$
a_F =
\left(\frac{3F_{n}R^{\ast}}{4E^{\ast}}\right)^{1/3}.
$$

The mechanically affected depth proxy is

$$
d_\mathrm{actual} = \frac{\chi a}{2},
$$

with the default reporting value `chi = 0.20`.  This is a mechanical proxy,
not a directly validated helium-release depth.

The velocity-based contact-radius estimate uses the normal relative velocity at
the contact point:

$$
\delta_\mathrm{max} =
\left(
\frac{15m^{\ast}v_n^{2}}{16E^{\ast}\sqrt{R^{\ast}}}
\right)^{2/5},
\qquad
a_{v} = \sqrt{R^{\ast}\delta_{\mathrm{max}}}.
$$

Here `m*` is the effective contact mass and `v_n` is the normal component of
the contact-point relative velocity.

## `flow_metrics`

This figure summarizes moving-bed behavior over time.  The shared x-axis is
simulation time in milliseconds:

$$
t = \mathrm{step}\times dt.
$$

Panel 1, `Particles`, reports:

- `mobile`: number of mobile grains in the simulation box.
- `in bank`: number of mobile grains whose positions satisfy the tube-bank
  z-range criterion.

Panel 2, `Flow (mg/s)`, now uses the outlet-deletion counter from
`global_metrics.dat`.  If `N_deleted(t)` is the cumulative number of grains
removed by the outlet deletion fix and `mbar` is the mean mobile-grain mass,
the processed mass is

$$
M_\mathrm{processed}(t) = N_\mathrm{deleted}(t)\bar{m}.
$$

The plotted processed-flow rate is

$$
\dot{M}_\mathrm{processed}(t)
=
\frac{dM_\mathrm{processed}}{dt}.
$$

This is more direct than estimating flow from the changing mass still present
inside the box.

Panel 3, `Energy`, reports the total mobile translational and rotational
kinetic energies from `ke/atom` and `erotate/sphere/atom`, summed over mobile
grains and plotted in microjoules.

Panel 4, `-v_z`, reports the negative of the mean z velocity of grains in the
tube-bank region:

$$
-\langle v_z\rangle_\mathrm{bank}.
$$

Since grains flow downward, `v_z` is typically negative; plotting `-v_z` makes
downward speed positive.

## `contact_metrics`

This figure reports contact statistics as a function of time.  The x-axis is
the contact-frame time in milliseconds.

Panel 1, `Contacts`, reports:

- `grain-grain`: contacts between two mobile grains.
- `grain-tube`: contacts between a mobile grain and a fixed tube bead.

Tube-tube contacts are ignored.

Panel 2, `Force (mN)`, reports the 50th and 90th percentiles of the LAMMPS
normal pair force `force` from `compute pair/local`.  The script no longer
replaces non-positive force values with a Hookean overlap estimate.  If the
LAMMPS force value is non-positive, it is set to zero for this diagnostic.

Panel 3, `d_actual`, reports the 50th and 90th percentiles of the force-based
affected-depth proxy, calculated from the LAMMPS normal force using the
force-based Hertz radius `a_F`.

Percentiles are calculated across all contacts present in a given sampled
contact frame.  They are not error bars and they are not statistics across
independent simulations.

## `grain_distributions`

This figure now contains two panels.

The left panel shows the distribution of `d_actual` values across sampled
contacts.  The force-based curve uses `a_F` from the LAMMPS normal force.  When
velocity-resolved contact data are available, the velocity-based curve uses
`a_v` from the normal relative velocity estimate.  The x-axis is logarithmic
because the affected-depth proxy spans a narrow but important small-length
range.

The right panel shows the distribution of tube-bank residence times.  Each time
a grain appears in a synchronized particle frame inside the tube-bank z-range,
the script adds

$$
\Delta t = dt \times \mathrm{sync\_every}
$$

to that grain's residence time.  The plotted value is therefore a
sampled-frame residence time, not a continuously tracked residence time.

## `snapshot_dactual`

If `trajectory_viz.dump` is available, the script plots the final visualization
frame in the x-z plane.  The x-axis and z-axis are positions in millimeters.
Fixed tube beads are plotted in dark gray.  Mobile grains are colored by the
largest force-based `d_actual` value that the grain experienced during the
sampled analysis window:

$$
d_{\mathrm{actual},i}^{\max}
=
\max_j d_{\mathrm{actual},ij}.
$$

This color is a history value assigned to the grain and then shown at its final
snapshot position.  It is not the instantaneous contact depth at the final
frame.

## Current Status of Area Coverage

`grain_summary.csv` still contains provisional columns named
`surface_coverage_force` and `surface_coverage_velocity`.  These are computed
by summing sampled contact patch areas and dividing by the grain surface area:

$$
C_i^\mathrm{sampled}=\frac{\sum_j \pi a_{ij}^2}{4\pi R_i^2}.
$$

This quantity is not currently plotted because it is not yet a reliable
physical surface-coverage result.

There are two main problems.

First, the result is strongly dependent on the contact-output interval.  The
current contact dump can miss short-lived contacts between sampled frames.  If
contacts were written every timestep and the same formula were used, the same
physical contact could be counted many times while it persists across multiple
timesteps.  Therefore the sampled sum can either undercount missed contacts or
overcount persistent contacts, depending on the output strategy.

Second, the current script does not know where on a rotating grain surface each
contact patch occurs.  It simply adds patch areas.  Real unique coverage
depends on whether new contacts hit new surface regions or repeatedly hit
previously contacted regions.

The right way to report area coverage will require additional simulations and
postprocessing:

1. Run a short high-frequency contact-output window, with `contact_every`
   small enough to resolve contact lifetimes.
2. Group repeated appearances of the same grain pair into contact episodes,
   from first contact to separation.
3. For each episode, record quantities such as maximum normal force, maximum
   contact radius, duration, normal impulse, and tangential work.
4. Count each episode once when estimating accumulated contact area, rather
   than counting every sampled frame as a new patch.
5. Develop a surface-mixing model.  A simple upper-bound/random-coverage model
   would use

$$
F_\mathrm{covered}=1-\exp\left(-\frac{\sum_j A_j}{4\pi R^2}\right),
$$

where `A_j` is the effective patch area of episode `j`.  This still assumes
that patches are randomly distributed over the grain surface; repeated contacts
on the same side of a grain would produce lower unique coverage.

Until this workflow is implemented, the area-coverage columns should be treated
only as diagnostic proxies and should not be used as final helium-release
metrics.
