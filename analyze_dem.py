#!/usr/bin/env python3
from pathlib import Path
import re
import sys

import numpy as np
from matplotlib import pyplot as plt

root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
figdir = root / "figures"
figdir.mkdir(exist_ok=True)

particle_dump = root / "particles_sync.dump"
contact_dump = root / "contact_events.local"
viz_dump = root / "trajectory_viz.dump"
global_metrics = root / "global_metrics.dat"

top_input = root / "in.heat_tube_bank"
body_input = root / "in.heat_tube_bank_body"
input_text = ""
if top_input.exists():
    input_text += top_input.read_text() + "\n"
if body_input.exists():
    input_text += body_input.read_text() + "\n"

def lmp_var(name, default):
    m = re.search(rf"^variable\s+{name}\s+equal\s+([0-9.eE+-]+)", input_text, re.M)
    if m:
        return float(m.group(1))
    return default

dt = lmp_var("dt", 2.0e-7)
sync_every = int(lmp_var("sync_every", 500))
contact_every = int(lmp_var("contact_every", 500))
npart = int(lmp_var("npart", 0))
rho = lmp_var("rho", 3300.0)
kn = lmp_var("kn", 1.0e3)
bank_zlo = lmp_var("bank_zlo", 5.2e-3)
bank_zhi = lmp_var("bank_zhi", 1.48e-2)
lx = lmp_var("lx", 6.0e-3)
lz = lmp_var("lz", 2.0e-2)
sample_dt = dt * contact_every

E_grain = 120.0e9
nu_grain = 0.25
E_tube = 200.0e9
nu_tube = 0.29
chi_values = np.array([0.05, 0.20, 1.00])
chi_ref = 0.20

def material(atom_type):
    if int(atom_type) == 2:
        return E_tube, nu_tube
    return E_grain, nu_grain

def effective_modulus(t1, t2):
    E1, nu1 = material(t1)
    E2, nu2 = material(t2)
    return 1.0 / ((1.0 - nu1 * nu1) / E1 + (1.0 - nu2 * nu2) / E2)

def read_dump(path):
    with path.open() as fh:
        while True:
            line = fh.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            step = int(fh.readline())
            fh.readline()
            n = int(fh.readline())
            fh.readline()
            for _ in range(3):
                fh.readline()
            cols = fh.readline().split()[2:]
            data = np.empty((n, len(cols)))
            for i in range(n):
                data[i] = [float(x) for x in fh.readline().split()]
            yield step, cols, data

def read_global_metrics(path):
    if not path.exists():
        return None, None
    rows = []
    cols = None
    name_map = {
        "TimeStep": "step",
        "v_nmobile": "nmobile",
        "v_nbank": "nbank",
        "v_noutlet": "noutlet",
        "f_outlet": "outlet_deleted",
        "c_sum_contact": "sum_contact",
        "c_max_contact": "max_contact",
        "c_sum_ke": "sum_ke_J",
        "c_sum_erot": "sum_erot_J",
        "c_bank_contact": "bank_contact",
        "c_bank_ke": "bank_ke_J",
        "c_bank_erot": "bank_erot_J",
    }
    with path.open() as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                tokens = stripped.lstrip("#").strip().split()
                if tokens and tokens[0] == "TimeStep":
                    cols = [name_map.get(tok, tok) for tok in tokens]
                continue
            rows.append([float(x) for x in stripped.split()])
    if not rows:
        return None, None
    data = np.asarray(rows)
    if cols is None or len(cols) != data.shape[1]:
        base_cols = [
            "step", "nmobile", "nbank", "sum_contact", "max_contact",
            "sum_ke_J", "sum_erot_J", "bank_contact", "bank_ke_J", "bank_erot_J",
        ]
        if data.shape[1] <= len(base_cols):
            cols = base_cols[:data.shape[1]]
        else:
            cols = base_cols + [f"extra_{i}" for i in range(data.shape[1] - len(base_cols))]
    return cols, data

def percentile(values, qs):
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return np.full(len(qs), np.nan)
    return np.percentile(arr, qs)

def write_csv(path, header, rows):
    arr = np.asarray(rows, dtype=float)
    np.savetxt(path, arr, delimiter=",", header=",".join(header), comments="")

radius = {}
mass = {}
atom_type = {}
bank_time = {}
bank_entries = {}
prev_bank_ids = set()
bed_rows = []
last_mobile = None

if not particle_dump.exists():
    raise SystemExit(f"Missing {particle_dump}")

for step, cols, data in read_dump(particle_dump):
    idx = {name: i for i, name in enumerate(cols)}
    ids = data[:, idx["id"]].astype(np.int64)
    types = data[:, idx["type"]].astype(np.int64)
    z = data[:, idx["z"]]
    vx = data[:, idx["vx"]]
    vy = data[:, idx["vy"]]
    vz = data[:, idx["vz"]]
    r = data[:, idx["radius"]]
    m = data[:, idx["mass"]]
    cnum = data[:, idx.get("c_contact_count", idx.get("c_contact_count[1]", 0))]
    ke = data[:, idx.get("c_ke_atom", idx.get("c_ke_atom[1]", 0))]
    erot = data[:, idx.get("c_erot_atom", idx.get("c_erot_atom[1]", 0))]

    for atom_id, typ, ri, mi in zip(ids, types, r, m):
        radius[int(atom_id)] = float(ri)
        mass[int(atom_id)] = float(mi)
        atom_type[int(atom_id)] = int(typ)

    bank = (z >= bank_zlo) & (z <= bank_zhi)
    bank_ids = set(int(v) for v in ids[bank])
    for atom_id in bank_ids:
        bank_time[atom_id] = bank_time.get(atom_id, 0.0) + dt * sync_every
        if atom_id not in prev_bank_ids:
            bank_entries[atom_id] = bank_entries.get(atom_id, 0) + 1
    prev_bank_ids = bank_ids

    speed = np.sqrt(vx * vx + vy * vy + vz * vz)
    in_bank = int(np.sum(bank))
    total_mass = float(np.sum(m))
    bank_mass = float(np.sum(m[bank])) if in_bank else 0.0
    mean_vz_bank = float(np.mean(vz[bank])) if in_bank else np.nan
    mean_speed_bank = float(np.mean(speed[bank])) if in_bank else np.nan
    z10, z50, z90 = np.percentile(z, [10, 50, 90]) if len(z) else (np.nan, np.nan, np.nan)
    bed_rows.append([
        step, step * dt, len(ids), in_bank, total_mass, bank_mass, mean_vz_bank,
        mean_speed_bank, float(np.sum(ke)), float(np.sum(erot)), float(np.sum(cnum)),
        float(np.max(cnum)) if len(cnum) else 0.0, z10, z50, z90,
    ])
    last_mobile = data.copy()
    last_cols = cols

bed_header = [
    "step", "time_s", "n_mobile", "n_bank", "mobile_mass_kg", "bank_mass_kg",
    "mean_vz_bank_m_s", "mean_speed_bank_m_s", "ke_trans_J", "ke_rot_J",
    "sum_contact_count", "max_contact_count", "z10_m", "z50_m", "z90_m",
]
write_csv(root / "bed_time_series.csv", bed_header, bed_rows)

bed_arr = np.asarray(bed_rows)

contact_cols = None
contacts_by_step = {}
force_samples = []
ft_samples = []
d_force_samples = []
contact_time_rows = []
contact_count = {}
grain_contact_count = {}
tube_contact_count = {}
coverage_force = {}
coverage_velocity = {}
max_d_force = {chi: {} for chi in chi_values}
max_d_velocity = {chi: {} for chi in chi_values}
contact_work = {}
normal_energy = {}
total_contacts = 0
total_grain_contacts = 0
total_tube_contacts = 0
missing_radius = 0

default_radius = {1: 0.5 * (lmp_var("dlo", 8.0e-5) + lmp_var("dhi", 1.2e-4)), 2: 0.5 * lmp_var("tube_diameter", 1.4e-4)}
default_mass = {1: 4.0 / 3.0 * np.pi * default_radius[1] ** 3 * rho, 2: 1.0e30}

if contact_dump.exists():
    for step, cols, data in read_dump(contact_dump):
        contact_cols = cols
        idx = {name: i for i, name in enumerate(cols)}
        contacts_by_step[step] = (cols, data)
        frame_forces = []
        frame_d = []
        frame_tube = 0
        frame_grain = 0

        for row in data:
            id1 = int(row[idx["c_cids[1]"]])
            id2 = int(row[idx["c_cids[2]"]])
            t1 = int(row[idx["c_cids[3]"]])
            t2 = int(row[idx["c_cids[4]"]])
            if t1 == 2 and t2 == 2:
                continue
            r1 = radius.get(id1, default_radius.get(t1, default_radius[1]))
            r2 = radius.get(id2, default_radius.get(t2, default_radius[1]))
            if id1 not in radius or id2 not in radius:
                missing_radius += 1
            dist = float(row[idx["c_cpair[1]"]])
            overlap = max(r1 + r2 - dist, 0.0)
            if overlap <= 0.0:
                continue
            reff = 1.0 / (1.0 / r1 + 1.0 / r2)
            estar = effective_modulus(t1, t2)
            fn = float(row[idx["c_cpair[5]"]])
            if fn <= 0.0:
                fn = 0.0
            ft = 0.0
            if "c_cpair[12]" in idx:
                ft = abs(float(row[idx["c_cpair[12]"]]))
            elif all(f"c_cpair[{i}]" in idx for i in (9, 10, 11)):
                ft = float(np.sqrt(sum(row[idx[f"c_cpair[{i}]"]] ** 2 for i in (9, 10, 11))))
            a_force = (3.0 * fn * reff / (4.0 * estar)) ** (1.0 / 3.0)
            d_force = chi_ref * a_force * 0.5
            force_samples.append(fn)
            ft_samples.append(ft)
            d_force_samples.append(d_force * 1.0e6)
            frame_forces.append(fn)
            frame_d.append(d_force * 1.0e6)
            total_contacts += 1
            if t1 == 1 and t2 == 1:
                total_grain_contacts += 1
                frame_grain += 1
            else:
                total_tube_contacts += 1
                frame_tube += 1

            for atom_id, typ, other_typ, ri in ((id1, t1, t2, r1), (id2, t2, t1, r2)):
                if typ != 1:
                    continue
                contact_count[atom_id] = contact_count.get(atom_id, 0) + 1
                if other_typ == 2:
                    tube_contact_count[atom_id] = tube_contact_count.get(atom_id, 0) + 1
                else:
                    grain_contact_count[atom_id] = grain_contact_count.get(atom_id, 0) + 1
                coverage_force[atom_id] = coverage_force.get(atom_id, 0.0) + np.pi * a_force * a_force
                for chi in chi_values:
                    val = chi * a_force * 0.5
                    max_d_force[chi][atom_id] = max(max_d_force[chi].get(atom_id, 0.0), val)

        fq = percentile(frame_forces, [50, 90, 99])
        dq = percentile(frame_d, [50, 90, 99])
        contact_time_rows.append([
            step, step * dt, len(frame_forces), frame_grain, frame_tube,
            fq[0], fq[1], fq[2], dq[0], dq[1], dq[2],
        ])

contact_header = [
    "step", "time_s", "n_contacts", "n_grain_contacts", "n_tube_contacts",
    "force_p50_N", "force_p90_N", "force_p99_N",
    "d_force_p50_um_chi0p2", "d_force_p90_um_chi0p2", "d_force_p99_um_chi0p2",
]
write_csv(root / "contact_time_series.csv", contact_header, contact_time_rows)

velocity_samples = []
vt_samples = []
d_vel_samples = []
work_samples = []
velocity_contacts = 0

for step, cols, data in read_dump(particle_dump):
    if step not in contacts_by_step:
        continue
    ccols, contacts = contacts_by_step[step]
    if len(contacts) == 0:
        continue
    pidx = {name: i for i, name in enumerate(cols)}
    cidx = {name: i for i, name in enumerate(ccols)}
    ids = data[:, pidx["id"]].astype(np.int64)
    id_to_i = {int(v): i for i, v in enumerate(ids)}
    pos = data[:, [pidx["x"], pidx["y"], pidx["z"]]]
    vel = data[:, [pidx["vx"], pidx["vy"], pidx["vz"]]]
    omega = data[:, [pidx["omegax"], pidx["omegay"], pidx["omegaz"]]]
    radii = data[:, pidx["radius"]]
    masses = data[:, pidx["mass"]]

    for row in contacts:
        id1 = int(row[cidx["c_cids[1]"]])
        id2 = int(row[cidx["c_cids[2]"]])
        t1 = int(row[cidx["c_cids[3]"]])
        t2 = int(row[cidx["c_cids[4]"]])
        if id1 not in id_to_i or id2 not in id_to_i or (t1 == 2 and t2 == 2):
            continue
        i = id_to_i[id1]
        j = id_to_i[id2]
        rij = pos[i] - pos[j]
        dist = np.linalg.norm(rij)
        if dist <= 0.0:
            continue
        n = rij / dist
        r1 = radii[i]
        r2 = radii[j]
        reff = 1.0 / (1.0 / r1 + 1.0 / r2)
        estar = effective_modulus(t1, t2)
        vi_c = vel[i] - r1 * np.cross(omega[i], n)
        vj_c = vel[j] + r2 * np.cross(omega[j], n)
        vrel = vi_c - vj_c
        vn = abs(float(np.dot(vrel, n)))
        vt = float(np.linalg.norm(vrel - np.dot(vrel, n) * n))
        if t1 == 1 and t2 == 1:
            meff = masses[i] * masses[j] / (masses[i] + masses[j])
        elif t1 == 1:
            meff = masses[i]
        else:
            meff = masses[j]
        if vn <= 0.0:
            continue
        delta = ((15.0 * meff * vn * vn) / (16.0 * estar * np.sqrt(reff))) ** 0.4
        a_vel = np.sqrt(reff * delta)
        d_vel = chi_ref * a_vel * 0.5
        velocity_samples.append(vn)
        vt_samples.append(vt)
        d_vel_samples.append(d_vel * 1.0e6)
        velocity_contacts += 1

        ft = abs(float(row[cidx["c_cpair[12]"]])) if "c_cpair[12]" in cidx else 0.0
        work = ft * vt * sample_dt
        e_norm = 0.5 * meff * vn * vn
        work_samples.append(work)
        for atom_id, typ, ri in ((id1, t1, r1), (id2, t2, r2)):
            if typ != 1:
                continue
            coverage_velocity[atom_id] = coverage_velocity.get(atom_id, 0.0) + np.pi * a_vel * a_vel
            contact_work[atom_id] = contact_work.get(atom_id, 0.0) + 0.5 * work
            normal_energy[atom_id] = normal_energy.get(atom_id, 0.0) + 0.5 * e_norm
            for chi in chi_values:
                val = chi * a_vel * 0.5
                max_d_velocity[chi][atom_id] = max(max_d_velocity[chi].get(atom_id, 0.0), val)

mobile_ids = sorted(i for i, t in atom_type.items() if t == 1)
grain_rows = []
for atom_id in mobile_ids:
    ri = radius.get(atom_id, default_radius[1])
    surface = 4.0 * np.pi * ri * ri
    grain_rows.append([
        atom_id, ri, mass.get(atom_id, default_mass[1]), bank_time.get(atom_id, 0.0),
        bank_entries.get(atom_id, 0), contact_count.get(atom_id, 0),
        grain_contact_count.get(atom_id, 0), tube_contact_count.get(atom_id, 0),
        min(1.0, coverage_force.get(atom_id, 0.0) / surface),
        min(1.0, coverage_velocity.get(atom_id, 0.0) / surface),
        max_d_force[0.05].get(atom_id, 0.0) * 1.0e6,
        max_d_force[0.20].get(atom_id, 0.0) * 1.0e6,
        max_d_force[1.00].get(atom_id, 0.0) * 1.0e6,
        max_d_velocity[0.05].get(atom_id, 0.0) * 1.0e6,
        max_d_velocity[0.20].get(atom_id, 0.0) * 1.0e6,
        max_d_velocity[1.00].get(atom_id, 0.0) * 1.0e6,
        contact_work.get(atom_id, 0.0), normal_energy.get(atom_id, 0.0),
    ])

grain_header = [
    "id", "radius_m", "mass_kg", "bank_residence_time_s", "bank_entries",
    "n_contacts", "n_grain_contacts", "n_tube_contacts",
    "surface_coverage_force", "surface_coverage_velocity",
    "max_d_force_um_chi0p05", "max_d_force_um_chi0p2", "max_d_force_um_chi1p0",
    "max_d_velocity_um_chi0p05", "max_d_velocity_um_chi0p2", "max_d_velocity_um_chi1p0",
    "tangential_work_proxy_J", "normal_collision_energy_proxy_J",
]
write_csv(root / "grain_summary.csv", grain_header, grain_rows)

grain_arr = np.asarray(grain_rows) if grain_rows else np.zeros((0, len(grain_header)))
d_force_q = percentile(d_force_samples, [10, 50, 90, 99])
d_vel_q = percentile(d_vel_samples, [10, 50, 90, 99])
force_q = percentile(force_samples, [10, 50, 90, 99])
vn_q = percentile(velocity_samples, [10, 50, 90, 99])
res_q = percentile(grain_arr[:, 3], [10, 50, 90, 99]) if len(grain_arr) else np.full(4, np.nan)
cov_f_q = percentile(grain_arr[:, 8], [10, 50, 90, 99]) if len(grain_arr) else np.full(4, np.nan)
cov_v_q = percentile(grain_arr[:, 9], [10, 50, 90, 99]) if len(grain_arr) else np.full(4, np.nan)
total_mobile_mass = float(np.sum(grain_arr[:, 2])) if len(grain_arr) else np.nan
total_work = float(np.sum(work_samples)) if work_samples else 0.0
total_normal = float(np.sum(list(normal_energy.values()))) if normal_energy else 0.0
if len(grain_arr):
    mean_grain_mass = float(np.nanmean(grain_arr[:, 2]))
else:
    dlo = lmp_var("dlo", 8.0e-5)
    dhi = lmp_var("dhi", 1.2e-4)
    if abs(dhi - dlo) > 0.0:
        mean_d3 = (dhi ** 4 - dlo ** 4) / (4.0 * (dhi - dlo))
    else:
        mean_d3 = dlo ** 3
    mean_grain_mass = rho * np.pi / 6.0 * mean_d3

summary_header = [
    "n_mobile_seen", "npart_requested", "n_contact_records", "n_velocity_contact_records",
    "n_grain_contacts", "n_tube_contacts", "missing_radius_records",
    "d_force_p10_um_chi0p2", "d_force_p50_um_chi0p2", "d_force_p90_um_chi0p2", "d_force_p99_um_chi0p2",
    "d_velocity_p10_um_chi0p2", "d_velocity_p50_um_chi0p2", "d_velocity_p90_um_chi0p2", "d_velocity_p99_um_chi0p2",
    "force_p10_N", "force_p50_N", "force_p90_N", "force_p99_N",
    "vn_p10_m_s", "vn_p50_m_s", "vn_p90_m_s", "vn_p99_m_s",
    "residence_p10_s", "residence_p50_s", "residence_p90_s", "residence_p99_s",
    "coverage_force_p10", "coverage_force_p50", "coverage_force_p90", "coverage_force_p99",
    "coverage_velocity_p10", "coverage_velocity_p50", "coverage_velocity_p90", "coverage_velocity_p99",
    "total_mobile_mass_seen_kg", "tangential_work_proxy_total_J", "normal_collision_energy_proxy_total_J",
    "tangential_work_proxy_J_per_kg", "normal_collision_energy_proxy_J_per_kg",
]
summary = np.array([
    len(mobile_ids), npart, total_contacts, velocity_contacts, total_grain_contacts, total_tube_contacts, missing_radius,
    *d_force_q, *d_vel_q, *force_q, *vn_q, *res_q, *cov_f_q, *cov_v_q,
    total_mobile_mass, total_work, total_normal,
    total_work / total_mobile_mass if total_mobile_mass and np.isfinite(total_mobile_mass) else np.nan,
    total_normal / total_mobile_mass if total_mobile_mass and np.isfinite(total_mobile_mass) else np.nan,
], dtype=float)
write_csv(root / "summary.csv", summary_header, summary.reshape(1, -1))

gcols, gdata = read_global_metrics(global_metrics)
if gdata is not None:
    np.savetxt(root / "global_metrics_clean.csv", gdata, delimiter=",", header=",".join(gcols), comments="")

flow_t_ms = np.array([])
processed_flow_mg_s = np.array([])
if gdata is not None and "step" in gcols and "outlet_deleted" in gcols:
    idx_g = {name: i for i, name in enumerate(gcols)}
    flow_t_s = gdata[:, idx_g["step"]] * dt
    processed_mass = gdata[:, idx_g["outlet_deleted"]] * mean_grain_mass
    if len(flow_t_s) > 1:
        processed_flow = np.gradient(processed_mass, flow_t_s, edge_order=1)
    else:
        processed_flow = np.zeros_like(flow_t_s)
    processed_flow[~np.isfinite(processed_flow)] = 0.0
    processed_flow[processed_flow < 0.0] = 0.0
    if len(bed_arr):
        keep = (flow_t_s >= bed_arr[0, 1]) & (flow_t_s <= bed_arr[-1, 1])
    else:
        keep = np.ones_like(flow_t_s, dtype=bool)
    flow_t_ms = flow_t_s[keep] * 1.0e3
    processed_flow_mg_s = processed_flow[keep] * 1.0e6

def finish_axis(ax):
    ax.minorticks_on()
    ax.tick_params(direction="in", top=True, right=True)

t_ms = bed_arr[:, 1] * 1.0e3
fig, ax = plt.subplots(4, 1, figsize=(5.6, 6.8), sharex=True)
ax[0].plot(t_ms, bed_arr[:, 2], "C0", lw=1.2, label="mobile")
ax[0].plot(t_ms, bed_arr[:, 3], "C1", lw=1.2, label="in bank")
if len(processed_flow_mg_s):
    ax[1].plot(flow_t_ms, processed_flow_mg_s, "C3", lw=1.2)
ax[2].plot(t_ms, bed_arr[:, 8] * 1.0e6, "C4", lw=1.2, label="trans")
ax[2].plot(t_ms, bed_arr[:, 9] * 1.0e6, "C5", lw=1.2, label="rot")
ax[3].plot(t_ms, -bed_arr[:, 6], "C6", lw=1.2)
ax[0].set_ylabel("Particles")
ax[1].set_ylabel("Flow (mg/s)")
ax[2].set_ylabel(r"Energy ($\mu$J)")
ax[3].set_ylabel(r"$-v_z$ (m/s)")
ax[3].set_xlabel("Time (ms)")
ax[0].legend(frameon=False, fontsize=9, loc="upper right")
ax[2].legend(frameon=False, fontsize=9, loc="upper right")
for a in ax:
    finish_axis(a)
fig.subplots_adjust(left=0.17, right=0.97, bottom=0.09, top=0.98, hspace=0.12)
fig.savefig(figdir / "flow_metrics.pdf", bbox_inches="tight")
fig.savefig(figdir / "flow_metrics.png", dpi=300, bbox_inches="tight")

contact_arr = np.asarray(contact_time_rows) if contact_time_rows else np.zeros((0, len(contact_header)))
if len(contact_arr):
    fig, ax = plt.subplots(3, 1, figsize=(5.6, 5.8), sharex=True)
    ct_ms = contact_arr[:, 1] * 1.0e3
    ax[0].plot(ct_ms, contact_arr[:, 3], "C0", lw=1.0, label="grain-grain")
    ax[0].plot(ct_ms, contact_arr[:, 4], "C1", lw=1.0, label="grain-tube")
    ax[1].plot(ct_ms, contact_arr[:, 5] * 1.0e3, "C2", lw=1.0, label="p50")
    ax[1].plot(ct_ms, contact_arr[:, 6] * 1.0e3, "C3", lw=1.0, label="p90")
    ax[2].plot(ct_ms, contact_arr[:, 8], "C4", lw=1.0, label="p50")
    ax[2].plot(ct_ms, contact_arr[:, 9], "C5", lw=1.0, label="p90")
    ax[0].set_ylabel("Contacts")
    ax[1].set_ylabel("Force (mN)")
    ax[2].set_ylabel(r"$d_{actual}$ ($\mu$m)")
    ax[2].set_xlabel("Time (ms)")
    ax[0].legend(frameon=False, fontsize=9, loc="upper right")
    ax[1].legend(frameon=False, fontsize=9, loc="upper right")
    ax[2].legend(frameon=False, fontsize=9, loc="upper right")
    for a in ax:
        finish_axis(a)
    fig.subplots_adjust(left=0.17, right=0.97, bottom=0.10, top=0.98, hspace=0.12)
    fig.savefig(figdir / "contact_metrics.pdf", bbox_inches="tight")
    fig.savefig(figdir / "contact_metrics.png", dpi=300, bbox_inches="tight")

if len(grain_arr):
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 3.0))
    if d_force_samples:
        bins_d = np.geomspace(max(1.0e-4, np.nanmin(d_force_samples)), max(1.0e-3, np.nanmax(d_force_samples + d_vel_samples)), 50)
        ax[0].hist(d_force_samples, bins=bins_d, histtype="step", lw=1.4, label="force")
        if d_vel_samples:
            ax[0].hist(d_vel_samples, bins=bins_d, histtype="step", lw=1.4, label="velocity")
        ax[0].set_xscale("log")
    ax[0].set_xlabel(r"$d_{actual}$ ($\mu$m)")
    ax[0].set_ylabel("Contacts")
    ax[0].legend(frameon=False, fontsize=9)
    ax[1].hist(grain_arr[:, 3] * 1.0e3, bins=50, histtype="stepfilled", alpha=0.35, color="C3")
    ax[1].set_xlabel("Residence time (ms)")
    ax[1].set_ylabel("Grains")
    for a in ax:
        finish_axis(a)
    fig.subplots_adjust(left=0.11, right=0.98, bottom=0.24, top=0.95, wspace=0.38)
    fig.savefig(figdir / "grain_distributions.pdf", bbox_inches="tight")
    fig.savefig(figdir / "grain_distributions.png", dpi=300, bbox_inches="tight")

if viz_dump.exists() and len(grain_arr):
    final = None
    for step, cols, data in read_dump(viz_dump):
        final = (step, cols, data)
    if final is not None:
        step, cols, data = final
        idx = {name: i for i, name in enumerate(cols)}
        mobile = data[data[:, idx["type"]] == 1]
        tubes = data[data[:, idx["type"]] == 2]
        maxd = {int(row[0]): row[11] for row in grain_arr}
        colors = np.array([maxd.get(int(v), 0.0) for v in mobile[:, idx["id"]]])
        fig, ax = plt.subplots(figsize=(4.8, 6.0))
        sc = ax.scatter(mobile[:, idx["x"]] * 1.0e3, mobile[:, idx["z"]] * 1.0e3, s=4, c=colors, cmap="viridis", lw=0, alpha=0.65)
        ax.scatter(tubes[:, idx["x"]] * 1.0e3, tubes[:, idx["z"]] * 1.0e3, s=5, c="0.15", alpha=0.45, lw=0)
        cb = fig.colorbar(sc, ax=ax, pad=0.02)
        cb.set_label(r"Max $d_{actual}$ ($\mu$m)")
        ax.set_xlim(0, lx * 1.0e3)
        ax.set_ylim(0, lz * 1.0e3)
        ax.set_aspect("equal")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("z (mm)")
        finish_axis(ax)
        fig.subplots_adjust(left=0.18, right=0.92, bottom=0.08, top=0.98)
        fig.savefig(figdir / "snapshot_dactual.pdf", bbox_inches="tight")
        fig.savefig(figdir / "snapshot_dactual.png", dpi=300, bbox_inches="tight")

print((root / "summary.csv").read_text())
