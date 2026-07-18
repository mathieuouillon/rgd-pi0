"""Kinematic-coverage figures for the note, reconstructed from Stage A slims.

The slim carries photons and the @DIS electron but no reconstructed pi0, so this
module repeats Stage B's greedy gamma-gamma pairing and SIDIS kinematics in
Python to draw:

* ``kinematics_pi0.pdf``     -- the 2x2 coverage: M(gg), z, pT^2, phi_h,
* ``phi_h_acceptance.pdf``   -- the phi_h acceptance for one target,
* ``binning_grid_factorized.pdf`` -- the committed edges over the pooled data.

phi_h follows the Trento convention and is reported in DEGREES over
``[-180, 180]`` -- the convention of ``core/Kinematics.hpp``. (A prior figure
plotted it in radians, which is the bug this replaces.) Matplotlib only.

    python -m pi0.plots_kinematics --slimdir /volatile/.../rgd-pi0 \
        --grids config/binning --outdir note/figures [--max-events N]
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

TARGETS = [("LD2", "LD2"), ("CxC", "C×C"), ("Cu", "Cu"), ("Sn", "Sn")]
COLORS = {"LD2": "#1f77b4", "CxC": "#17a2b8", "Cu": "#d62728", "Sn": "#2ca02c"}

# --- cut values (config/cuts.json); passed in, defaulted here for clarity ---
MPI0 = 0.1349768
MASS_WINDOW = 0.2
MIN_MGG = 0.001
E_GAMMA_MIN_DEG = 8.0
OA_A, OA_B, OA_OFF = 17.561, 0.756, 1.0  # theta_min(p) = A exp(-B p) + offset
EBEAM = 10.53


def _opening_min_deg(p: float) -> float:
    return OA_A * math.exp(-OA_B * p) + OA_OFF


def reconstruct(ev: dict, cuts: dict) -> dict[str, np.ndarray]:
    """Greedy-exclusive gamma-gamma pairing + SIDIS kinematics, per event.

    ``ev`` is the slim's ``events`` arrays (numpy, jagged branches as object
    arrays). Returns arrays of the per-pi0 M(gg), z, pT^2 and phi_h (deg).
    """
    p = cuts.get("pairing", {})
    win = float(p.get("mass_window_gev", MASS_WINDOW))
    mmin = float(p.get("min_mgg_gev", MIN_MGG))
    ega_min = float(p.get("e_gamma_min_angle_deg", E_GAMMA_MIN_DEG))
    oa = p.get("opening_angle", {})
    a_deg = float(oa.get("a_deg", OA_A)); b_inv = float(oa.get("b_inv_gev", OA_B))
    off = float(oa.get("offset_deg", OA_OFF))
    zmin = float(p.get("z_min", 0.0)); zmax = float(p.get("z_max", 1.0))
    ebeam = float(cuts.get("beam", {}).get("energy_gev", EBEAM))

    def theta_min(pp: float) -> float:
        return a_deg * math.exp(-b_inv * pp) + off

    mgg_o, z_o, pt2_o, phi_o = [], [], [], []
    q2 = ev["q2"]; nu = ev["nu"]
    ex, ey, ez = ev["ex"], ev["ey"], ev["ez"]
    gpx, gpy, gpz, gega = ev["gpx"], ev["gpy"], ev["gpz"], ev["g_e_gamma_deg"]
    n = len(nu)
    for i in range(n):
        # e-gamma angle cut, applied per photon before pairing
        keep = np.asarray(gega[i], dtype=float) > ega_min
        if keep.sum() < 2:
            continue
        px = np.asarray(gpx[i], dtype=float)[keep]
        py = np.asarray(gpy[i], dtype=float)[keep]
        pz = np.asarray(gpz[i], dtype=float)[keep]
        E = np.sqrt(px * px + py * py + pz * pz)
        # virtual-photon frame from this event's electron
        qx, qy, qz = -ex[i], -ey[i], ebeam - ez[i]
        qn = math.sqrt(qx * qx + qy * qy + qz * qz)
        qhx, qhy, qhz = qx / qn, qy / qn, qz / qn
        # y_hat = (q x k)/|.|, k = (0,0,Ebeam);  x_hat = y_hat x q_hat
        yx, yy, yz = qy * ebeam, -qx * ebeam, 0.0
        yn = math.hypot(yx, yy) or 1.0
        yx, yy, yz = yx / yn, yy / yn, 0.0
        xhx = yy * qhz - yz * qhy
        xhy = yz * qhx - yx * qhz
        xhz = yx * qhy - yy * qhx

        m = len(E)
        cand = []
        for a in range(m):
            for b in range(a + 1, m):
                cos12 = (px[a] * px[b] + py[a] * py[b] + pz[a] * pz[b]) / (E[a] * E[b])
                mass = math.sqrt(max(0.0, 2.0 * E[a] * E[b] * (1.0 - cos12)))
                if mass < mmin or abs(mass - MPI0) >= win:
                    continue
                ppx, ppy, ppz = px[a] + px[b], py[a] + py[b], pz[a] + pz[b]
                pmag = math.sqrt(ppx * ppx + ppy * ppy + ppz * ppz)
                theta12 = math.degrees(math.acos(max(-1.0, min(1.0, cos12))))
                if theta12 <= theta_min(pmag):
                    continue
                cand.append((abs(mass - MPI0), a, b, mass, ppx, ppy, ppz, E[a] + E[b]))
        cand.sort(key=lambda c: c[0])
        used: set[int] = set()
        for _, a, b, mass, ppx, ppy, ppz, epi in cand:
            if a in used or b in used:
                continue
            z = epi / nu[i]
            if not (zmin < z < zmax):
                continue
            used.add(a); used.add(b)
            pdotq = ppx * qhx + ppy * qhy + ppz * qhz
            perpx, perpy, perpz = ppx - pdotq * qhx, ppy - pdotq * qhy, ppz - pdotq * qhz
            pt2 = perpx * perpx + perpy * perpy + perpz * perpz
            phi = math.degrees(math.atan2(perpx * yx + perpy * yy + perpz * yz,
                                          perpx * xhx + perpy * xhy + perpz * xhz))
            mgg_o.append(mass); z_o.append(z); pt2_o.append(pt2); phi_o.append(phi)
    return {"mgg": np.array(mgg_o), "z": np.array(z_o),
            "pt2": np.array(pt2_o), "phih": np.array(phi_o)}


def _load_slims(slimdir: Path, target: str, max_events: int | None) -> dict:
    import uproot
    files = sorted((slimdir / target).glob("*/slim_*.root"))
    cols = ["q2", "xb", "nu", "ex", "ey", "ez", "ee", "gpx", "gpy", "gpz", "g_e_gamma_deg"]
    acc = {c: [] for c in cols}
    seen = 0
    for f in files:
        a = uproot.open(f)["events"].arrays(cols, library="np")
        acc = {c: acc[c] + [a[c]] for c in cols}
        seen += len(a["nu"])
        if max_events and seen >= max_events:
            break
    return {c: np.concatenate(acc[c]) if acc[c] else np.array([]) for c in cols}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="pi0 kinematic-coverage figures.")
    p.add_argument("--slimdir", type=Path, required=True,
                   help="dir holding <target>/<run>/slim_*.root")
    p.add_argument("--grids", type=Path, required=True, help="config/binning directory")
    p.add_argument("--config", type=Path, default=Path("config"), help="config directory")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--max-events", type=int, default=400000,
                   help="cap events read per target (coverage plot; default 400k)")
    args = p.parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    cuts = json.load(open(args.config / "cuts.json"))

    import matplotlib as mpl
    mpl.use("Agg")
    mpl.rcParams.update({"savefig.dpi": 120, "savefig.bbox": "tight", "font.size": 9.5})
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    recon = {}
    events = {}  # per-event (x_B, Q^2) for the Grid A density overlay
    for t, _ in TARGETS:
        ev = _load_slims(args.slimdir, t, args.max_events)
        recon[t] = reconstruct(ev, cuts)
        events[t] = {"q2": np.asarray(ev["q2"], float),
                     "xb": np.asarray(ev["xb"], float)}

    def step(ax, data, lo, hi, nb, t, lbl):
        h, edges = np.histogram(data[(data >= lo) & (data <= hi)], bins=nb, range=(lo, hi),
                                density=True)
        ax.step(edges, np.r_[h, h[-1]], where="post", color=COLORS[t], lw=1.1, label=lbl)

    # --- Figure 1: 2x2 coverage --------------------------------------------
    fig, ax = plt.subplots(2, 2, figsize=(8.2, 5.6))
    fig.suptitle(r"$\pi^0$ kinematics --- unit-area normalized", y=0.98)
    for t, lbl in TARGETS:
        r = recon[t]
        step(ax[0, 0], r["mgg"], 0.0, 0.5, 100, t, lbl)
        step(ax[0, 1], r["z"], 0.0, 1.0, 60, t, lbl)
        step(ax[1, 0], r["pt2"], 0.0, 2.0, 80, t, lbl)
        step(ax[1, 1], r["phih"], -180.0, 180.0, 36, t, lbl)
    ax[0, 0].set(title=r"$M(\gamma\gamma)$", xlabel=r"$M(\gamma\gamma)$ [GeV]")
    ax[0, 1].set(title=r"$z$", xlabel=r"$z$")
    ax[1, 0].set(title=r"$p_T^2$", xlabel=r"$p_T^2$ [GeV$^2$]")
    ax[1, 1].set(title=r"$\phi_h$", xlabel=r"$\phi_h$ [deg]")
    for a in ax.flat:
        a.set_ylabel("Normalized density")
        a.grid(alpha=0.2)
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.legend(h, l, ncol=4, loc="lower center", frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    fig.savefig(args.outdir / "kinematics_pi0.pdf")
    plt.close(fig)

    # --- phi_h acceptance (LD2) --------------------------------------------
    fig, a = plt.subplots(figsize=(6.4, 3.6))
    ph = recon["LD2"]["phih"]
    a.hist(ph, bins=36, range=(-180, 180), histtype="step", color=COLORS["LD2"], lw=1.2)
    a.set(title=r"$\phi_h$ yield, LD$_2$ (raw acceptance)",
          xlabel=r"$\phi_h$ [deg]", ylabel=r"$N^+ + N^-$", xlim=(-180, 180))
    a.grid(alpha=0.2)
    fig.savefig(args.outdir / "phi_h_acceptance.pdf")
    plt.close(fig)

    # --- factorized grid, drawn over the pooled data it partitions ----------
    # The density overlay is what makes the figure show its own caption: the
    # equal-statistics packing (bins narrow where data is dense) and the empty
    # corner cells left by the Q^2-x_B correlation. Pooled over all four
    # targets, matching how make_grid fit the edges -- Grid A per event
    # (x_B, Q^2), Grid B per pi0 (p_T^2, z).
    ga = json.load(open(args.grids / "grid_A_q2_xb.json"))
    gb = json.load(open(args.grids / "grid_B_z_pt2.json"))
    (q2e, xbe) = (ga["axes"][0]["edges"], ga["axes"][1]["edges"])
    (ze, pt2e) = (gb["axes"][0]["edges"], gb["axes"][1]["edges"])
    xb_all = np.concatenate([events[t]["xb"] for t, _ in TARGETS])
    q2_all = np.concatenate([events[t]["q2"] for t, _ in TARGETS])
    pt2_all = np.concatenate([recon[t]["pt2"] for t, _ in TARGETS])
    z_all = np.concatenate([recon[t]["z"] for t, _ in TARGETS])

    def _grid(ax, xd, yd, xe, ye, xhi, title, xlabel, ylabel):
        good = np.isfinite(xd) & np.isfinite(yd)
        if good.any():
            ax.hist2d(xd[good], yd[good], bins=(140, 140),
                      range=[(xe[0], xhi), (ye[0], ye[-1])],
                      cmap="Blues", cmin=1, norm=LogNorm())
        for x in xe:
            ax.axvline(x, color="0.25", lw=0.6)
        for y in ye:
            ax.axhline(y, color="0.25", lw=0.6)
        ax.set(title=title, xlabel=xlabel, ylabel=ylabel,
               xlim=(xe[0], xhi), ylim=(ye[0], ye[-1]))

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(9.0, 4.0))
    _grid(axa, xb_all, q2_all, xbe, q2e, xbe[-1],
          r"Grid A  $(x_B, Q^2)$", r"$x_B$", r"$Q^2$ [GeV$^2$]")
    _grid(axb, pt2_all, z_all, pt2e, ze, min(pt2e[-1], 1.5),
          r"Grid B  $(p_T^2, z)$", r"$p_T^2$ [GeV$^2$]", r"$z$")
    fig.tight_layout()
    fig.savefig(args.outdir / "binning_grid_factorized.pdf")
    plt.close(fig)

    print(f"wrote kinematics_pi0.pdf, phi_h_acceptance.pdf, binning_grid_factorized.pdf "
          f"to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
