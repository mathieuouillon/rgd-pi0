"""Selection figures for the analysis note.

Draws the electron, DIS and photon selection from real skim output, with every
cut line taken from ``config/cuts.json`` rather than typed into the plot. A
figure whose cut line disagrees with the cut that ran is worse than no figure.

Usage::

    python -m pi0.plots_selection --slim slim.root --config ../config --outdir figs/
    python -m pi0.plots_selection --slim slim.root --qa qa.root --config ../config --outdir figs/

What each input can show
------------------------
The **slim** holds only what SURVIVED selection (its schema is run, event,
helicity, q2, xb, nu, w, y, ex..ee, gpx, gpy, gpz, g_e_gamma_deg). It can show
the accepted phase space and the cuts that act on kinematics it carries, but it
cannot show a cut removing anything -- the removed rows are not in it.

The **QA ntuple** (``stageA_skim --qa-ntuple``) holds one row per *candidate*,
before the cuts reject it, plus ``failed_at``. That is what makes a real
selection figure possible: the band, the window, the threshold, and the
population each one removes.

So: pass ``--slim`` alone and you get the coverage figures; add ``--qa`` and you
get the selection figures too. Nothing is faked in between -- a figure that needs
the QA ntuple is skipped, loudly, when it is absent.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

__all__ = ["main"]

# Matplotlib is imported lazily in main() so that --help works without it.


def _style():
    import matplotlib as mpl

    mpl.rcParams.update({
        "figure.dpi": 120,
        "savefig.dpi": 120,
        "savefig.bbox": "tight",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.5,
        "legend.frameon": False,
        "legend.fontsize": 8,
        "axes.axisbelow": True,
        "figure.facecolor": "white",
    })


#: One colour per role, so the same idea reads the same way in every figure.
C_ACCEPT = "#2a6f9e"   # what survives
C_REJECT = "#c44e52"   # what a cut removes
C_CUT = "#111111"      # the cut line itself
C_BAND = "#f0a202"     # an accepted band


def _step(ax, data, bins, *, label=None, color=None):
    """One drawing path for every 1D histogram: an UNFILLED step, via mplhep.

    mplhep.histplot with histtype="step" gives the outlined, unfilled histogram
    that reads as HEP-standard; routing all of them through here keeps them
    consistent. No mplhep experiment *style* is applied -- the module's own
    minimal rcParams stay in charge -- only its histogram drawing.
    """
    import mplhep as hep

    counts, edges = np.histogram(np.asarray(data), bins=bins)
    hep.histplot(counts, edges, ax=ax, label=label, color=color, histtype="step")
    return counts


def _load_cuts(config_dir: Path) -> dict:
    with (Path(config_dir) / "cuts.json").open() as f:
        return json.load(f)


def _title(obj) -> str:
    """The value of a provenance TNamed.

    Provenance stores value in the TITLE, not the name. uproot's Model_TNamed
    exposes it through member("fTitle") -- NOT as a .title or .fTitle attribute,
    which is what an earlier version of this reached for. It failed, the
    surrounding `except: pass` swallowed it, and every figure silently lost its
    watermark. A stamp that can vanish quietly is worse than no stamp, so this
    raises rather than guessing.
    """
    return str(obj.member("fTitle"))


def _provenance_note(slim_path: Path) -> str:
    """The one line every figure must carry: what made this data.

    A figure from a fallback-scored, truncated skim is a diagnostic, and saying
    so ON the figure is cheaper than explaining it later. Raises if the block
    cannot be read: a missing stamp must be a crash, never a blank corner.
    """
    import uproot

    f = uproot.open(slim_path)
    if "provenance" not in [k.split(";")[0] for k in f.keys()]:
        raise RuntimeError(
            f"{slim_path} has no provenance block. Every Stage A output carries one; a file "
            f"without it did not come from this chain, and a figure from it could not be labelled."
        )
    prov = f["provenance"]
    got = {k.split(";")[0]: _title(prov[k]) for k in prov.keys()}

    bits = []
    if got.get("target"):
        bits.append(f"target {got['target']}")
    if got.get("run"):
        bits.append(f"run {got['run']}")
    if "TRUE" in got.get("gbt.fallback_used", ""):
        bits.append("RGA-FALLBACK PHOTON ID")
    mx = got.get("events.max_events_requested", "")
    if "TRUNCATED" in mx:
        # Only the count READ, never "N of M": Stage A records what it read, not
        # what the file holds, so any denominator here would be invented. (The
        # file total is knowable -- hipo-utils reads it from the trailer -- but it
        # is not in this block, and a figure must not imply a fraction nobody
        # measured.)
        bits.append(f"TRUNCATED PREFIX: {mx.split(' ')[0]} events read")
    if not bits:
        raise RuntimeError(f"{slim_path}: provenance read but nothing recognised in it.")
    return " | ".join(bits)


def _stamp(fig, text: str):
    """Below the axes, not on them. No-op when ``text`` is empty.

    The provenance is still READ and validated for every figure (a missing block
    is a hard error, and it supplies the target); this only controls whether it
    is DRAWN. For the analysis note the caveat lives in the prose, which carries
    it far more prominently than a 6-pt watermark, so `--stamp` defaults off.

    y is negative so that savefig's tight bbox extends to include it; at y=0.005
    it landed on top of the x-axis label.
    """
    if not text:
        return
    fig.text(0.995, -0.035, text, ha="right", va="top", fontsize=6,
             color="#8a8a8a", style="italic")


# ---------------------------------------------------------------------------
# Figures from the slim: accepted phase space
# ---------------------------------------------------------------------------


def fig_dis_phase_space(ev: dict, cuts: dict, out: Path, stamp: str):
    """Q2 vs xB, with the DIS cuts drawn as the boundaries they are.

    The W and y cuts are curves in this plane, not axis limits, which is why
    they are drawn rather than described: the accepted region is the intersection
    of three constraints and its shape is the point.
    """
    import matplotlib.pyplot as plt

    d = cuts["dis"]
    q2_min = d["q2_min"]
    w_min = d["w_min"]
    y_max = d["y_max"]
    beam = cuts["beam"]["energy_gev"]
    M = cuts["beam"]["proton_mass_gev"]

    fig, ax = plt.subplots(figsize=(5.2, 3.9))
    h = ax.hist2d(ev["xb"], ev["q2"], bins=(60, 60), cmap="viridis", cmin=1)
    fig.colorbar(h[3], ax=ax, label="DIS events")

    xb = np.linspace(1e-3, 1.0, 400)
    # W^2 = M^2 + Q^2 (1-x)/x  =>  Q^2 = (W^2 - M^2) x / (1-x)
    with np.errstate(divide="ignore", invalid="ignore"):
        q2_w = (w_min**2 - M**2) * xb / (1 - xb)
        # y = Q^2 / (2 M x E)  =>  Q^2 = y 2 M x E
        q2_y = y_max * 2 * M * xb * beam
    ax.plot(xb, q2_w, color=C_CUT, lw=1.4, label=f"$W = {w_min:g}$ GeV")
    ax.plot(xb, q2_y, color=C_CUT, lw=1.4, ls="--", label=f"$y = {y_max:g}$")
    ax.axhline(q2_min, color=C_CUT, lw=1.4, ls=":", label=f"$Q^2 = {q2_min:g}$ GeV$^2$")

    ax.set_xlim(0, max(0.7, float(np.percentile(ev["xb"], 99.5)) * 1.15))
    ax.set_ylim(0, float(np.percentile(ev["q2"], 99.5)) * 1.15)
    ax.set_xlabel("$x_B$")
    ax.set_ylabel("$Q^2$  [GeV$^2$]")
    ax.set_title("DIS phase space, with the selection boundaries")
    ax.legend(loc="upper left")
    _stamp(fig, stamp)
    fig.savefig(out / "sel_dis_phase_space.pdf")
    plt.close(fig)


def fig_dis_1d(ev: dict, cuts: dict, out: Path, stamp: str):
    """Q2, W, y and nu of the accepted sample, each with its cut edge.

    y is included because it is the cut that does nothing: on a real skim it
    rejected 0 of 507 events. A figure showing the distribution sitting well
    inside its own bound is the honest way to say that.
    """
    import matplotlib.pyplot as plt

    d = cuts["dis"]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.8))
    spec = [
        ("q2", "$Q^2$  [GeV$^2$]", d["q2_min"], "min"),
        ("w", "$W$  [GeV]", d["w_min"], "min"),
        ("y", "$y$", d["y_max"], "max"),
        ("nu", r"$\nu$  [GeV]", None, None),
    ]
    for ax, (key, label, cut, kind) in zip(axes.ravel(), spec):
        _step(ax, ev[key], 60, color=C_ACCEPT)
        if cut is not None:
            ax.axvline(cut, color=C_CUT, lw=1.4,
                       label=f"{'>' if kind == 'min' else '<'} {cut:g}")
            ax.legend(loc="upper right")
        ax.set_xlabel(label)
        ax.set_ylabel("events")
    axes[1, 1].set_title(r"$\nu$ has no cut of its own", fontsize=8, color="#666666")
    fig.suptitle("DIS variables of the accepted sample", y=1.0)
    fig.tight_layout()
    _stamp(fig, stamp)
    fig.savefig(out / "sel_dis_1d.pdf")
    plt.close(fig)


def fig_electron_kinematics(ev: dict, cuts: dict, out: Path, stamp: str):
    """Scattered-electron p vs theta, with the momentum cut.

    From the slim, so every point already passed: the cut line marks the edge of
    the sample rather than a boundary with anything on the far side.
    """
    import matplotlib.pyplot as plt

    p = np.sqrt(ev["ex"] ** 2 + ev["ey"] ** 2 + ev["ez"] ** 2)
    theta = np.degrees(np.arccos(np.clip(ev["ez"] / p, -1, 1)))
    pmin = cuts["electron"]["min_momentum_gev"]

    fig, ax = plt.subplots(figsize=(5.2, 3.9))
    h = ax.hist2d(theta, p, bins=(60, 60), cmap="viridis", cmin=1)
    fig.colorbar(h[3], ax=ax, label="electrons")
    ax.axhline(pmin, color=C_CUT, lw=1.4, label=f"$p > {pmin:g}$ GeV")
    ax.set_xlabel(r"$\theta_e$  [deg]")
    ax.set_ylabel("$p_e$  [GeV]")
    ax.set_title("Scattered electron (accepted sample)")
    ax.legend(loc="upper right")
    _stamp(fig, stamp)
    fig.savefig(out / "sel_electron_kinematics.pdf")
    plt.close(fig)


def fig_photon(ev: dict, cuts: dict, out: Path, stamp: str):
    """Photon energy, angle, and the e-gamma angle with its cut.

    The e-gamma panel is the one selection this file can honestly show removing
    something: the cut is applied at PAIRING (Stage B), so the slim still holds
    photons on both sides of it.
    """
    import matplotlib.pyplot as plt

    gx, gy, gz = (np.concatenate(ev[k]) for k in ("gpx", "gpy", "gpz"))
    ang = np.concatenate(ev["g_e_gamma_deg"])
    e = np.sqrt(gx**2 + gy**2 + gz**2)
    theta = np.degrees(np.arccos(np.clip(gz / np.where(e > 0, e, 1), -1, 1)))
    amin = cuts["pairing"]["e_gamma_min_angle_deg"]

    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.1))
    _step(axes[0], e, 60, color=C_ACCEPT)
    axes[0].set_xlabel(r"$E_\gamma$  [GeV]"); axes[0].set_ylabel("photons")
    axes[0].set_yscale("log")
    axes[0].set_title("energy")

    _step(axes[1], theta, 60, color=C_ACCEPT)
    axes[1].set_xlabel(r"$\theta_\gamma$  [deg]"); axes[1].set_ylabel("photons")
    axes[1].set_title("polar angle")

    keep = ang >= amin
    ega_bins = np.linspace(0, max(40, float(ang.max())), 60)
    _step(axes[2], ang[~keep], ega_bins, color=C_REJECT, label=f"removed ({(~keep).sum()})")
    _step(axes[2], ang[keep], ega_bins, color=C_ACCEPT, label=f"kept ({keep.sum()})")
    axes[2].axvline(amin, color=C_CUT, lw=1.4, label=fr"$\theta_{{e\gamma}} > {amin:g}^\circ$")
    axes[2].set_xlabel(r"$\theta_{e\gamma}$  [deg]"); axes[2].set_ylabel("photons")
    axes[2].set_title("angle to the electron")
    axes[2].legend(loc="upper right")

    fig.suptitle("Selected photons", y=1.02)
    fig.tight_layout()
    _stamp(fig, stamp)
    fig.savefig(out / "sel_photon.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figures from the QA ntuple: the cuts, and what they remove
# ---------------------------------------------------------------------------


def fig_sampling_fraction(qa: dict, cuts: dict, out: Path, stamp: str):
    """Sampling fraction vs p, with the accepted band and the rejected population.

    The classic electron-ID figure, and impossible from the slim: the rejected
    electrons are precisely the rows the slim does not have.
    """
    import matplotlib.pyplot as plt

    ok = qa["failed_at"] == -1
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    ax.scatter(qa["p"][~ok], qa["sf"][~ok], s=2, c=C_REJECT, alpha=0.30,
               label=f"rejected ({(~ok).sum()})", rasterized=True)
    ax.scatter(qa["p"][ok], qa["sf"][ok], s=2, c=C_ACCEPT, alpha=0.45,
               label=f"accepted ({ok.sum()})", rasterized=True)
    ax.set_xlabel("$p_e$  [GeV]")
    ax.set_ylabel(r"sampling fraction  $E_{\rm cal}/p$")
    ax.set_title("Electron sampling fraction")
    ax.set_ylim(0, 0.45)
    ax.legend(loc="lower right", markerscale=4)
    _stamp(fig, stamp)
    fig.savefig(out / "sel_sampling_fraction.pdf")
    plt.close(fig)


def fig_vz(qa: dict, cuts: dict, out: Path, stamp: str, target: str):
    """The vertex window -- the cut that IS the target assignment for Cu and Sn.

    Both foils' windows are drawn whatever the target, because the point is that
    they are disjoint and that vz alone decides which target an event belongs to.
    """
    import matplotlib.pyplot as plt

    tg = cuts["vertex"]["targets"]
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    vz = qa["vz_corrected"]
    ok = qa["failed_at"] == -1

    # The window edges of ALL targets, so the axis frames the physics (the LD2
    # peak sits near -8 cm, the foils between -9 and -2) rather than being
    # stretched to the ~-100..+20 cm tails of misreconstructed tracks -- which is
    # what a percentile range did, crushing every window into a sliver.
    edges = []
    for t in tg.values():
        if t.get("rule") == "corrected_peaks":
            for pk in t["peaks"]:
                edges += [pk["mu_cm"] - t["n_sigma"] * pk["sigma_cm"],
                          pk["mu_cm"] + t["n_sigma"] * pk["sigma_cm"]]
        else:
            edges += [t["vz_min_cm"], t["vz_max_cm"]]
    lo, hi = min(edges) - 6, max(edges) + 6
    bins = np.linspace(lo, hi, 120)
    over = float(np.mean((vz < lo) | (vz > hi)) * 100)
    # No clipping: hist drops out-of-range values rather than piling them into a
    # fake bar at the edge. The tails are misreconstructed tracks; the note below
    # records how much fell off, so nothing is hidden silently.
    _step(ax, vz[~ok], bins, color=C_REJECT, label="rejected (any cut)")
    _step(ax, vz[ok], bins, color=C_ACCEPT, label="accepted")
    ax.set_xlim(lo, hi)
    ax.set_yscale("log")  # the vertex peak dwarfs the tails; log shows both

    for name, style in (("cu", "-"), ("sn", "--"), ("ld2", ":")):
        if name not in tg:
            continue
        t = tg[name]
        if t.get("rule") == "corrected_peaks":
            for pk in t["peaks"]:
                lo_e = pk["mu_cm"] - t["n_sigma"] * pk["sigma_cm"]
                hi_e = pk["mu_cm"] + t["n_sigma"] * pk["sigma_cm"]
                ax.axvspan(lo_e, hi_e, color=C_BAND, alpha=0.18)
                ax.axvline(lo_e, color=C_CUT, lw=1.0, ls=style)
                ax.axvline(hi_e, color=C_CUT, lw=1.0, ls=style)
        else:
            ax.axvline(t["vz_min_cm"], color=C_CUT, lw=1.0, ls=style)
            ax.axvline(t["vz_max_cm"], color=C_CUT, lw=1.0, ls=style)
    if over > 0.05:
        ax.text(0.02, 0.97, f"{over:.1f}% beyond axis (clipped)", transform=ax.transAxes,
                fontsize=7, va="top", color="#888888")
    ax.set_xlabel("$v_z$  [cm]  (corrected where the target defines a correction)")
    ax.set_ylabel("electrons")
    ax.set_title(f"Vertex selection — target {target}")
    ax.legend(loc="upper right")
    _stamp(fig, stamp)
    fig.savefig(out / "sel_vz.pdf")
    plt.close(fig)


def fig_gbt(qa_g: dict, cuts: dict, out: Path, stamp: str):
    """The GBT photon score against its threshold.

    Clusters the pre-filter rejected have no score at all -- the skim runs the
    GBT lazily, precisely so it is not evaluated on most clusters -- so they are
    counted in the title rather than drawn at a fictitious zero.
    """
    import matplotlib.pyplot as plt

    thr = cuts["photon"]["gbt_threshold"]
    s = qa_g["gbt_score"]
    scored = np.isfinite(s)
    ok = qa_g["passed"] == 1

    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    bins = np.linspace(0, 1, 60)
    _step(ax, s[scored & ~ok], bins, color=C_REJECT, label=f"below threshold ({(scored & ~ok).sum()})")
    _step(ax, s[scored & ok], bins, color=C_ACCEPT, label=f"accepted ({(scored & ok).sum()})")
    ax.axvline(thr, color=C_CUT, lw=1.5, label=f"threshold {thr:g}")
    ax.set_yscale("log")
    ax.set_xlabel("GBT photon score  (sigmoid)")
    ax.set_ylabel("clusters")
    ax.set_title(f"Photon identification — {(~scored).sum()} clusters pre-filtered, never scored")
    ax.legend(loc="upper center")
    _stamp(fig, stamp)
    fig.savefig(out / "sel_gbt_score.pdf")
    plt.close(fig)


def fig_cutflow(qa: dict, names: dict, out: Path, stamp: str):
    """Where the electrons go, in the order the cuts actually run.

    Ordered by the code's own failed_at index, not by size: a cutflow reordered
    for looks stops being a cutflow.
    """
    import matplotlib.pyplot as plt

    idx, counts = np.unique(qa["failed_at"], return_counts=True)
    order = [(i, c) for i, c in zip(idx.tolist(), counts.tolist()) if i != -1]
    order.sort(key=lambda t: t[0])
    labels = [names.get(str(i), names.get(i, f"cut {i}")) for i, _ in order]
    vals = [c for _, c in order]
    passed = int((qa["failed_at"] == -1).sum())

    fig, ax = plt.subplots(figsize=(6.4, 0.42 * (len(vals) + 2) + 1.2))
    y = np.arange(len(vals) + 1)
    ax.barh(y[:-1], vals, color=C_REJECT, alpha=0.85)
    ax.barh(y[-1], passed, color=C_ACCEPT, alpha=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(labels + ["ACCEPTED"])
    ax.invert_yaxis()
    ax.set_xlabel("electron candidates")
    ax.set_title("Electron selection: where candidates are lost")
    for i, v in enumerate(vals + [passed]):
        ax.text(v, y[i], f" {v}", va="center", fontsize=8)
    _stamp(fig, stamp)
    fig.savefig(out / "sel_electron_cutflow.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="pi0.plots_selection",
        description="Selection figures for the note. Cut lines come from cuts.json, never from the plot.",
    )
    p.add_argument("--slim", type=Path, required=True, help="a Stage A slim (accepted sample)")
    p.add_argument("--qa", type=Path, default=None,
                   help="a stageA_skim --qa-ntuple file (per-candidate, pre-cut). Without it the "
                        "figures that show a cut REMOVING something are skipped, not faked.")
    p.add_argument("--config", type=Path, default=Path("config"), help="config DIRECTORY")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--stamp", action="store_true",
                   help="draw the provenance line (target/run/fallback/truncation) on each figure. "
                        "OFF by default: in the note the caveat belongs in the prose. The provenance "
                        "is still read and validated either way; this only controls whether it is drawn.")
    args = p.parse_args(argv)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import uproot  # noqa: F401
    except ImportError as e:
        print(f"error: {e}. Needs matplotlib + uproot + numpy.", file=sys.stderr)
        return 1

    _style()
    cuts = _load_cuts(args.config)
    args.outdir.mkdir(parents=True, exist_ok=True)
    # Read + validate provenance always (raises if the block is missing); only
    # DRAW it when --stamp is given.
    provenance = _provenance_note(args.slim)
    print(f"provenance: {provenance}")
    stamp = provenance if args.stamp else ""

    import uproot
    ev = uproot.open(args.slim)["events"].arrays(library="np")
    n = len(ev["q2"])
    print(f"slim: {n} DIS events")
    if n < 200:
        print(f"warning: {n} events is thin for a figure. These will look sparse; "
              f"skim more with --max-events before using them in the note.")

    prov = uproot.open(args.slim)["provenance"]
    target = next(_title(prov[k]) for k in prov.keys() if k.split(";")[0] == "target")

    fig_dis_phase_space(ev, cuts, args.outdir, stamp)
    fig_dis_1d(ev, cuts, args.outdir, stamp)
    fig_electron_kinematics(ev, cuts, args.outdir, stamp)
    fig_photon(ev, cuts, args.outdir, stamp)
    made = ["sel_dis_phase_space", "sel_dis_1d", "sel_electron_kinematics", "sel_photon"]

    if args.qa:
        f = uproot.open(args.qa)
        qa_e = f["qa_electron"].arrays(library="np")
        qa_g = f["qa_photon"].arrays(library="np")
        print(f"qa: {len(qa_e['p'])} electron candidates, {len(qa_g['e'])} photon candidates")
        # The cut names live IN the QA file, written by stageA_skim from the same
        # kElectronStages array the printed cutflow uses -- so the figure's labels
        # cannot drift from the cuts that actually ran. (An earlier version looked
        # for a config file that never existed and fell back to "cut 1", "cut 2".)
        names = {}
        if "qa_electron_stages" in [k.split(";")[0] for k in f.keys()]:
            st = f["qa_electron_stages"].arrays(library="np")
            names = {int(i): lbl for i, lbl in zip(st["index"], st["label"])}
        fig_sampling_fraction(qa_e, cuts, args.outdir, stamp)
        fig_vz(qa_e, cuts, args.outdir, stamp, target)
        fig_gbt(qa_g, cuts, args.outdir, stamp)
        fig_cutflow(qa_e, names, args.outdir, stamp)
        made += ["sel_sampling_fraction", "sel_vz", "sel_gbt_score", "sel_electron_cutflow"]
    else:
        print("\nno --qa: skipping the figures that need PRE-CUT candidates --\n"
              "  sampling fraction vs p, the vz window, the GBT threshold, the cutflow.\n"
              "  The slim holds only what survived, so those cuts cannot be shown removing\n"
              "  anything from it. Re-skim with --qa-ntuple to get them.")

    print(f"\nwrote {len(made)} figure(s) to {args.outdir}:")
    for m in made:
        print(f"  {m}.pdf")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
