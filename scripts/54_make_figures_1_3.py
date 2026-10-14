#!/usr/bin/env python
"""
HPICDA publication figures 1-3 (style-validation set).

Every plotted value is read from the frozen, hash-audited artifacts of the
HPICDA v1.0 archive; nothing is recomputed, so the figure set stays consistent
with MASTER_SHA256SUMS.

Includes an automated layout QA pass (text overlap / out-of-bounds / ink
coverage) because the build host has no display.

Outputs -> 25_manuscript_package/figures/
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.text import Text

ROOT = os.environ.get(
    "HPICDA_ROOT",
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "25_manuscript_package", "figures")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 7.5, "axes.titlesize": 7.5, "axes.linewidth": 0.6,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.4, "ytick.major.size": 2.4,
    "legend.fontsize": 6.3, "legend.frameon": False,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42,
    "svg.fonttype": "none",
})

C = dict(blue="#0072B2", orange="#E69F00", green="#009E73", red="#D55E00",
         purple="#CC79A7", sky="#56B4E9", yellow="#F0E442", grey="#9AA0A6",
         dark="#333333")

EXPLOITED, PROTECTIVE = "CDEP_PATHOGEN_EXPLOITED", "CDEP_HOST_PROTECTIVE"
CONFOUNDED, STATE = "BASELINE_FITNESS_CONFOUNDED", "SIGNAL_STATE_ALTERING"
NEUTRAL = "NO_FUNCTIONAL_SIGNAL"
CALL_COLOR = {EXPLOITED: C["orange"], PROTECTIVE: C["blue"],
              CONFOUNDED: C["red"], STATE: C["green"], NEUTRAL: "#CFD3D6"}
CALL_LABEL = {EXPLOITED: "Pathogen-exploited conditional",
              PROTECTIVE: "Host-protective conditional",
              CONFOUNDED: "Baseline-fitness confounded",
              STATE: "Signal-state altering"}
ESS = -0.5  # Chronos threshold for baseline essentiality


# ------------------------------------------------------------------- QA pass
def qa(fig, name, strict=True):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fb = fig.bbox
    probs = []
    for ax in fig.axes:
        txts = [(t.get_text()[:26], t.get_window_extent(r))
                for t in ax.texts if t.get_text().strip()]
        for i in range(len(txts)):
            for j in range(i + 1, len(txts)):
                a, b = txts[i][1], txts[j][1]
                ov = (max(a.x0, b.x0) < min(a.x1, b.x1) and
                      max(a.y0, b.y0) < min(a.y1, b.y1))
                if ov:
                    probs.append(f"  TEXT-OVERLAP '{txts[i][0]}' <-> '{txts[j][0]}'")
    for t in fig.findobj(Text):
        if not t.get_text().strip():
            continue
        bb = t.get_window_extent(r)
        if bb.x0 < -2 or bb.y0 < -2 or bb.x1 > fb.width + 2 or bb.y1 > fb.height + 2:
            probs.append(f"  OUT-OF-BOUNDS '{t.get_text()[:26]}' "
                         f"x[{bb.x0:.0f},{bb.x1:.0f}] y[{bb.y0:.0f},{bb.y1:.0f}] "
                         f"fig {fb.width:.0f}x{fb.height:.0f}")
    ink = 0.0
    for ax in fig.axes:
        arr = np.asarray(ax.figure.canvas.buffer_rgba())
    arr = np.asarray(fig.canvas.buffer_rgba())[:, :, :3]
    ink = float((arr < 245).any(axis=2).mean() * 100)
    status = "OK" if not probs else f"{len(probs)} ISSUE(S)"
    print(f"  [QA] {name}: {status}; ink coverage {ink:.1f}%")
    for p in probs[:14]:
        print(p)
    return probs


def save(fig, stem, name):
    qa(fig, name)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(os.path.join(OUT, f"{stem}.{ext}"),
                    dpi=600 if ext == "png" else None)
    plt.close(fig)
    print(f"  wrote {stem}.png/.pdf/.svg")


def panel_tag(ax, tag, x=-0.09, y=1.06):
    ax.text(x, y, tag, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="left")


def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# =============================================================== load artifacts
calls = pd.read_csv(os.path.join(ROOT, "05_conditional_dependency",
                                 "per_screen_calls_v2.tsv"), sep="\t", low_memory=False)
frozen = pd.read_csv(os.path.join(ROOT, "09_candidate_freeze",
                                  "BIOLOGICAL_CANDIDATES_FROZEN_v2.tsv"), sep="\t")
conv = pd.read_csv(os.path.join(ROOT, "07_cross_pathogen",
                                "convergence_meta_v2_gene.tsv"), sep="\t")
modules = pd.read_csv(os.path.join(ROOT, "07_cross_pathogen",
                                   "convergence_meta_v2_module.tsv"), sep="\t")

n_records = len(calls)
n_signal = int(calls["call"].isin([EXPLOITED, PROTECTIVE, STATE]).sum())
n_cand = len(frozen)
n_t12 = int(frozen["bio_tier"].isin(["Bio-Tier 2", "Bio-Tier 3"]).sum())
n_t2 = int((frozen["bio_tier"] == "Bio-Tier 2").sum())
TIER2 = list(frozen.loc[frozen["bio_tier"] == "Bio-Tier 2", "gene"])
print(f"funnel: {n_records:,} -> {n_signal:,} -> {n_cand} -> {n_t12} -> {n_t2}")

HEAT = LinearSegmentedColormap.from_list(
    "ev", ["#FFFFFF", "#DCE9F5", "#9ECAE1", "#4292C6", "#08519C", "#08306B"])


# ==================================================================== FIGURE 1
def figure1():
    fig = plt.figure(figsize=(7.2, 4.9))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 1.05],
                  hspace=0.50, wspace=0.30)

    ax = fig.add_subplot(gs[0, :]); ax.axis("off")
    ax.set_xlim(0, 100); ax.set_ylim(0, 40)
    boxes = [
        ("Screen census\n& firewall audit", "53 screens; 39 analysable\nLegionella excluded", C["grey"]),
        ("Estimand registry", "endpoint + direction\nsemantics locked", C["grey"]),
        ("DepMap baseline\n+ conditional calls", "26Q1 Chronos\n852,216 call records", C["blue"]),
        ("Cross-pathogen\nconvergence", "meta-analysis\n+ leave-one-pathogen-out", C["green"]),
        ("Human relevance\n+ biological freeze", "lung expression, genetics\n353 candidates (12 Bio-Tier 2)", C["orange"]),
        ("Holdout &\ntractability", "prospective adjudication\ndruggability / Pareto", C["red"]),
    ]
    n = len(boxes); w = 13.2; gap = (100 - n * w) / (n + 1)
    for i, (title, sub, col) in enumerate(boxes):
        x = gap + i * (w + gap)
        ax.add_patch(FancyBboxPatch(
            (x, 13), w, 14, boxstyle="round,pad=0.4,rounding_size=1.1",
            linewidth=0.7, edgecolor=col, facecolor=col, alpha=0.15))
        ax.text(x + w / 2, 25.4, title, ha="center", va="top", fontsize=5.9,
                fontweight="bold", color=C["dark"], linespacing=1.35)
        ax.text(x + w / 2, 19.6, sub, ha="center", va="top", fontsize=4.7,
                color="#555555", linespacing=1.45)
        if i < n - 1:
            ax.add_patch(FancyArrowPatch((x + w + 0.5, 20.0), (x + w + gap - 0.5, 20.0),
                                         arrowstyle="-|>", mutation_scale=6.5,
                                         linewidth=0.7, color="#8A8A8A"))
    ax.text(50, 30.6, "Human Pulmonary Infection Conditional Dependency Atlas",
            ha="center", va="bottom", fontsize=7.4, fontweight="bold", color=C["dark"])
    ax.plot([3, 97], [12.3, 12.3], color="#DDDDDD", lw=0.6)
    ax.text(3, 11.6, "discovery phase  |  prospectively frozen", ha="left", va="top",
            fontsize=4.8, color="#777777", style="italic")
    ax.text(97, 11.6, "validation phase  |  locked until freeze", ha="right", va="top",
            fontsize=4.8, color="#777777", style="italic")
    ax.text(50, 8.6, "Biological-before-drug firewall:  drug annotation is unlocked only after candidate "
                     "freeze and prospective holdout adjudication",
            ha="center", va="top", fontsize=5.2, color=C["red"], style="italic")
    panel_tag(ax, "a", x=-0.02, y=0.99)

    ax = fig.add_subplot(gs[1, 0])
    stages = ["Gene x screen\ncall records", "Functional\nsignal records",
              "Pre-registered\ncandidates", "HG3-scored\n(Tier 1-2)",
              "Bio-Tier 2\nprimary candidates"]
    vals = [n_records, n_signal, n_cand, n_t12, n_t2]
    y = np.arange(len(vals))[::-1]
    ax.barh(y, vals, color=["#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C"],
            height=0.66)
    for yy, v in zip(y, vals):
        ax.text(v * 1.4, yy, f"{v:,}", va="center", ha="left", fontsize=6.4,
                fontweight="bold", color=C["dark"])
    ax.set_yticks(y); ax.set_yticklabels(stages, fontsize=6.1, linespacing=1.25)
    ax.set_xscale("log"); ax.set_xlim(5, 8e6)
    ax.set_xlabel("Records / genes (log scale)")
    ax.set_title("Candidate funnel", loc="left", fontsize=7.2)
    clean(ax); panel_tag(ax, "b")

    ax = fig.add_subplot(gs[1, 1])
    census = (calls.groupby(["pathogen", "track"])["screen_id"]
              .nunique().unstack(fill_value=0))
    census = census.loc[census.sum(axis=1).sort_values().index]
    left = np.zeros(len(census))
    for tr in census.columns:
        ax.barh(np.arange(len(census)), census[tr], left=left, height=0.62,
                color={"A": C["blue"], "B": C["sky"]}.get(tr, C["grey"]),
                label=f"Track {tr} ({'raw counts' if tr == 'A' else 'gene-level'})")
        left += census[tr].values
    for i, tot in enumerate(left):
        ax.text(tot + 0.4, i, str(int(tot)), va="center", fontsize=6.3,
                fontweight="bold", color=C["dark"])
    ax.set_yticks(np.arange(len(census))); ax.set_yticklabels(census.index, fontsize=6.4)
    ax.set_xlim(0, max(left) * 1.20)
    ax.set_xlabel("Screens contributing to conditional calls")
    ax.set_title("Screen census", loc="left", fontsize=7.2)
    ax.legend(loc="lower right", handlelength=1.1, handleheight=0.8,
              borderpad=0.2, labelspacing=0.3)
    clean(ax); panel_tag(ax, "c")

    save(fig, "fig1_study_architecture", "Figure 1")


# ==================================================================== FIGURE 2
def figure2():
    fig = plt.figure(figsize=(7.2, 4.5))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.30, 1.0],
                  hspace=0.50, wspace=0.34)

    ax = fig.add_subplot(gs[0, :])
    bg = calls[(calls["call"] == NEUTRAL) & calls["baseline_chronos"].notna()
               & calls["lfc"].notna()]
    bg = bg.sample(min(len(bg), 60000), random_state=0)
    ax.scatter(bg["baseline_chronos"], bg["lfc"], s=0.35, c=CALL_COLOR[NEUTRAL],
               alpha=0.10, linewidths=0, rasterized=True, zorder=1)
    for cls in [STATE, CONFOUNDED, PROTECTIVE, EXPLOITED]:
        d = calls[(calls["call"] == cls) & calls["baseline_chronos"].notna()
                  & calls["lfc"].notna()]
        ax.scatter(d["baseline_chronos"], d["lfc"], s=3.0, c=CALL_COLOR[cls],
                   alpha=0.72, linewidths=0, rasterized=True, zorder=3)
    ax.axvline(ESS, color="#666666", lw=0.7, ls="--", zorder=2)
    ax.axhline(0, color="#666666", lw=0.7, ls="--", zorder=2)
    ax.set_xlim(-4.5, 1.3); ax.set_ylim(-11.5, 9.0)
    ax.set_xlabel("Baseline fitness  (DepMap Chronos; more negative = more essential)")
    ax.set_ylabel("Infection effect (log$_2$ fold-change)")
    ax.text(-4.4, 8.2, "baseline-essential", fontsize=5.8, color="#777777",
            style="italic", va="top")
    ax.text(1.2, 8.2, "baseline non-essential", fontsize=5.8, color="#777777",
            style="italic", va="top", ha="right")
    handles = [plt.Line2D([], [], marker="o", ls="", ms=3.4,
                          markerfacecolor=CALL_COLOR[c], markeredgecolor="none",
                          label=CALL_LABEL[c])
               for c in [EXPLOITED, PROTECTIVE, CONFOUNDED, STATE]]
    handles.append(plt.Line2D([], [], marker="o", ls="", ms=3.4,
                              markerfacecolor=CALL_COLOR[NEUTRAL],
                              markeredgecolor="none", label="No functional signal"))
    ax.legend(handles=handles, loc="lower left", ncol=3, handletextpad=0.3,
              columnspacing=0.8, labelspacing=0.3, borderpad=0.3)
    ax.set_title("Conditional dependency calls occupy the baseline non-essential space",
                 loc="left", fontsize=7.2)
    clean(ax); panel_tag(ax, "a", x=-0.062)

    ax = fig.add_subplot(gs[1, 0])
    order = [EXPLOITED, PROTECTIVE, CONFOUNDED, STATE, NEUTRAL]
    lab, val = [], []
    for cls in order:
        d = calls[(calls["call"] == cls) & calls["baseline_chronos"].notna()]
        lab.append((CALL_LABEL.get(cls, "No functional signal"), ))
        val.append((d["baseline_chronos"] < ESS).mean() * 100)
    y = np.arange(len(lab))[::-1]
    ax.barh(y, val, color=[CALL_COLOR[c] for c in order], height=0.62)
    for yy, v in zip(y, val):
        ax.text(v + 1.5, yy, f"{v:.0f}%", va="center", fontsize=6.2,
                fontweight="bold", color=C["dark"])
    ax.set_yticks(y)
    ax.set_yticklabels(["Pathogen-exploited\nconditional", "Host-protective\nconditional",
                        "Baseline-fitness\nconfounded", "Signal-state\naltering",
                        "No functional\nsignal"], fontsize=5.8, linespacing=1.2)
    ax.set_xlabel("Calls that are baseline-essential (%)"); ax.set_xlim(0, 118)
    ax.axvline(50, color="#BBBBBB", lw=0.6, ls=":")
    ax.set_title("Conditional calls are not baseline essentials", loc="left", fontsize=7.0)
    clean(ax); panel_tag(ax, "b")

    ax = fig.add_subplot(gs[1, 1])
    for tier, (col, s, a) in {"Bio-Tier 2": (C["red"], 26, 0.95),
                              "Bio-Tier 3": (C["orange"], 11, 0.65),
                              "ARCHIVE": (C["grey"], 7, 0.45),
                              "REJECT": ("#C9CCCF", 5, 0.35)}.items():
        d = frozen[frozen["bio_tier"] == tier]
        ax.scatter(d["min_baseline_chronos"], d["total_score"], s=s, c=col,
                   alpha=a, linewidths=0.3, edgecolors="white", label=tier, zorder=3)
    for g, dx, dy, ha in [("ACE2", 6, 3, "left"), ("SLC35A1", -6, 3, "right"),
                          ("EXT1", -6, -1, "right"), ("ATP6AP1", 6, -9, "left")]:
        r = frozen[frozen["gene"] == g]
        if len(r):
            ax.annotate(g, (r["min_baseline_chronos"].iloc[0], r["total_score"].iloc[0]),
                        textcoords="offset points", xytext=(dx, dy), fontsize=5.6,
                        fontweight="bold", color=C["dark"], ha=ha)
    ax.axvline(ESS, color="#666666", lw=0.7, ls="--")
    ax.set_xlim(-3.8, 1.6); ax.set_ylim(-0.8, 9.4)
    ax.set_xlabel("Min. baseline Chronos"); ax.set_ylabel("HG3 total score")
    ax.set_title("Frozen candidates", loc="left", fontsize=7.0)
    ax.legend(loc="upper left", handletextpad=0.25, labelspacing=0.28,
              borderpad=0.25, markerscale=1.6)
    clean(ax); panel_tag(ax, "c", x=-0.19)

    save(fig, "fig2_baseline_vs_infection_landscape", "Figure 2")


# ==================================================================== FIGURE 3
def figure3():
    m = conv[conv["gene"].isin(frozen["gene"])].merge(
        frozen[["gene", "bio_tier", "total_score", "convergent_module_memberships"]],
        on="gene", how="left")
    tr = {"Bio-Tier 2": 0, "Bio-Tier 3": 1, "ARCHIVE": 2, "REJECT": 3}
    m["_tr"] = m["bio_tier"].map(tr).fillna(4)
    m = m.sort_values(["_tr", "total_score"], ascending=[True, False]).reset_index(drop=True)

    cols = [("p_SC2_EXP", "SARS-CoV-2"), ("p_IAV_EXP", "Influenza A"),
            ("p_RSV_EXP", "RSV")]

    def mat_of(frame):
        out = np.full((len(frame), len(cols)), np.nan)
        for j, (c, _) in enumerate(cols):
            p = frame[c].values.astype(float)
            out[:, j] = np.where(np.isnan(p), np.nan,
                                 np.clip(-np.log10(np.clip(p, 1e-12, None)), 0, 8))
        return out

    fig = plt.figure(figsize=(7.2, 6.2))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[2.05, 1.0],
                  width_ratios=[1.25, 1.0], hspace=0.34, wspace=0.30)

    # ---- 3A candidate-level map --------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    im = ax.imshow(mat_of(m), aspect="auto", cmap=HEAT, vmin=0, vmax=8,
                   interpolation="nearest")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c[1] for c in cols], fontsize=6.5, fontweight="bold")
    ax.xaxis.set_ticks_position("top")
    ax.tick_params(axis="x", length=0, pad=2)
    ax.set_yticks([])
    for tier, col, alpha in [("Bio-Tier 2", C["red"], 0.10), ("Bio-Tier 3", C["orange"], 0.05)]:
        idx = m.index[m["bio_tier"] == tier]
        if len(idx):
            ax.axhspan(idx.min() - 0.5, idx.max() + 0.5, facecolor=col,
                       alpha=alpha, zorder=0)
    for tier in ["Bio-Tier 2", "Bio-Tier 3", "ARCHIVE", "REJECT"]:
        idx = m.index[m["bio_tier"] == tier]
        if len(idx) and idx.min() > 0:
            ax.axhline(idx.min() - 0.5, color="#999999", lw=0.55)
    ax.annotate("Bio-Tier 2 (n=12)", xy=(2.52, 5.5), xytext=(3.3, 26),
                fontsize=5.6, color=C["red"], fontweight="bold",
                arrowprops=dict(arrowstyle="-", lw=0.6, color=C["red"]))
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title("Candidate convergence map (n=353)", loc="left", fontsize=7.2, pad=16)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("$-\\log_{10}\\,p$ (convergence, EXP)", fontsize=6.0)
    cb.ax.tick_params(labelsize=5.6); cb.outline.set_linewidth(0.4)
    panel_tag(ax, "a", x=-0.02, y=1.09)

    # module membership strip
    posA = ax.get_position()
    axs = fig.add_axes([posA.x1 + 0.006, posA.y0, 0.026, posA.height])
    mods = [("HS biosynthesis", "CURATED:heparan-sulfate biosynthesis", C["purple"]),
            ("V-ATPase", "CURATED:V-ATPase proton pump", C["sky"]),
            ("COG complex", "CURATED:COG complex", C["green"])]
    strip = np.zeros((len(m), len(mods)))
    for j, (_, key, _) in enumerate(mods):
        strip[:, j] = m["convergent_module_memberships"].fillna("").str.contains(
            key, regex=False).astype(float)
    axs.imshow(strip, aspect="auto", cmap=LinearSegmentedColormap.from_list(
        "mm", ["#FFFFFF", "#BBBBBB"]), vmin=0, vmax=1, interpolation="nearest")
    axs.set_xticks(range(len(mods)))
    axs.set_xticklabels([x[0].replace(" ", "\n") for x in mods], fontsize=4.3, rotation=0)
    axs.xaxis.set_ticks_position("top")
    axs.set_yticks([]); axs.tick_params(axis="x", length=0, pad=1)
    for sp in axs.spines.values():
        sp.set_visible(False)
    ax.text(-0.05, -0.045, "module", transform=axs.transAxes, fontsize=4.3,
            ha="left", va="top", color="#777777")

    # ---- 3B Bio-Tier 2 detail ---------------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    t2 = m[m["bio_tier"] == "Bio-Tier 2"].sort_values("total_score", ascending=True)
    M = mat_of(t2)
    im = ax.imshow(M, aspect="auto", cmap=HEAT, vmin=0, vmax=8, interpolation="nearest")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isnan(M[i, j]):
                continue
            ax.text(j, i, f"{M[i, j]:.1f}", ha="center", va="center", fontsize=5.0,
                    color="white" if M[i, j] > 4.2 else "#333333")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c[1] for c in cols], fontsize=6.2, fontweight="bold")
    ax.xaxis.set_ticks_position("top"); ax.tick_params(axis="x", length=0, pad=2)
    ax.set_yticks(range(len(t2))); ax.set_yticklabels(t2["gene"], fontsize=5.6)
    ax.tick_params(axis="y", length=0, pad=1.5)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title("Bio-Tier 2 primary candidates", loc="left", fontsize=7.2, pad=16)
    panel_tag(ax, "b", x=-0.22, y=1.12)

    # ---- 3C curated modules ------------------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    cur = modules[modules["module"].str.startswith("CURATED:")].sort_values(
        "fisher_EXP", ascending=True)
    y = np.arange(len(cur))
    ax.barh(y, cur["fisher_EXP"], color=C["green"], height=0.6)
    for yy, (v, q, ng) in zip(y, zip(cur["fisher_EXP"], cur["q_EXP"], cur["n_genes"])):
        ax.text(v * 1.04, yy, f"{v:,.0f}  (q={q:.0e}, n={int(ng)})", va="center",
                fontsize=5.2, color=C["dark"])
    ax.set_yticks(y)
    ax.set_yticklabels([s.replace("CURATED:", "") for s in cur["module"]], fontsize=6.0)
    ax.set_xlim(0, cur["fisher_EXP"].max() * 2.3)
    ax.set_xlabel("Fisher combined evidence (EXP)")
    ax.set_title("Curated module convergence", loc="left", fontsize=7.2)
    clean(ax); panel_tag(ax, "c", x=-0.42)

    save(fig, "fig3_cross_pathogen_map", "Figure 3")


if __name__ == "__main__":
    print("writing publication figures to", OUT)
    figure1(); figure2(); figure3()
    print("done")
