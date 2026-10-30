"""
Collate a grid made by run_shell.py into summary tables, check the shell conditions
against the inputs, and plot the ion fractions, optionally against Cloudy. For a grid run
with T_e from thermal balance, also plot T_e and the heating and cooling.

summary.ecsv has the ion fractions and shell conditions, thermal.ecsv the temperature and the
total and component heating and cooling rates divided by n_H^2 (erg cm^3 s^-1):
    heating  phot (photoionization, with Auger), ff, comp (Compton), line, other
    cooling  line, ff, rr (radiative recombination), comp, other
where other is the total minus the rest (for Sirocco's cooling, mostly dielectronic
recombination and collisional ionization).

Usage:
    python collate.py <folder> [--cloudy <cloudy_dir> [--suffix bb]]
"""
import argparse
import json
import os

import matplotlib.pyplot as plt
import numpy as np
import roman
from astropy.io import ascii
from astropy.table import Table

ELEMENTS = ["H", "He", "C", "N", "O", "Fe"]
CLOUDY_NAMES = dict(H="hydrogen", He="helium", C="carbon", N="nitrogen", O="oxygen", Fe="iron")
# input quantity: matching column in the windsave2table master file
CHECKS = dict(U="ip", xi="xi", nh="nh", t_e="t_e")
XLABELS = {"t_r": r"$\log T_{\rm bb}$ (K)", "U": r"$\log U$", "xi": r"$\log \xi$"}
HEAT = ["phot", "ff", "comp", "line"]
COOL = ["line", "ff", "rr", "comp"]
# Sirocco heat.txt columns making up each component
SIROCCO_HEAT = dict(phot=["heat_phot", "heat_auge"], ff=["heat_ff"], comp=["heat_comp"], line=["heat_line"])
SIROCCO_COOL = dict(line=["lum_lines"], ff=["lum_ff"], rr=["cool_rr"], comp=["cool_comp"])


def collate(folder):
    """make folder/summary.ecsv, one row per model, with inputs, shell conditions and ion fractions"""
    info = json.load(open(os.path.join(folder, "inputs.json")))
    rows = []
    for m in info["models"]:
        row = dict(logval=m["logval"], root=m["root"], norm=m["norm"],
                   **{"in_" + k: m[k] for k in CHECKS})
        try:
            master = ascii.read(os.path.join(folder, m["root"] + ".master.txt"))
        except FileNotFoundError:
            print("No outputs for {}".format(m["root"]))
            continue
        cell = master[master["inwind"] == 0][0]
        row.update({k: cell[k] for k in ("converge", "t_e", "t_r", "ne", "nh", "ip", "xi")})
        for el in ELEMENTS:
            frac = ascii.read(os.path.join(folder, "{}.{}.frac.txt".format(m["root"], el)))
            fcell = frac[frac["inwind"] == 0][0]
            row.update({"{}_{}".format(el, c): fcell[c] for c in frac.colnames if c[0] == "i" and c[1:].isdigit()})
        rows.append(row)
    if not rows:
        raise RuntimeError("no completed models in {}".format(folder))
    tab = Table(rows=rows)
    tab.meta = {k: v for k, v in info.items() if k != "models"}
    tab.write(os.path.join(folder, "summary.ecsv"), format="ascii.ecsv", overwrite=True)
    collate_thermal(folder)
    return tab


def thermal_row(t_e, nh, heat_tot, cool_tot, heat, cool):
    """a row of a thermal table from rates per unit volume: the components are dicts of rates"""
    row = dict(t_e=t_e, heat_tot=heat_tot / nh**2, cool_tot=cool_tot / nh**2)
    row.update({"heat_" + k: heat[k] / nh**2 for k in HEAT})
    row.update({"cool_" + k: cool[k] / nh**2 for k in COOL})
    row["heat_other"] = row["heat_tot"] - sum(row["heat_" + k] for k in HEAT)
    row["cool_other"] = row["cool_tot"] - sum(row["cool_" + k] for k in COOL)
    return row


def collate_thermal(folder):
    """make folder/thermal.ecsv from Sirocco's heat.txt files"""
    info = json.load(open(os.path.join(folder, "inputs.json")))
    rows = []
    for m in info["models"]:
        fname = os.path.join(folder, m["root"] + ".heat.txt")
        if not os.path.exists(fname):
            continue
        heat = ascii.read(fname)
        c = heat[heat["inwind"] == 0][0]
        vol = c["vol"]
        rows.append(dict(logval=m["logval"], root=m["root"], converge=c["converge"], **thermal_row(
            c["t_e"], m["nh"], c["heat_tot"] / vol, c["cool_tot"] / vol,
            {k: sum(c[col] for col in cols) / vol for k, cols in SIROCCO_HEAT.items()},
            {k: sum(c[col] for col in cols) / vol for k, cols in SIROCCO_COOL.items()})))
    tab = Table(rows=rows)
    tab.meta = {k: v for k, v in info.items() if k != "models"}
    tab.write(os.path.join(folder, "thermal.ecsv"), format="ascii.ecsv", overwrite=True)
    return tab


def load(folder, name="summary"):
    """read a summary table made by collate(): name is summary or thermal"""
    return Table.read(os.path.join(folder, name + ".ecsv"))


def check(tab, rtol=0.05, umin=1e-4):
    """print the ratio of Sirocco's shell conditions to the inputs; return True if all within rtol.

    U and xi are only checked where the input U > umin, since at very low U the few
    ionizing photons are not sampled by the Monte Carlo
    """
    ok = True
    checks = CHECKS if tab.meta.get("fixed_te", True) else {k: v for k, v in CHECKS.items() if k != "t_e"}
    if "t_e" not in checks:
        print("T_e from thermal balance: {} of {} models not converged".format(np.sum(tab["converge"] != 0), len(tab)))
    sampled = np.array(tab["in_U"] > umin)
    if not sampled.all():
        print("U, xi not checked for {} models with U < {:g}".format(np.sum(~sampled), umin))
    for name, col in checks.items():
        use = sampled if name in ("U", "xi") else np.ones(len(tab), dtype=bool)
        ratio = np.array(tab[col] / tab["in_" + name])[use]
        iworst = np.argmax(np.abs(ratio - 1))
        ok &= abs(ratio[iworst] - 1) < rtol
        print("{:4s} sirocco/input: mean {:.4f}, worst {:.4f} at log {} = {:.3f}".format(
            name, ratio.mean(), ratio[iworst], tab.meta["loop"], tab["logval"][use][iworst]))
    return ok


def ion_fracs(tab, el):
    return np.array([tab[c] for c in tab.colnames if c.startswith(el + "_i")]).T


def mean_ion(fracs):
    """average ionization state, 0 for all neutral and 1 for fully stripped"""
    states = np.arange(1, fracs.shape[1] + 1)
    return (np.sum(states * fracs, axis=1) - 1) / (fracs.shape[1] - 1)


def read_cloudy(cloudy_dir, el, suffix):
    """Cloudy log10(loop parameter) and ion fractions for an element, from a folder made by
    cloudy_shell.py (cloudy_summary.ecsv) or one of cl_<element>_<suffix>.dat files"""
    summary = os.path.join(cloudy_dir, "cloudy_summary.ecsv")
    if os.path.exists(summary):
        tab = Table.read(summary)
        return np.array(tab["logval"]), ion_fracs(tab, el)
    fname = os.path.join(cloudy_dir, "cl_{}_{}.dat".format(CLOUDY_NAMES[el], suffix))
    if not os.path.exists(fname):
        return None, None
    data = np.loadtxt(fname)
    return data[:, 0], data[:, 1:]


def savefig(fig, prefix, folder):
    os.makedirs("figures", exist_ok=True)
    fig.savefig(os.path.join("figures", "{}_{}.png".format(prefix, os.path.basename(os.path.normpath(folder)))), dpi=150)
    plt.close(fig)


def plot(tab, folder, cloudy_dir=None, suffix="bb"):
    """2x3 panel plot of ion fractions (H, He, C, O, Fe) and mean ionization state"""
    xlabel = XLABELS[tab.meta["loop"]]
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    x = tab["logval"]
    for ax, el in zip(axes.flat, ["H", "He", "C", "O", "Fe"]):
        fracs = ion_fracs(tab, el)
        xc, cfracs = read_cloudy(cloudy_dir, el, suffix) if cloudy_dir else (None, None)
        for i in range(fracs.shape[1]):
            color = "C{}".format(i % 10)
            ax.semilogy(x, fracs[:, i], "-o", ms=3, color=color, label="{} {}".format(el, roman.toRoman(i + 1)))
            if xc is not None:
                ax.semilogy(xc, cfracs[:, i], "-", lw=3, alpha=0.5, color=color)
        if el != "Fe":
            ax.legend(fontsize=8)
        ax.set(xlabel=xlabel, ylabel="Fractional abundance", ylim=(1e-4, 1.2), xlim=(x.min(), x.max()))

    ax = axes.flat[-1]
    for i, el in enumerate(ELEMENTS):
        ax.plot(x, mean_ion(ion_fracs(tab, el)), "--", c="C{}".format(i))
        xc, cfracs = read_cloudy(cloudy_dir, el, suffix) if cloudy_dir else (None, None)
        if xc is not None:
            ax.plot(xc, mean_ion(cfracs), "-", lw=3, alpha=0.5, c="C{}".format(i))
        ax.plot([], [], c="C{}".format(i), label=el)
    ax.set(xlabel=xlabel, ylabel=r"$\chi$", xlim=(x.min(), x.max()))
    ax.legend()
    axes.flat[1].set_title(os.path.basename(os.path.normpath(folder)), fontsize=12)
    fig.tight_layout()
    savefig(fig, "summary", folder)


def plot_thermal(tab, folder, cloudy_dir=None):
    """T_e, heating and cooling against the loop parameter: Sirocco as points (crossed where its
    T_e has not converged), Cloudy (from cloudy_thermal.ecsv in cloudy_dir) as thick lines"""
    fname = os.path.join(cloudy_dir, "cloudy_thermal.ecsv") if cloudy_dir else ""
    ctab = Table.read(fname) if os.path.exists(fname) else None
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for t, style in ((tab, dict(ls="--", marker="o", ms=3)), (ctab, dict(ls="-", lw=3, alpha=0.5))):
        if t is None:
            continue
        axes[0].semilogy(t["logval"], t["t_e"], c="k", **style)
        for ax, kind, parts in ((axes[1], "heat", HEAT), (axes[2], "cool", COOL)):
            for i, part in enumerate(["tot"] + parts + ["other"]):
                y = np.abs(t["{}_{}".format(kind, part)])  # Cloudy's net cooling by an agent can be negative
                ax.semilogy(t["logval"], y, c="k" if part == "tot" else "C{}".format(i), **style)
    for ax, kind, parts in ((axes[1], "heat", HEAT), (axes[2], "cool", COOL)):
        for i, part in enumerate(["tot"] + parts + ["other"]):
            ax.plot([], [], c="k" if part == "tot" else "C{}".format(i), label=part)
        ax.set(ylabel=r"{}ing / $n_H^2$ (erg cm$^3$ s$^{{-1}}$)".format(kind.capitalize()))
        ax.legend(fontsize=8)
    bad = np.array(tab["converge"] != 0)
    axes[0].plot(tab["logval"][bad], tab["t_e"][bad], "rx", ms=8, label="not converged")
    axes[0].plot([], [], "k--o", ms=3, label="Sirocco")
    axes[0].plot([], [], "k-", lw=3, alpha=0.5, label="Cloudy")
    axes[0].set(ylabel=r"$T_e$ (K)")
    axes[0].legend()
    for ax in axes:
        ax.set(xlabel=XLABELS[tab.meta["loop"]], xlim=(tab["logval"].min(), tab["logval"].max()))
    axes[1].set_title(os.path.basename(os.path.normpath(folder)), fontsize=10)
    ymax = max(np.max(tab["heat_tot"]), np.max(tab["cool_tot"]))
    for ax in axes[1:]:
        ax.set_ylim(1e-6 * ymax, 3 * ymax)
    fig.tight_layout()
    savefig(fig, "thermal", folder)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("folder")
    p.add_argument("--cloudy", default=None, help="Cloudy folder from cloudy_shell.py, or with cl_<element>_<suffix>.dat files")
    p.add_argument("--suffix", default="bb", help="suffix of the Cloudy files")
    p.add_argument("--rtol", type=float, default=0.05, help="tolerance for the input check")
    p.add_argument("--umin", type=float, default=1e-4, help="smallest input U at which U and xi are checked")
    args = p.parse_args()
    tab = collate(args.folder)
    print("all within tolerance" if check(tab, args.rtol, args.umin) else "some values outside tolerance")
    plot(tab, args.folder, args.cloudy, args.suffix)
    if not tab.meta.get("fixed_te", True):
        plot_thermal(load(args.folder, "thermal"), args.folder, args.cloudy)
