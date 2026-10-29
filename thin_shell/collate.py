"""
Collate a grid made by run_shell.py into one summary table, check the shell conditions
against the inputs, and plot the ion fractions, optionally against Cloudy.

Usage:
    python collate.py <folder> [--cloudy <cloudy_dir> --suffix bb]
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
    return tab


def check(tab, rtol=0.05, umin=1e-4):
    """print the ratio of Sirocco's shell conditions to the inputs; return True if all within rtol.

    U and xi are only checked where the input U > umin, since at very low U the few
    ionizing photons are not sampled by the Monte Carlo
    """
    ok = True
    sampled = np.array(tab["in_U"] > umin)
    if not sampled.all():
        print("U, xi not checked for {} models with U < {:g}".format(np.sum(~sampled), umin))
    for name, col in CHECKS.items():
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
    """Cloudy ion fractions, first column log10 of the loop parameter, then one column per ion"""
    fname = os.path.join(cloudy_dir, "cl_{}_{}.dat".format(CLOUDY_NAMES[el], suffix))
    if not os.path.exists(fname):
        return None, None
    data = np.loadtxt(fname)
    return data[:, 0], data[:, 1:]


def plot(tab, folder, cloudy_dir=None, suffix="bb"):
    """2x3 panel plot of ion fractions (H, He, C, O, Fe) and mean ionization state"""
    xlabel = {"t_r": r"$\log T_{\rm bb}$ (K)", "U": r"$\log U$", "xi": r"$\log \xi$"}[tab.meta["loop"]]
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
    os.makedirs("figures", exist_ok=True)
    fig.savefig(os.path.join("figures", "summary_{}.png".format(os.path.basename(os.path.normpath(folder)))), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("folder")
    p.add_argument("--cloudy", default=None, help="folder with Cloudy cl_<element>_<suffix>.dat files")
    p.add_argument("--suffix", default="bb", help="suffix of the Cloudy files")
    p.add_argument("--rtol", type=float, default=0.05, help="tolerance for the input check")
    p.add_argument("--umin", type=float, default=1e-4, help="smallest input U at which U and xi are checked")
    args = p.parse_args()
    tab = collate(args.folder)
    print("all within tolerance" if check(tab, args.rtol, args.umin) else "some values outside tolerance")
    plot(tab, args.folder, args.cloudy, args.suffix)
