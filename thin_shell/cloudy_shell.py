"""
Run Cloudy for the same thin shell models as run_shell.run_grid.

run_cloudy_grid takes the same keywords as run_shell.run_grid, so one grid definition can be
passed to both; the Sirocco-only ones (version, banding, ...) are ignored. Each model is set by its
ionization parameter (xi for a xi loop, U otherwise), using Cloudy's own continua where they match
the Sirocco SED: blackbody, table power law (Sirocco's cloudy SED) and bremsstrahlung (alpha = -0.2
only). Other SEDs (model, or brem with another alpha) are written once per grid as a table SED over
Sirocco's photon energy range. A plain power law is refused: Cloudy has no matching continuum, so
use the cloudy SED, which both codes share. The abundances come from the Sirocco masterfile's
element file, with the other elements turned off, and charge exchange is turned off for a *_noce masterfile unless charge_exchange is given.

The results are collated into <folder>/cloudy_summary.ecsv, with the same columns as Sirocco's
summary.ecsv, for collate.plot. Cloudy is run as $CLOUDY_EXE (default cloudy.exe) with -r <root>.

Example:
    python cloudy_shell.py --loop U --lmin -4 --lmax 3 --npts 29 --sed cloudy --alpha -0.9 --nproc 4
"""
import inspect
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from astropy.table import Table
from tqdm import tqdm

import run_shell
import sed as seds

cld_data_path = os.environ["CLOUDY_DATA_PATH"]
cld_source_folder = "{}/source".format(os.path.dirname(cld_data_path.split(os.pathsep)[0]))

# Cloudy's elements, H to Zn
CLOUDY_ELEMENTS = ["hydrogen", "helium", "lithium", "beryllium", "boron", "carbon", "nitrogen", "oxygen", "fluorine",
                   "neon", "sodium", "magnesium", "aluminium", "silicon", "phosphorus", "sulphur", "chlorine", "argon",
                   "potassium", "calcium", "scandium", "titanium", "vanadium", "chromium", "manganese", "iron",
                   "cobalt", "nickel", "copper", "zinc"]
SYMBOLS = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K",
           "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"]


RYD = 13.605693  # eV
SED_FILE = "sed_shape.sed"  # the table SED written in the grid folder, when needed


def sed_table(sed, emin, emax, npts=2000):
    """a Cloudy SED file: F_nu against photon energy in eV, from emin to emax"""
    energy = np.logspace(np.log10(emin), np.log10(emax), npts)
    fnu = seds.lnu(sed, energy / seds.HEV)
    fnu = np.maximum(fnu, 1e-40 * fnu.max())  # Cloudy needs positive values
    lines = ["{:.6e} {:.6e}".format(e, f) for e, f in zip(energy, fnu)]
    lines[0] += " Fnu units eV"
    return "# {}\n{}\n***\n".format(json.dumps(sed), "\n".join(lines))


def log(x):
    """log10(x) for an input line, without printing -0.000000"""
    return "{:.6f}".format(round(np.log10(x), 6) + 0.0)


def continuum(sed):
    """the Cloudy command for the shape of a Sirocco SED, or None if it needs a table SED"""
    if sed["type"] == "bb":
        return "blackbody {:.6e}".format(sed["t"])
    if sed["type"] == "cloudy":
        return "table power law spectral index {}, low = {:.6g}, hi = {:.6g}".format(
            sed["alpha"], sed["low"] / RYD, sed["hi"] / RYD)
    if sed["type"] == "brem" and sed["alpha"] == -0.2:
        return "bremsstrahlung {:.6e}".format(sed["t"])
    return None


def make_input(model, loop, rshell, elements, emin, emax, charge_exchange=False, extra=()):
    """the Cloudy input for one thin shell model; elements is {symbol: (Z, abundance, weight)}"""
    present = {z: abun for z, abun, weight in elements.values() if z <= len(CLOUDY_ELEMENTS)}
    lines = [continuum(model["sed"]) or 'table SED "{}"'.format(SED_FILE),
             "xi " + log(model["xi"]) if loop == "xi" else "ionization parameter " + log(model["U"]),
             "radius {} {}".format(log(rshell), log(rshell * 1e-11)),  # 1e-11 r thick, as in Sirocco
             "sphere",
             "hden " + log(model["nh"]),
             "constant temperature {} log".format(log(model["t_e"]))]
    for z, name in enumerate(CLOUDY_ELEMENTS[1:], start=2):
        lines.append("element {} abundance {:.4f}".format(name, np.log10(present[z])) if z in present
                     else "element {} off".format(name))
    if not charge_exchange:
        lines.append("no charge exchange")
    lines += list(extra)
    root = model["root"]
    lines += ['save last element {} "{}.{}.ele"'.format(CLOUDY_ELEMENTS[z - 1], root, z) for z in sorted(present)]
    lines += ['save last physical conditions "{}.pcon"'.format(root),
              'save last heating "{}.heat"'.format(root),
              'save last cooling "{}.cool"'.format(root)]
    return "\n".join(lines) + "\n"


def run_model(root, folder, exe):
    """run Cloudy on root.in in folder and return True if it finished. Cloudy's exit code is nonzero
    if there were warnings, so finishing is judged from the output; check root.out for warnings.
    Cloudy only looks for a table SED on its data path, so the model folder (.) is put first on it"""
    path = os.pathsep.join([".", os.environ.get("CLOUDY_DATA_PATH", "+")])
    with open(os.path.join(folder, root + ".log"), "w") as log:
        subprocess.run([exe, "-r", root], cwd=folder, stdout=log, stderr=subprocess.STDOUT,
                       env=dict(os.environ, CLOUDY_DATA_PATH=path))
    out = os.path.join(folder, root + ".out")
    return os.path.exists(out) and " Cloudy ends:" in open(out).read()


def run_cloudy_grid(loop="t_r", lmin=3.5, lmax=6.0, npts=61, sed=None, te=1e4, nh=1e7, rshell=1e11, rstar=1e10,
                    masterfile=run_shell.DEFAULT_MASTERFILE, emin=0.1, emax=1e5, nproc=1, outdir=".", folder=None,
                    charge_exchange=None, extra=(), cloudy_exe=None, **sirocco_only):
    """run Cloudy for the models of run_shell.run_grid with the same keywords, collate, and return the folder"""
    unknown = set(sirocco_only) - set(inspect.signature(run_shell.run_grid).parameters)
    if unknown:
        raise TypeError("unknown keywords {}".format(sorted(unknown)))
    sed = sed or {"type": "bb", "t": 1e5}
    if sed["type"] == "power":
        raise ValueError("Cloudy has no continuum matching Sirocco's power law: use the cloudy SED")
    default_ce = "noce" not in os.path.basename(masterfile)
    charge_exchange = default_ce if charge_exchange is None else charge_exchange
    emin, emax = run_shell.photon_range(sed, emin, emax)
    models = run_shell.make_models(loop, lmin, lmax, npts, sed, te, nh, rshell, rstar, emax)
    elements = seds.elements(masterfile, run_shell.find_data())

    ce_label = "" if charge_exchange == default_ce else ("_ce" if charge_exchange else "_noce")
    folder = folder or os.path.join(outdir, "cloudy_{}_{}_{}{}_te{:.0e}_nh{:.0e}".format(
        loop, run_shell.sed_label(loop, sed), os.path.basename(masterfile).replace(".dat", ""), ce_label, te, nh))
    os.makedirs(folder, exist_ok=True)
    if continuum(sed) is None:
        with open(os.path.join(folder, SED_FILE), "w") as f:
            f.write(sed_table(sed, emin, emax))
    for m in models:
        with open(os.path.join(folder, m["root"] + ".in"), "w") as f:
            f.write(make_input(m, loop, rshell, elements, emin, emax, charge_exchange, extra))
    with open(os.path.join(folder, "inputs.json"), "w") as f:
        json.dump(dict(loop=loop, code="cloudy", rshell=rshell, emin=emin, emax=emax,
                       charge_exchange=charge_exchange, models=models), f, indent=1)

    exe = "{}/cloudy.exe".format(cld_source_folder) if cloudy_exe is None else cloudy_exe
    print("running {}".format(folder))
    with ThreadPoolExecutor(nproc) as pool:
        ok = list(tqdm(pool.map(lambda m: run_model(m["root"], folder, exe), models), total=npts))
    for m, success in zip(models, ok):
        if not success:
            print("Failed on {}".format(m["root"]))
    collate_cloudy(folder)
    return folder


def read_element(fname, z):
    """the ion fractions (Z + 1 stages) in the last zone of a Cloudy 'save element' file"""
    data = np.loadtxt(fname, ndmin=2)
    return data[-1, 1:z + 2]


def collate_cloudy(folder):
    """make folder/cloudy_summary.ecsv, one row per model, with the same columns as Sirocco's summary.ecsv"""
    info = json.load(open(os.path.join(folder, "inputs.json")))
    rows = []
    for m in info["models"]:
        row = dict(logval=m["logval"], root=m["root"], norm=m["norm"], in_U=m["U"], in_xi=m["xi"],
                   in_nh=m["nh"], in_t_e=m["t_e"])
        for z, symbol in enumerate(SYMBOLS, start=1):
            fname = os.path.join(folder, "{}.{}.ele".format(m["root"], z))
            if os.path.exists(fname):
                row.update({"{}_i{:02d}".format(symbol, i + 1): f for i, f in enumerate(read_element(fname, z))})
        rows.append(row)
    tab = Table(rows=rows)
    tab.meta = {k: v for k, v in info.items() if k != "models"}
    tab.write(os.path.join(folder, "cloudy_summary.ecsv"), format="ascii.ecsv", overwrite=True)
    return tab



if __name__ == "__main__":
    a = run_shell.parse_arguments(sirocco=False)
    run_cloudy_grid(a.loop, a.lmin, a.lmax, a.npts, run_shell.sed_from_args(a), a.te, a.nh, a.rshell, a.rstar,
                    a.masterfile, nproc=a.nproc, outdir=a.outdir)
