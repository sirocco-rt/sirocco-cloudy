"""
Run a grid of thin shell Sirocco models, looping over one parameter.

The loop parameter is one of
    t_r   blackbody temperature, at fixed central object radius (--rstar)
    U     ionization parameter, setting the source normalisation for the chosen SED
    xi    ionization parameter, as above
The physical parameters are the electron temperature (fixed with sirocco -f), the hydrogen
density and the SED; the others control the Sirocco run.

Each grid is written to its own folder (inside --outdir), which holds the parameter files, the
windsave2table outputs and inputs.json, a record of the inputs for collate.py. The scripts can be
run from any folder: the atomic data is taken from ./data if it exists, otherwise $SIROCCO/xdata.

Example:
    python run_shell.py --version 1.4 --loop U --lmin -4 --lmax 3 --npts 15 --sed power --alpha -0.9
"""
import argparse
import json
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from tqdm import tqdm

import sed as seds

PF_TEMPLATE = """System_type(star,cv,bh,agn,previous)                   bh
Central_object.mass(msol)                  0.8
Central_object.radius(cm)                  {rstar:.6e}
Binary.mass_sec(msol)                      0.4
Binary.period(hr)                          1000
Disk.type(none,flat,vertically.extended,rmin>central.obj.rad)                 none
Central_object.radiation(yes,no)                  yes
Central_object.rad_type_to_make_wind(bb,models,power,cloudy,brems,mono)       {rad_type}
{sed_block}
Central_object.geometry_for_source(sphere,lamp_post,bubble)               sphere
Wind.number_of_components                  1
Wind.type(SV,star,hydro,corona,kwd,homologous,shell,imported)                shell
Photons_per_cycle                          {nphotons}
Ionization_cycles                          {cycles}
Spectrum_cycles                            0
Wind.ionization(on.the.spot,ML93,LTE_tr,LTE_te,fixed,matrix_bb,matrix_pow,matrix_est)           matrix_pow
Line_transfer(pure_abs,pure_scat,sing_scat,escape_prob,thermal_trapping,macro_atoms_escape_prob,macro_atoms_thermal_trapping)     {rt_mode}
Matom_transition_mode(mc_jumps,matrix) matrix
Wind.radiation(yes,no)                           no
Atomic_data                                {masterfile}
Surface.reflection.or.absorption(reflect,absorb,thermalized.rerad)              reflect
Wind_heating.extra_processes(none,adiabatic,nonthermal,both)                 none
Shell.wind_mdot(msol/yr)                   {mdot:.6e}
Shell.wind.radmin(cm)                      {rshell:.6e}
Shell.wind.radmax(cm)                      {rshell_max:.12e}
Shell.wind_v_at_rmin(cm)                   1
Shell.wind.v_at_rmax(cm)                   1.000010
Shell.wind.acceleration_exponent           1
Wind.t.init                                {te:.6e}
Wind.filling_factor(1=smooth,<1=clumped)   1
Reverb.type(none,photon,wind,matom)                 none
Photon_sampling.approach(T_star,cv,yso,AGN,tde_bb,min_max_freq,user_bands,cloudy_test,wide,logarithmic)          {banding}
Photon_sampling.nbands      {nbands}
Photon_sampling.low_energy_limit(eV)    {emin}
Photon_sampling.high_energy_limit(eV)  {emax}
"""

SED_BLOCKS = {
    "bb": "Central_object.blackbody_temp(K)           {t:.6e}",
    "power": "Central_object.luminosity(ergs/s)          {lum:.6e}\nCentral_object.power_law_index             {alpha}",
    "cloudy": "Central_object.luminosity(ergs/s)          {lum:.6e}\nCentral_object.power_law_index             {alpha}\n"
              "Central_object.cloudy.low_energy_break(ev)  {low}\nCentral_object.cloudy.high_energy_break(ev) {hi}",
    "brem": "Central_object.luminosity(ergs/s)          {lum:.6e}\nCentral_object.bremsstrahlung_temp(K)       {t:.6e}\n"
            "Central_object.bremsstrahlung_alpha        {alpha}",
    "model": "Input_spectra.model_file                   {ls}\nCentral_object.luminosity(ergs/s)          {lum:.6e}",
}
DEFAULT_MASTERFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "standard80_noce.dat")
RAD_TYPES = {"bb": "bb", "power": "power", "cloudy": "cloudy", "brem": "brems", "model": "models"}
MODEL_LS = "model_sed.ls"  # the models list file written in each grid folder for a model SED


def make_pf(sed, norm, te, mdot, rshell, masterfile, macro_atoms=False, banding="logarithmic", nbands=10,
            nphotons=100000, cycles=2, emin=0.1, emax=1e5):
    """the parameter file text for one thin shell model. norm is R_star for a blackbody, L(2-10 keV) otherwise"""
    rstar = norm if sed["type"] == "bb" else 1e10
    if rstar >= rshell:
        raise ValueError("central object radius {:.3e} is outside the shell at {:.3e}".format(rstar, rshell))
    return PF_TEMPLATE.format(
        rstar=rstar, rad_type=RAD_TYPES[sed["type"]], sed_block=SED_BLOCKS[sed["type"]].format(lum=norm, ls=MODEL_LS, **sed),
        nphotons=nphotons, cycles=cycles, masterfile=masterfile, mdot=mdot, rshell=rshell,
        rshell_max=rshell * (1 + 1e-11), te=te, banding=banding, nbands=nbands, emin=emin, emax=emax,
        rt_mode="macro_atoms_thermal_trapping" if macro_atoms else "thermal_trapping")


def run_model(root, folder, version, flags=""):
    """run sirocco and windsave2table on root.pf in folder, tidy up, and return True on success"""
    def run(cmd):
        with open(os.path.join(folder, root + ".out.txt"), "a") as out:
            subprocess.run(cmd, cwd=folder, stdout=out, stderr=subprocess.STDOUT, check=True)

    run(["sirocco-{}".format(version), "-f"] + flags.split() + [root])
    with open(os.path.join(folder, root + ".sig")) as f:
        if "COMPLETE" not in f.readlines()[-1]:
            return False
    run(["windsave2table-{}".format(version), root])
    for suffix in (".wind_save", ".sig", ".out.pf", ".out.pf.old", ".spec_tot", ".log_spec_tot",
                   ".spec_tot_wind", ".log_spec_tot_wind", ".spec", ".log_spec"):
        if os.path.exists(os.path.join(folder, root + suffix)):
            os.remove(os.path.join(folder, root + suffix))
    shutil.rmtree(os.path.join(folder, "diag_" + root), ignore_errors=True)
    return True


def find_data():
    """the Sirocco atomic data folder: ./data if it exists, otherwise $SIROCCO/xdata"""
    for d in ("data", os.path.join(os.environ.get("SIROCCO", ""), "xdata")):
        if os.path.isdir(d):
            return os.path.realpath(d)
    raise FileNotFoundError("no atomic data: make a data link here or set $SIROCCO")


def grid_name(loop, sed, version, macro_atoms, masterfile, banding, nbands, nphotons, te, nh):
    def label(k, v):
        return os.path.splitext(os.path.basename(v))[0] if k == "file" else "{}{:g}".format(k, v)
    sed_label = "bb" if loop == "t_r" else "_".join([sed["type"]] + [label(k, v) for k, v in sed.items() if k != "type"])
    band = "log{:02d}".format(nbands) if banding == "logarithmic" else banding
    return "{}_{}_v{}_{}_{}_{}_np{:.0e}_te{:.0e}_nh{:.0e}".format(
        loop, sed_label, version, "matom" if macro_atoms else "basic",
        os.path.basename(masterfile).replace(".dat", ""), band, nphotons, te, nh)


def run_grid(loop="t_r", lmin=3.5, lmax=6.0, npts=61, sed=None, te=1e4, nh=1e7, rshell=1e11, rstar=1e10,
             version="1.4", flags="", macro_atoms=False, masterfile=DEFAULT_MASTERFILE, banding="logarithmic",
             nbands=10, nphotons=100000, cycles=2, emin=0.1, emax=1e5, nproc=1, outdir=".", folder=None):
    """run a grid of models with log10(loop parameter) from lmin to lmax and return the folder name.

    The folder is outdir/<name from the parameters> unless given. Sirocco is run inside it, with
    a data link to the atomic data.
    """
    sed = sed or {"type": "bb", "t": 1e5}
    if loop == "t_r" and sed["type"] != "bb":
        raise ValueError("looping over t_r needs a blackbody SED")
    if sed["type"] == "cloudy" and banding != "cloudy_test":
        print("the cloudy SED needs cloudy_test photon sampling: using it instead of {}".format(banding))
        banding = "cloudy_test"
    if sed["type"] == "model":
        # Sirocco can only make photons within the wavelength range of the model
        elo, ehi = seds.model_energy_range(sed["file"])
        if emin < elo or emax > ehi:
            emin, emax = max(emin, 1.001 * elo), min(emax, ehi / 1.001)
            print("limiting photon energies to the model SED range: {:.3g} to {:.3g} eV".format(emin, emax))
    folder = folder or os.path.join(outdir, grid_name(loop, sed, version, macro_atoms, masterfile, banding, nbands,
                                                      nphotons, te, nh))
    os.makedirs(folder, exist_ok=True)
    data_dir = find_data()
    if not os.path.exists(os.path.join(folder, "data")):
        os.symlink(data_dir, os.path.join(folder, "data"))
    shutil.copy(masterfile, folder)
    if sed["type"] == "model":
        shutil.copy(sed["file"], folder)
        with open(os.path.join(folder, MODEL_LS), "w") as f:
            f.write("{} -1\n".format(os.path.basename(sed["file"])))
    mdot = seds.shell_mdot(nh, rshell, 1.0, seds.rho2nh(masterfile, data_dir))

    models = []
    for logval in np.linspace(lmin, lmax, npts):
        model_sed = {**sed, "t": 10**logval} if loop == "t_r" else sed
        norm = rstar if loop == "t_r" else seds.norm_from_ip(sed, 10**logval, loop, rshell, nh, emax)
        u, xi = seds.ionization_parameters(model_sed, norm, rshell, nh, emax)
        root = "shell_{}{:+.3f}".format(loop, logval)
        with open(os.path.join(folder, root + ".pf"), "w") as f:
            f.write(make_pf(model_sed, norm, te, mdot, rshell, os.path.basename(masterfile), macro_atoms,
                            banding, nbands, nphotons, cycles, emin, emax))
        models.append(dict(root=root, logval=logval, sed=model_sed, norm=norm, U=u, xi=xi, t_e=te, nh=nh))

    with open(os.path.join(folder, "inputs.json"), "w") as f:
        json.dump(dict(loop=loop, version=version, rshell=rshell, mdot=mdot, models=models), f, indent=1)

    print("running {}".format(folder))
    with ThreadPoolExecutor(nproc) as pool:
        ok = list(tqdm(pool.map(lambda m: run_model(m["root"], folder, version, flags), models), total=npts))
    for m, success in zip(models, ok):
        if not success:
            print("Failed on {}".format(m["root"]))
    return folder


def parse_arguments():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", required=True, help="suffix of the sirocco-<version> executable")
    p.add_argument("--loop", choices=["t_r", "U", "xi"], default="t_r", help="parameter to loop over")
    p.add_argument("--lmin", type=float, default=3.5, help="log10 of the first loop value")
    p.add_argument("--lmax", type=float, default=6.0, help="log10 of the last loop value")
    p.add_argument("--npts", type=int, default=61, help="number of models")
    p.add_argument("--sed", choices=seds.SED_TYPES, default="bb", help="SED type")
    p.add_argument("--t_sed", type=float, default=1e5, help="blackbody or bremsstrahlung temperature (K)")
    p.add_argument("--alpha", type=float, default=-0.9, help="power law or bremsstrahlung index")
    p.add_argument("--low", type=float, default=1.0, help="cloudy SED low energy break (eV)")
    p.add_argument("--hi", type=float, default=1e4, help="cloudy SED high energy break (eV)")
    p.add_argument("--sed_file", help="model SED file (wavelength in A, F_lambda), e.g. from data/seds/model_files")
    p.add_argument("--te", type=float, default=1e4, help="electron temperature (K)")
    p.add_argument("--nh", type=float, default=1e7, help="hydrogen density (cm^-3)")
    p.add_argument("--rshell", type=float, default=1e11, help="shell radius (cm)")
    p.add_argument("--rstar", type=float, default=1e10, help="blackbody radius for --loop t_r (cm)")
    p.add_argument("--flags", default="", help="extra flags for sirocco")
    p.add_argument("--macro_atoms", action="store_true", help="use macro atom line transfer")
    p.add_argument("--masterfile", default=DEFAULT_MASTERFILE, help="atomic data masterfile")
    p.add_argument("--banding", default="logarithmic", help="photon sampling approach")
    p.add_argument("--nbands", type=int, default=10, help="number of bands for logarithmic banding")
    p.add_argument("--nphotons", type=int, default=100000, help="photons per cycle")
    p.add_argument("--cycles", type=int, default=2, help="ionization cycles")
    p.add_argument("--nproc", type=int, default=1, help="number of models to run at once")
    p.add_argument("--outdir", default=".", help="folder in which to make the grid folder")
    return p.parse_args()


def sed_from_args(args):
    return {"bb": {"type": "bb", "t": args.t_sed},
            "power": {"type": "power", "alpha": args.alpha},
            "cloudy": {"type": "cloudy", "alpha": args.alpha, "low": args.low, "hi": args.hi},
            "brem": {"type": "brem", "t": args.t_sed, "alpha": args.alpha},
            "model": {"type": "model", "file": args.sed_file}}[args.sed]


if __name__ == "__main__":
    a = parse_arguments()
    run_grid(a.loop, a.lmin, a.lmax, a.npts, sed_from_args(a), a.te, a.nh, a.rshell, a.rstar, a.version, a.flags,
             a.macro_atoms, a.masterfile, a.banding, a.nbands, a.nphotons, a.cycles, nproc=a.nproc, outdir=a.outdir)
