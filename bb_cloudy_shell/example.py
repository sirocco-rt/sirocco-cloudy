"""
Thin shell blackbody test: grids over blackbody temperature with different photon banding
strategies, each compared to the Cloudy results in cld_data_noce.

    python example.py           run the Sirocco grids, collate them into sirocco/, and plot
    python example.py --plot    remake the plots from the stored summaries in sirocco/

For each banding configuration it
1. runs the grid of sirocco models with run_shell.run_grid, in sirocco/<grid name>/,
2. collates the outputs into sirocco/<grid name>/summary.ecsv and checks U, xi, nH and T_e against the inputs,
3. plots the ion fractions against Cloudy in figures/.
Running needs the Sirocco atomic data, from a data link in this folder or $SIROCCO/xdata.
"""
import glob
import os
import sys

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "thin_shell"))
import collate
import run_shell

# define the banding for each model
nbands = [0, 0, 0, 0, 1, 5, 10, 20]
banding = ["cv", "tde_bb", "AGN", "T_star", "logarithmic",
           "logarithmic", "logarithmic", "logarithmic"]

# number of photons
nphotons = 100000

# location of cloudy data, and where the sirocco grids go
cloudy_dir = "cld_data_noce"
sirocco_dir = "sirocco"

# version of sirocco
version = "1.4"

# number of points to sample
npts = 61

# if you want to add flags for sirocco
flags = ""

# number of models to run at once
nproc = 4

if "--plot" in sys.argv:
    for folder in sorted(glob.glob(os.path.join(sirocco_dir, "*", "summary.ecsv"))):
        folder = os.path.dirname(folder)
        collate.plot(collate.load(folder), folder, cloudy_dir=cloudy_dir, suffix="bb")
    sys.exit()

for i in range(8):
    # run the sirocco models: 10^3.5 to 10^6 K blackbodies of radius 1e10 cm, shell at 1e11 cm
    folder = run_shell.run_grid(loop="t_r", lmin=3.5, lmax=6.0, npts=npts, te=1e4, nh=1e7, version=version,
                                flags=flags, macro_atoms=False, banding=banding[i], nbands=nbands[i],
                                nphotons=nphotons, nproc=nproc, outdir=sirocco_dir)

    # summarize the data and check the shell conditions
    tab = collate.collate(folder)
    collate.check(tab)

    # make plots compared to cloudy
    collate.plot(tab, folder, cloudy_dir=cloudy_dir, suffix="bb")
