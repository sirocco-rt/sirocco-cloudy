"""
The bb_cloudy_shell example, using the thin_shell scripts: blackbody temperature grids with
different photon banding strategies, each compared to the Cloudy results in bb_cloudy_shell.
Run from this folder with a data link to the Sirocco atomic data (e.g. via Setup_Sirocco_Dir).

For each banding configuration it
1. runs the grid of sirocco models with run_shell.run_grid,
2. collates the outputs into <folder>/summary.ecsv and checks U, xi, nH and T_e against the inputs,
3. plots the ion fractions against Cloudy in figures/.
"""
import sys 
sys.path.append("../thin_shell/")  # add the bb_cloudy_shell folder to the path
import collate
import run_shell

# define the banding for each model
nbands = [0, 0, 0, 0, 1, 5, 10, 20]
banding = ["cv", "tde_bb", "AGN", "T_star", "logarithmic",
           "logarithmic", "logarithmic", "logarithmic"]

# number of photons
nphotons = 100000

# location of cloudy data
cloudy_dir = "cld_data_noce"

# version of sirocco
version = "1.4"

# number of points to sample
npts = 61

# if you want to add flags for sirocco
flags = ""

# number of models to run at once
nproc = 4

for i in range(8):
    # run the sirocco models: 10^3.5 to 10^6 K blackbodies of radius 1e10 cm, shell at 1e11 cm
    sirocco_dir = run_shell.run_grid(loop="t_r", lmin=3.5, lmax=6.0, npts=npts, te=1e4, nh=1e7, version=version,
                                     flags=flags, macro_atoms=False, masterfile="standard80_noce.dat",
                                     banding=banding[i], nbands=nbands[i], nphotons=nphotons, nproc=nproc)

    # summarize the data and check the shell conditions
    tab = collate.collate(sirocco_dir)
    collate.check(tab)

    # make plots compared to cloudy
    collate.plot(tab, sirocco_dir, cloudy_dir=cloudy_dir, suffix="bb")
