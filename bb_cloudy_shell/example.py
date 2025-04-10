import bb_run_cloudy_banding as run
import summarize_bb
import PlotsTemp as plot
"""
This script runs a series of simulations using the `sirocco` model with different banding configurations.
It performs the following steps for each configuration:
1. Defines the banding and number of bands for each model.
2. Runs the `sirocco` model using the `run_multi_sirocco` function with specified parameters.
3. Summarizes the simulation data using the `summarize` function from the `summarize_bb` module.
4. Generates summary plots comparing the `sirocco` results to `cloudy` data using the `summary_plot` function from the `PlotsTemp` module.
"""


# define the banding for each model
nbands = [0, 0, 0, 0, 1, 5, 10, 20]
banding = ["cv", "tde_bb", "AGN", "T_star", "logarithmic",
           "logarithmic", "logarithmic", "logarithmic"]

# number of photons
nphotons = 100000

# location of cloudy data
cloudy_dir = "cld_data_noce"

# version of sirocco
version = "1.0"

# number of points to sample
npts = 61

# if you want to add flags for sirocco
flags = ""


for i in range(8):
    # run the sirocco model
    sirocco_dir = run.run_multi_sirocco(log_tmin=3.5, log_tmax=6, npts=npts, version=version, flags="", macro_atoms=False,
                                        data_string="standard80_noce", data_folder="./", banding=banding[i], nbands=nbands[i], nphotons=nphotons)

    # summarize the data
    summarize_bb.summarize(sirocco_dir)

    # Make plots compared to cloudy
    plot.summary_plot(sirocco_dir, cloudy_dir)
