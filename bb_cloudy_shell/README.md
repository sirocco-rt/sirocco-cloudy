## Readme: blackbody ion banding tests

This folder contains tests in which a thin shell

It is possible to vary the number of photons, the banding strategy, the atomic data, the line transfer mode, any flags for sirocco, and the blackbody temperature. The basic workflow is to run the Sirocco models with a steering script, then summarize the results into a table before plotting against cloudy. The `example.py` script provides an illustration of how to do this. 

Each Sirocco model takes perhaps 7 seconds to run with 100,000 photons, so a grid with 61 points can be done in under 10 minutes on one core.

Please note, by default this has no charge exchange in Cloudy and Sirocco - this produces the best agreement for testing. The masterfile `standard80_noce.dat` achieves this.

## example.py 
This script runs a series of simulations using the sirocco model with different banding configurations.
It performs the following steps for each configuration:
1. Defines the banding and number of bands for each model.
2. Runs the `sirocco` model using the `run_multi_sirocco` function with specified parameters.
3. Summarizes the simulation data using the `summarize` function from the `summarize_bb` module.
4. Generates summary plots comparing the `sirocco` results to `cloudy` data using the `summary_plot` function from the `PlotsTemp` module.

## Alternative workflow outside of example.py
* Run the sirocco models with `python bb_run_cloudy_banding.py --version <sirocco_version> [options]` -- the full set of options can be seen below or in the bb_run_cloudy_banding.py docstring. 
* look at the name of the folder just created, something like `bb1.0_basic_standard80_noce_cv_np1e+05`
* summarize the data and plot using `summarize_bb.py bb1.0_basic_standard80_noce_cv_np1e+05`

## Script documentation:
also in files.

**bb_run_cloudy_banding.py:** 

```This script is designed to run multiple Sirocco models with specified parameters for blackbody temperature 
and photon sampling. It provides flexibility in configuring the models through various command-line arguments.

Usage:
    python bb_run_cloudy_banding.py --version <sirocco_version> [options]

Required Arguments:
    --version <str>          Version of the Sirocco executable to use.

Optional Arguments:
    --log_tmin <float>       Minimum logarithm of the blackbody temperature. Defaults to 3.5.
    --log_tmax <float>       Maximum logarithm of the blackbody temperature. Defaults to 6.0.
    --npts <int>             Number of temperature points to sample between log_tmin and log_tmax. Defaults to 61.
    --flags <str>            Additional flags to pass to Sirocco at runtime. Defaults to an empty string.
    --macro_atoms            Use macro atom line transfer. Defaults to False.
    --data_string <str>      Atomic data masterfile identifier. Defaults to "standard80_noce".
    --data_folder <str>      Path to the folder containing atomic data. Defaults to "data/".
    --banding <str>          Banding strategy for photon sampling. Options include "logarithmic" and "cloudy_test".
                             Defaults to "logarithmic".
    --nbands <int>           Number of bands to use for photon sampling. Defaults to 10.
    --nphotons <int>         Number of photons per cycle. Defaults to 100,000.

Example:
    python bb_run_cloudy_banding.py --version 1.0 --log_tmin 4.0 --log_tmax 5.5 --npts 50 --macro_atoms \
        --data_string "standard80" --data_folder "data/" --banding "logarithmic" --nbands 8 --nphotons 200000
```

**PlotTemps.py:** 

```Provides functions to read, process, and visualize ionization data from Cloudy and Sirocco outputs.
Functions:
1. read_cloudy: Reads and processes Cloudy output files into a table.
2. plot_one: Plots fractional abundances of ionization states for a given element.
3. make_chi_panel: Computes and plots the average ionization state ($\chi$) for multiple elements.
4. summary_plot: Generates a summary plot comparing ionization states from Sirocco and Cloudy data for various elements.
```

**summarize_bb.py:** 
```
processes and summarizes data from Sirocco, including plotting

Key Functions:
- `get_roots`: Extracts and sorts file roots based on temperature values from filenames.
- `get_frac`: Reads fractional abundance data for a specific element from a file.
- `make_summary_table`: Creates a summary table combining master data and fractional abundance data for a given element.
- `summarize`: Generates summary tables for multiple ions (H, He, C, N, O, Fe).

Main Execution:
- Reads a directory path from command-line arguments.
- Summarizes ion data and generates summary tables.
- Calls a plotting function from `PlotsTemp` to visualize results.
```
