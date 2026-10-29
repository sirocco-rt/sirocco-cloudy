## Thin shell tests

These scripts generalise `bb_cloudy_shell`. A single optically thin shell, at a fixed electron temperature and hydrogen density, is illuminated by a chosen SED. A grid of Sirocco models loops over one parameter:

| `--loop` | varies | source normalisation |
|---|---|---|
| `t_r` | blackbody temperature | fixed radius `--rstar` (the original `bb_cloudy_shell` test) |
| `U` | ionization parameter U | set from U for the chosen SED |
| `xi` | ionization parameter ξ | set from ξ for the chosen SED |

The SED types (`--sed`) are `bb`, `power`, `cloudy` (Cloudy's "table power law"), `brem` and `model`. A `model` SED is read from a file in Sirocco's model format (`--sed_file`): two columns, wavelength in Å and F_λ in any units, such as the files in `seds/model_files` of the [sirocco-rt/data](https://github.com/sirocco-rt/data) repository. It is normalised by its 2–10 keV luminosity. The file is copied into the grid folder with a models list file pointing to it, and photon energies are limited to the wavelength range it covers. Sirocco is run with `-f`, so the electron temperature stays at `--te`.

The definitions follow Sirocco's estimators. The shell radius `--rshell` (cm) is fixed.

* U = Q(13.6 eV to the high energy limit) / (4π r² n_H c)
* ξ = L(13.6 eV to 13.6 keV) / (n_H r²)
* For a blackbody, the source normalisation is the central object radius. For the other SEDs it is the 2–10 keV luminosity, which is how Sirocco normalises them.
* For a blackbody at fixed temperature, U can only go up to the value reached when the source radius equals the shell radius.

### Files
* `sed.py`: SED shapes, conversions between source normalisation and U or ξ, n_H to shell mdot.
* `run_shell.py`: writes the parameter files and runs Sirocco and windsave2table for a grid (command line or `run_grid`).
* `collate.py`: collects a grid into `<folder>/summary.ecsv`, checks Sirocco's U, ξ, n_H and T_e against the inputs, and plots ion fractions, optionally against Cloudy.
* `example.py`: a few example grids.

### Usage
Run from this folder, with a `data` link to the Sirocco atomic data (`$SIROCCO/bin/Setup_Sirocco_Dir`). Each grid is written to its own folder, named after its parameters.

```
python run_shell.py --version 1.4 --loop U --lmin -4 --lmax 3 --npts 29 --sed power --alpha -0.9 --nproc 4
python collate.py U_power_alpha-0.9_v1.4_basic_standard80_noce_log10_np1e+05_te1e+04_nh1e+07
python run_shell.py --version 1.4 --loop t_r --npts 61
python collate.py t_r_bb_v1.4_basic_standard80_noce_log10_np1e+05_te1e+04_nh1e+07 --cloudy ../bb_cloudy_shell/cld_data_noce
python run_shell.py --version 1.4 --loop U --lmin -4 --lmax 3 --npts 29 --sed model --sed_file ../../data/seds/model_files/mean_qso.dat
```

`python run_shell.py -h` lists all options. The masterfile `standard80_noce.dat` turns off charge exchange, as in `bb_cloudy_shell`.
