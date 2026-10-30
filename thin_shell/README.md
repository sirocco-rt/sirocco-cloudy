## Thin shell tests

These are the driving scripts for the thin shell tests; each test setup has its own folder (e.g. `../bb_cloudy_shell`). A single optically thin shell, at a fixed electron temperature and hydrogen density, is illuminated by a chosen SED. A grid of Sirocco models loops over one parameter:

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
* `cloudy_shell.py`: runs Cloudy for the same models (command line or `run_cloudy_grid`, which takes the same keywords as `run_grid`).
* `collate.py`: collects a grid into `<folder>/summary.ecsv`, checks Sirocco's U, ξ, n_H and T_e against the inputs, and plots ion fractions, optionally against Cloudy.
* `example.py`: a few example grids.
* `standard80_noce.dat`: the default atomic data masterfile, which turns off charge exchange.

### Usage
The scripts can be run from any folder. The atomic data comes from a `data` link in the current folder (`$SIROCCO/bin/Setup_Sirocco_Dir`), or `$SIROCCO/xdata` if there is none. Each grid is written to its own folder, named after its parameters, inside `--outdir` (default: the current folder). From a test folder, add this folder to the Python path to use `run_shell.run_grid`, `collate.collate`, `collate.load` and `collate.plot` (see `../bb_cloudy_shell/example.py`).

```
python run_shell.py --version 1.4 --loop U --lmin -4 --lmax 3 --npts 29 --sed power --alpha -0.9 --nproc 4
python collate.py U_power_alpha-0.9_v1.4_basic_standard80_noce_log10_np1e+05_te1e+04_nh1e+07
python run_shell.py --version 1.4 --loop t_r --npts 61
python collate.py t_r_bb_v1.4_basic_standard80_noce_log10_np1e+05_te1e+04_nh1e+07 --cloudy ../bb_cloudy_shell/cld_data_noce
python run_shell.py --version 1.4 --loop U --lmin -4 --lmax 3 --npts 29 --sed model --sed_file ../../data/seds/model_files/mean_qso.dat
```

`python run_shell.py -h` lists all options.

### Cloudy
`cloudy_shell.py` sets up and runs Cloudy for the same models as `run_shell.py`, from the same keywords, so a grid can be passed to both:
```
folder = run_shell.run_grid(**grid)
cloudy_folder = cloudy_shell.run_cloudy_grid(**grid)   # Sirocco-only keywords are ignored
collate.plot(collate.collate(folder), folder, cloudy_dir=cloudy_folder)
```
The command line takes the same options as `run_shell.py`, without `--version`. Each model is set up as follows.

* **SED:** written as a `table SED` file over exactly Sirocco's photon energy range, with zero intensity outside it.
* **Normalisation:** `luminosity total`, at the same shell radius, so U and ξ match Sirocco's definitions.
* **Physical conditions:** `constant temperature` at T_e, `hden` n_H, and a shell 10⁻¹¹ r thick.
* **Abundances:** from the masterfile's element file, with the other elements turned off. Charge exchange is off for a `*_noce` masterfile unless `charge_exchange` says otherwise.

`run_cloudy_grid` saves the ion fractions of each element, the physical conditions, and the heating and cooling for the last zone. It collates the ion fractions into `cloudy_summary.ecsv`, which has the same columns as Sirocco's `summary.ecsv`. Cloudy is run as `$CLOUDY_EXE` (default `cloudy.exe`). The scripts were written against Cloudy c23.01.
