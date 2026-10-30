## Thin shell tests

These are the driving scripts for the thin shell tests; each test setup has its own folder (e.g. `../bb_cloudy_shell`). A single optically thin shell, at a fixed electron temperature and hydrogen density, is illuminated by a chosen SED. A grid of Sirocco models loops over one parameter:

| `--loop` | varies | source normalisation |
|---|---|---|
| `t_r` | blackbody temperature | fixed radius `--rstar` (the original `bb_cloudy_shell` test) |
| `U` | ionization parameter U | set from U for the chosen SED |
| `xi` | ionization parameter ξ | set from ξ for the chosen SED |

The SED types (`--sed`) are `bb`, `power`, `cloudy` (Cloudy's "table power law"), `brem` and `model`. For comparisons with Cloudy, use `cloudy` rather than `power`: it is the same SED in both codes. Its default breaks, 0.124 eV and 50 keV, are Cloudy's defaults, and Sirocco needs `cloudy_test` photon sampling for it, which `run_grid` sets. A `model` SED is read from a file in Sirocco's model format (`--sed_file`): two columns, wavelength in Å and F_λ in any units, such as the files in `seds/model_files` of the [sirocco-rt/data](https://github.com/sirocco-rt/data) repository. It is normalised by its 2–10 keV luminosity. The file is copied into the grid folder with a models list file pointing to it, and photon energies are limited to the wavelength range it covers. Sirocco is run with `-f`, so the electron temperature stays at `--te`.

The definitions follow Sirocco's estimators. The shell radius `--rshell` (cm) is fixed.

* U = Q(13.6 eV to the high energy limit) / (4π r² n_H c)
* ξ = L(13.6 eV to 13.6 keV) / (n_H r²)
* For a blackbody, the source normalisation is the central object radius. For the other SEDs it is the 2–10 keV luminosity, which is how Sirocco normalises them.
* For a blackbody at fixed temperature, U can only go up to the value reached when the source radius equals the shell radius.

### Files
* `sed.py`: SED shapes, conversions between source normalisation and U or ξ, n_H to shell mdot.
* `run_shell.py`: writes the parameter files and runs Sirocco and windsave2table for a grid (command line or `run_grid`).
* `cloudy_shell.py`: runs Cloudy for the same models (command line or `run_cloudy_grid`, which takes the same keywords as `run_grid`).
* `collate.py`: collects a grid into `<folder>/summary.ecsv` (ion fractions) and `thermal.ecsv` (T_e, heating and cooling), checks Sirocco's U, ξ, n_H and T_e against the inputs, and plots ion fractions, and T_e, heating and cooling, optionally against Cloudy.
* `example.py`: a few example grids.
* `standard80_noce.dat`: the default atomic data masterfile, which turns off charge exchange.

### Usage
The scripts can be run from any folder. The atomic data comes from a `data` link in the current folder (`$SIROCCO/bin/Setup_Sirocco_Dir`), or `$SIROCCO/xdata` if there is none. Each grid is written to its own folder, named after its parameters, inside `--outdir` (default: the current folder). From a test folder, add this folder to the Python path to use `run_shell.run_grid`, `collate.collate`, `collate.load` and `collate.plot` (see `../bb_cloudy_shell/example.py`).

```
python run_shell.py --version 1.4 --loop U --lmin -4 --lmax 3 --npts 29 --sed cloudy --alpha -0.9 --nproc 4
python collate.py U_cloudy_alpha-0.9_low0.124_hi50000_v1.4_basic_standard80_noce_cloudy_test_np1e+05_te1e+04_nh1e+07
python run_shell.py --version 1.4 --loop t_r --npts 61
python collate.py t_r_bb_v1.4_basic_standard80_noce_log10_np1e+05_te1e+04_nh1e+07 --cloudy ../bb_cloudy_shell/cld_data_noce
python run_shell.py --version 1.4 --loop U --lmin -4 --lmax 3 --npts 29 --sed model --sed_file ../../data/seds/model_files/mean_qso.dat
```

`python run_shell.py -h` lists all options.

### Thermal balance
With `fixed_te=False` (`--thermal` on the command line), T_e is found from thermal balance instead of being held fixed: Sirocco is run without `-f` and Cloudy without `constant temperature`. Sirocco changes T_e by at most about 25 per cent per ionization cycle, so each model starts from an estimate of its temperature (`te_start`, see below), and Sirocco is run with `-p` (logarithmic photon stepping) so the early cycles are cheap. The number of cycles defaults to 15 in this mode (2 for a fixed T_e). `collate.check` then reports how many models have not converged, rather than checking T_e.

`te_start` (`--te_start`) sets the starting temperature:
* `"estimate"` (default): 3×10⁴ K below U = 0.01, the Compton temperature ⟨hν⟩/4k of the SED above U = 100, and a log-log interpolation between. This is a rough guess, not a prediction.
* a number: the same starting temperature for every model.
* a folder from `run_cloudy_grid(..., fixed_te=False)`: Cloudy's T_e at the same loop value. This converges fastest, but where there is more than one equilibrium temperature it starts Sirocco near Cloudy's, which can hide a disagreement.
* `None`: `te`.

`thermal.ecsv` and `cloudy_thermal.ecsv` hold T_e and the total and component heating and cooling rates, divided by n_H² (erg cm³ s⁻¹):

| | components |
|---|---|
| heating | `phot` (photoionization, including Auger), `ff`, `comp` (Compton), `line`, `other` |
| cooling | `line`, `ff`, `rr` (recombination), `comp`, `other` |

`other` is the total minus the named components: for Sirocco's cooling mostly dielectronic recombination and collisional ionization, and for Cloudy collisional ionization, charge transfer and agents too small to be listed. Cloudy's heating and cooling agents are assigned to components by their labels (`HEAT_GROUPS` and `COOL_GROUPS` in `cloudy_shell.py`). `collate.plot_thermal` plots T_e, heating and cooling against the loop parameter (see `example.py`).

### Cloudy
`cloudy_shell.py` sets up and runs Cloudy for the same models as `run_shell.py`, from the same keywords, so a grid can be passed to both:
```
folder = run_shell.run_grid(**grid)
cloudy_folder = cloudy_shell.run_cloudy_grid(**grid)   # Sirocco-only keywords are ignored
collate.plot(collate.collate(folder), folder, cloudy_dir=cloudy_folder)
```
The command line takes the same options as `run_shell.py`, without `--version`. Each model is set up as follows.

* **Normalisation:** `ionization parameter` (log U), or `xi` (log ξ) for a ξ loop. Cloudy's ξ uses the same 1–1000 Ryd band as Sirocco. The t_r loop uses the U of each blackbody at the shell.
* **SED:** Cloudy's own continua where they match the Sirocco SED: `blackbody T`; `table power law` for `cloudy` (the same index and breaks); `bremsstrahlung T` for `brem` with α = −0.2. Any other SED (`model`, or `brem` with another α) is written once per grid as a `table SED` file over Sirocco's photon energy range. A `power` SED is refused, since Cloudy has no matching continuum.
* **Physical conditions:** `constant temperature` at T_e, `hden` n_H, and a shell 10⁻¹¹ r thick.
* **Abundances:** from the masterfile's element file, with the other elements turned off. Charge exchange is off for a `*_noce` masterfile unless `charge_exchange` says otherwise.

`run_cloudy_grid` saves the ion fractions of each element, the physical conditions, and the heating and cooling for the last zone. It collates them into `cloudy_summary.ecsv` and `cloudy_thermal.ecsv`, which have the same columns as Sirocco's `summary.ecsv` and `thermal.ecsv`. Cloudy is run as `$CLOUDY_EXE` (default `cloudy.exe`). The scripts were written against Cloudy c23.01.
