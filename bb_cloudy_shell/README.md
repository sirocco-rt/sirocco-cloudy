## Blackbody thin shell test

A single thin shell is illuminated by a blackbody of radius 10¹⁰ cm at 10¹¹ cm. The blackbody temperature runs from 10^3.5 to 10^6 K, with n_H = 10⁷ cm⁻³ and T_e fixed at 10⁴ K. The test compares Sirocco's ion fractions with Cloudy's for several photon banding strategies.

There is no charge exchange in either code, since that gives the best agreement for testing. The default masterfile, `thin_shell/standard80_noce.dat`, turns it off.

The Sirocco runs use the scripts in `../thin_shell` (see the README there).

### Files
* `example.py`: runs the Sirocco grids for each banding strategy, collates and checks them, and plots them against Cloudy.
* `cld_data_noce/`: the Cloudy results (`cl_<element>_bb.dat`, first column log T_bb) and the Cloudy driver `cloudy_bb_loop.py`.
* `sirocco/`: one folder per Sirocco grid, made by `example.py`. Only `summary.ecsv` and `inputs.json` are kept under version control.
* `figures/`: the comparison plots.

### Usage
```
python example.py           # run the Sirocco grids, collate them into sirocco/, and plot
python example.py --plot    # remake the plots from the summaries in sirocco/
```
Running needs `sirocco-<version>` and `windsave2table-<version>` on the path, and the atomic data, from a `data` link here or `$SIROCCO/xdata`. Each Sirocco model takes roughly 5–25 s with 10⁵ photons; `nproc` in `example.py` sets how many run at once.
