"""
Example grids: the blackbody temperature loop of bb_cloudy_shell, then U and xi loops
for a power law and a blackbody of fixed temperature, each run with Sirocco and Cloudy
from the same keywords and plotted against each other. Needs the Sirocco atomic data (a
data link here or $SIROCCO/xdata) and Cloudy ($CLOUDY_EXE, default cloudy.exe).
"""
import cloudy_shell
import collate
import run_shell

meta = dict(version="1.4", nphotons=100000, banding="logarithmic", nbands=10, nproc=4)

grids = [
    dict(loop="t_r", lmin=3.5, lmax=6.0, npts=61, sed={"type": "bb", "t": 1e5}),
    dict(loop="U", lmin=-4.0, lmax=3.0, npts=29, sed={"type": "power", "alpha": -0.9}),
    dict(loop="xi", lmin=-2.0, lmax=5.0, npts=29, sed={"type": "power", "alpha": -0.9}),
    dict(loop="U", lmin=-4.0, lmax=3.0, npts=29, sed={"type": "bb", "t": 1e5}),
]

for grid in grids:
    folder = run_shell.run_grid(te=1e4, nh=1e7, **grid, **meta)
    cloudy_folder = cloudy_shell.run_cloudy_grid(te=1e4, nh=1e7, **grid, **meta)
    tab = collate.collate(folder)
    collate.check(tab)
    collate.plot(tab, folder, cloudy_dir=cloudy_folder)
