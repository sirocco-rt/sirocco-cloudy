# sirocco-cloudy

Scripts and data for comparing Sirocco with Cloudy: ionization state, heating and cooling, temperature, and so on. At the moment, this code is not packaged as a module but could be.

The repository has one folder of driving scripts, and one folder per test setup:
* `thin_shell`: the driving scripts. They run grids of single thin shell Sirocco models, looping over blackbody temperature or ionization parameter (U or ξ) for a choice of SED (blackbody, power law, Cloudy table power law, bremsstrahlung or a tabulated SED), at fixed electron temperature and hydrogen density. They also collate the results, check them against the inputs, and plot them.
* `bb_cloudy_shell`: a single shell illuminated by blackbodies of varying temperature, run with several photon banding strategies.

Each test folder holds a script that runs its setup, the Cloudy results, and the collated Sirocco results, so its plots can be remade without rerunning either code. The Python requirements are in `requirements.txt`.
