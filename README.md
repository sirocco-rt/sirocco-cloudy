# sirocco-cloudy

This repository will contain separate folders, each of which contain scripts (and in some cases data) for running comparisons between Sirocco and Cloudy. The script can be used for running both pieces of software, and making plots to compare ionization state, heating and cooling, temperature, and so on. At the moment, this code is not packaged as a module but could be. 

Each folder should contain scripts for running Sirocco and plotting/analysing the results, as well as ancillary data. Each folder should contain a README and/or an example script. 

Currently, the included tests are:
* bb_cloudy_shell: thin shell tests of Sirocco where a single shell is illuminated by a varying temperature blackbody. Includes the capability to vary the number
* thin_shell: generalised thin shell tests, looping over blackbody temperature or ionization parameter (U or xi) for a choice of SED (blackbody, power law, Cloudy table power law, bremsstrahlung), at fixed electron temperature and hydrogen density.
