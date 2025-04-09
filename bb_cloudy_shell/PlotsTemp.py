#!/usr/bin/env python
# coding: utf-8
import sys
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import roman


def read_cloudy(filename='cl_helium_PL.dat'):
    """
    Reads a Cloudy output file for a given ion and processes its data into a Table

    Args:
        filename (str, optional): The name of the file to read. Defaults to 'cl_helium_PL.dat'.

    Returns:
        astropy.table.Table: A table containing the processed data from the file, or None if the file could not be read.
    """
    xfile = '%s' % (filename)
    try:
        x = ascii.read(xfile)
    except:
        print('Error: Could not find %s' % xfile)
        return None

    names = x.colnames
    x.rename_column('col1', 't_r')
    j = 1
    while j < len(names):
        x.rename_column(names[j], 'i%02d' % (j))
        j += 1
    return x


def plot_one(ax=None, cloudy='cl_helium_PL.dat', sirocco='sirocco_he', element='He', imin=0, imax=0, ipmin=-8, ipmax=4, outname=''):
    print('Comparing: ', cloudy, sirocco)
    pfile = '%s' % (sirocco)
    try:
        xp = ascii.read(pfile)
    except:
        print('Error: Could not read %s' % pfile)
        return None

    if ax == None:
        fig, ax = plt.subplots()

    xc = read_cloudy(cloudy)
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    colname = xc.colnames
    if imin <= 0:
        i = 1
    else:
        i = imin

    if imax <= 0 or imax+1 > len(colname):
        imax = len(colname)
    else:
        imax = imax+1
    color_index = 0
    while i < imax:
        xlabel = '%s %s' % (element, roman.toRoman(i))
        ax.semilogy(xc['t_r'], xc[colname[i]], '-', lw=3,
                    alpha=0.7, color=color_cycle[color_index])
        ax.semilogy(np.log10(xp['t_r']), xp[colname[i]],
                    "-o", color=color_cycle[color_index], label=xlabel)
        i += 1
        color_index += 1
        if color_index == len(color_cycle):
            color_index = 0

    if element != "Fe":
        ax.legend()
    ax.set_xlabel('log (Radiation T)')
    ax.set_ylabel('Fractional Abundance')
    ax.set_ylim(1e-4, 1.2)

    ax.set_xlim(3, 6)
    plt.tight_layout()


def make_chi_panel(ax, sirocco_dir, cl):
    for i, ion in enumerate(["H", "He", "C", "N", "O", "Fe"]):
        sirocco_fname = '{}/sirocco_{}.txt'.format(sirocco_dir, ion)
        xp = np.genfromtxt(sirocco_fname, unpack=True, skip_header=2)
        xp2 = ascii.read(sirocco_fname)
        ion_densities = xp[6:]
        # print (ion_densities)
        states = np.arange(1, len(ion_densities)+1, 1)
        xx = (np.sum(states * ion_densities.T, axis=1)-1) / \
            (len(ion_densities)-1)
        ax.plot(np.log10(xp2["t_r"]),  xx, ls="--", c="C"+str(i))

        xc = np.genfromtxt(cl[ion], unpack=True, skip_header=2)
        cion_densities = xc[1:]
        xcip = xc[0]
        cxx = (np.sum(states * cion_densities.T, axis=1)-1) / \
            (len(cion_densities)-1)
        ax.plot(xcip,  cxx, ls='-', lw=3, alpha=0.7, c="C"+str(i), label=ion)
    ax.set_ylabel("$\chi$")
    ax.legend()
    ax.set_xlim(3, 6)


def summary_plot(sirocco_dir, cloudy_dir):
    """
    Generates a summary plot comparing ionization states from sirocco and Cloudy outputs.

    The function performs the following:
    - Reads ionization data for elements (H, He, C, N, O, Fe) from the specified directories.
    - Creates a 2x3 grid of subplots to visualize ionization states for each element.
    - For each element, plots the ionization states from both sirocco and Cloudy data.
    - Computes and plots the average ionization state ($\chi$) for all elements in the last subplot.
    - Saves the resulting figure as a PNG file in the "figures" directory, with a filename based on `sirocco_dir`.

    Args:
        sirocco_dir (str): Path to the directory containing sirocco-generated ionization data files.
        cloudy_dir (str): Path to the directory containing Cloudy-generated ionization data files.

    Returns:
        None
    """
    cl = dict()
    suffix = "bb"
    cl["H"] = '%s/cl_hydrogen_%s.dat' % (cloudy_dir, suffix)
    cl["He"] = '%s/cl_helium_%s.dat' % (cloudy_dir, suffix)
    cl["C"] = '%s/cl_carbon_%s.dat' % (cloudy_dir, suffix)
    cl["N"] = '%s/cl_nitrogen_%s.dat' % (cloudy_dir, suffix)
    cl["O"] = '%s/cl_oxygen_%s.dat' % (cloudy_dir, suffix)
    cl["Fe"] = '%s/cl_iron_%s.dat' % (cloudy_dir, suffix)

    xstart = "figures/summary_{}".format(sirocco_dir.replace('/', ''))

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(15, 9))
    for i, ion in enumerate(["H", "He", "C", "O", "Fe"]):
        ax = axes.flatten()[i]
        sirocco_fname = '{}/sirocco_{}.txt'.format(sirocco_dir, ion)

        if ion == "Fe":
            plot_one(ax, cloudy=cl[ion], sirocco=sirocco_fname,
                     element=ion, imin=1, imax=27, ipmin=-8, ipmax=4)
        else:
            plot_one(ax, cloudy=cl[ion], sirocco=sirocco_fname,
                     element=ion, imin=0, imax=0, ipmin=-8, ipmax=4)

    ax = axes.flatten()[-1]
    make_chi_panel(ax, sirocco_dir, cl)
    aa = axes.flatten()[1]
    aa.set_title(sirocco_dir, fontsize=26)

    fig.tight_layout(pad=0.1)
    fig.savefig(xstart+".png", dpi=200)
