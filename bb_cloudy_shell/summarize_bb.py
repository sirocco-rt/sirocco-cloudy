"""
Processes and summarizes data from Sirocco, including plotting

Key Functions:
- `get_roots`: Extracts and sorts file roots based on temperature values from filenames.
- `get_frac`: Reads fractional abundance data for a specific element from a file.
- `make_summary_table`: Creates a summary table combining master data and fractional abundance data for a given element.
- `summarize`: Generates summary tables for multiple ions (H, He, C, N, O, Fe).

Main Execution:
- Reads a directory path from command-line arguments.
- Summarizes ion data and generates summary tables.
- Calls a plotting function from `PlotsTemp` to visualize results.
"""
import PlotsTemp
import numpy as np
import sys
import os
import astropy.io.ascii as io
from astropy.table import Table, join, vstack
sys.path.append("../")

def get_roots(folder="./", prefix='xx'):
    # xfiles=glob('%s/%s*master.txt' % (folder, prefix))
    xfiles = [f for f in os.listdir(
        folder) if prefix in f and ".master.txt" in f]
    print(xfiles)
    xorder = []
    xroot = []
    for one in xfiles:
        foo = one.replace('.master.txt', '')
        xroot.append(foo)
        foo = foo.replace('%s_' % prefix, '')
        foo = foo.replace('t', '+')
        # foo=foo.replace('m','-')
        xorder.append(eval(foo))
    xtab = Table([xroot, xorder], names=['Root', 'Temp'])
    xtab['Temp'].format = '.3f'
    xtab.sort('Temp')
    return xtab


def get_frac(folder, root, element='H'):
    filename = '%s/%s.%s.frac.txt' % (folder, root, element)
    try:
        xtab = io.read(filename)
    except:
        print('Error: could not read %s', filename)
        return None

    xtab.remove_column('r')
    xtab.remove_column('inwind')
    return xtab


def make_summary_table(folder, roots, element='H'):

    record = []
    for one in roots:
        xtab = io.read("{}/{}.master.txt".format(folder, one))
        ytab = xtab['i', 'converge', 't_e', 't_r', 'ip', 'xi']
        ztab = ytab[ytab['i'] == 1]
        htab = get_frac(folder, one, element)
        zztab = join(ztab, htab, join_type='left')
        record.append(zztab)

    xxx = vstack(record)
    xxx['ip'].format = '.3e'
    xxx['xi'].format = '.3e'
    xxx['t_e'].format = '.3e'
    xxx['t_r'].format = '.3e'
    xxx.write('%s/sirocco_%s.txt' % (folder, element),
              format='ascii.fixed_width_two_line', overwrite=True)
    return xxx


def summarize(folder):
    """make the summary table for each ion

    Args:
        folder (string): folder containing the data
    """

    ions = ["H", "He", "C", "N", "O", "Fe"]
    root_tab = get_roots(folder)

    for ion in ions:
        make_summary_table(folder, root_tab['Root'], ion)


if __name__ == "__main__":
    sirocco_dir = sys.argv[1]
    summarize(sirocco_dir)
    cloudy_dir = "cld_data_noce"
    PlotsTemp.summary_plot(sirocco_dir, cloudy_dir)
