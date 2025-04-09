import subprocess
import numpy as np
import shutil
from glob import glob
import os
import sys
from tqdm import tqdm
import alpro.parallel as parallel


def get_pf_blackbody(bb_temp, macro_atoms=True, data_string="standard80", data_folder="data/", banding="cloudy_test", nbands=5, nphotons=100000):
    """Get the parameter file string for a blackbody SED, thin shell model

    Args:
        bb_temp (float): Blackbody temperature
        macro_atoms (bool, optional): decide whether to use macro atom line transfer. Defaults to True.
        data_string (str, optional): atomic data masterfile. Defaults to "standard80".
        data_folder (str, optional): where the atomic data is stored. Defaults to "data/".
        banding (str, optional): what banding strategy to use. Defaults to "cloudy_test".
        nbands (int, optional): how many bands to use. Defaults to 5, not always used.

    Returns:
        _type_: _description_
    """
    # decide whether to use macro atoms and change the line transfer mode accordingly
    if macro_atoms:
        rt_mode = "macro_atoms_thermal_trapping"
    else:
        rt_mode = "thermal_trapping"

    base = '''System_type(star,cv,bh,agn,previous)                   bh

    ### Parameters for the Central Object
    Central_object.mass(msol)                  0.8
    Central_object.radius(cm)                  1e10
    Binary.mass_sec(msol)                           0.4
    Binary.period(hr)                             1000

    ### Parameters for the Disk (if there is one)
    Disk.type(none,flat,vertically.extended,rmin>central.obj.rad)                 none

    ### Parameters for Boundary Layer or the compact object in an X-ray Binary or AGN
    Central_object.radiation(yes,no)                  yes
    Central_object.rad_type_to_make_wind(bb,models,power,cloudy,brems,mono)               bb
    Central_object.blackbody_temp(K)           {:8.4e}
    Central_object.geometry_for_source(sphere,lamp_post,bubble)               sphere

    ### Parameters describing the various winds or coronae in the system
    Wind.number_of_components                  1
    Wind.type(SV,star,hydro,corona,kwd,homologous,shell,imported)                shell

    ### Parameters associated with photon number, cycles,ionization and radiative transfer options
    Photons_per_cycle                          {}
    Ionization_cycles                          2
    Spectrum_cycles                            0
    Wind.ionization(on.the.spot,ML93,LTE_tr,LTE_te,fixed,matrix_bb,matrix_pow,matrix_est)           matrix_pow
    Line_transfer(pure_abs,pure_scat,sing_scat,escape_prob,thermal_trapping,macro_atoms_escape_prob,macro_atoms_thermal_trapping)     {}
    Matom_transition_mode(mc_jumps,matrix) matrix
    Wind.radiation(yes,no)                           no
    Atomic_data                                {}{}.dat
    Surface.reflection.or.absorption(reflect,absorb,thermalized.rerad)              reflect
    Wind_heating.extra_processes(none,adiabatic,nonthermal,both)                 none

    ### Parameters for Domain 0
    Shell.wind_mdot(msol/yr)                   4.7e-20
    Shell.wind.radmin(cm)                      1e11
    Shell.wind.radmax(cm)                      1.00000000001e11
    Shell.wind_v_at_rmin(cm)                   1
    Shell.wind.v_at_rmax(cm)                   1.000010
    Shell.wind.acceleration_exponent           1
    Wind.t.init                                10000
    Wind.filling_factor(1=smooth,<1=clumped)   1

    ### Parameters for Reverberation Modeling (if needed)
    Reverb.type(none,photon,wind,matom)                 none

    ### Other parameters
    Photon_sampling.approach(T_star,cv,yso,AGN,tde_bb,min_max_freq,user_bands,cloudy_test,wide,logarithmic)          {}
    Photon_sampling.nbands      {}
    Photon_sampling.low_energy_limit(eV)    0.1
    Photon_sampling.high_energy_limit(eV)  100000
    '''.format(bb_temp, nphotons, rt_mode, data_folder, data_string, banding, nbands)

    return base


def run_sirocco_model(log_t=4, version='1.0', flags="", macro_atoms=True, data_string="standard80", data_folder="data/", banding="cloudy_test", nbands=5, nphotons=100000):
    """Run a Sirocco model with specified parameters.

    Args:
        log_t (int, optional): Logarithm of the blackbody temperature. Defaults to 4.
        version (str, optional): Version of the Sirocco executable to use. Defaults to '1.0'.
        flags (str, optional): any additional flags to pass to Sirocco at runtime
        macro_atoms (bool, optional): Whether to use macro atom line transfer. Defaults to True.
        data_string (str, optional): Atomic data masterfile identifier. Defaults to "standard80".
        data_folder (str, optional): Path to the folder containing atomic data. Defaults to "data/".
        banding (str, optional): Banding strategy to use for photon sampling. Defaults to "cloudy_test".
        nbands (int, optional): Number of bands to use for photon sampling. Defaults to 5.

    Returns:
        int: 0 if the model ran successfully, 1 otherwise.
    """
    bb_temp = 10.0 ** log_t
    pfname = 'xx_t{:.03f}'.format(log_t)

    base = get_pf_blackbody(bb_temp, macro_atoms=macro_atoms, data_string=data_string,
                            data_folder=data_folder, banding=banding, nbands=nbands, nphotons=nphotons)
    # print(log_t)

    f = open(pfname+'.pf', 'w')
    f.write(base)
    f.close()

    sirocco = "sirocco-{}".format(version)

    xcommand = "{} -f {} {} > {}.out.txt".format(
        sirocco, flags, pfname, pfname)
    # print(xcommand)
    subprocess.run(xcommand, shell=True, check=True)

    try:
        signame = '%s.sig' % pfname
        f = open(signame)
    except:
        print('Error: Could not open %s ' % signame)
        return 1

    lines = f.readlines()
    last_line = lines[-1]
    # print(last_line)
    if last_line.count('COMPLETE'):
        windsave2table = 'windsave2table-%s' % version
        xcommand = '%s %s >> %s.out.txt' % (windsave2table, pfname, pfname)
        # print(xcommand)
        subprocess.run(xcommand, shell=True, check=True)
        return 0
    return 1


def run_multi_sirocco(log_tmin=4, log_tmax=6, npts=10, version='1.0', flags="", macro_atoms=True, data_string="standard80", data_folder="data/", banding="cloudy_test", nbands=5,  nphotons=100000):
    """Run multiple Sirocco models over a range of blackbody temperatures.

    Args:
        log_tmin (int, optional): Minimum logarithm of the blackbody temperature. Defaults to 4.
        log_tmax (int, optional): Maximum logarithm of the blackbody temperature. Defaults to 6.
        npts (int, optional): Number of temperature points to sample between log_tmin and log_tmax. Defaults to 10.
        version (str, optional): Version of the Sirocco executable to use. Defaults to '1.0'.
        flags (str, optional): any additional flags to pass to Sirocco at runtime
        macro_atoms (bool, optional): Whether to use macro atom line transfer. Defaults to True.
        data_string (str, optional): Atomic data masterfile identifier. Defaults to "standard80".
        data_folder (str, optional): Path to the folder containing atomic data. Defaults to "data/".
        banding (str, optional): Banding strategy to use for photon sampling. Defaults to "cloudy_test".
        nbands (int, optional): Number of bands to use for photon sampling. Defaults to 5.

    Returns:
        directory (str): Directory where the results are stored.
    """

    log_ts = np.linspace(log_tmin, log_tmax, npts)
    # print(log_ip)
    i = 1

    if macro_atoms:
        matom_string = "matom"
    else:
        matom_string = "basic"

    if banding == "logarithmic":
        directory = "bb{}_{}_{}_log{:02d}_np{:d}".format(
            version, matom_string, data_string, nbands, nphotons)
    else:
        directory = "bb{}_{}_{}_{}_np{:d}".format(
            version, matom_string, data_string, banding, nphotons)
    try:
        os.mkdir(directory)
    except FileExistsError:
        pass

    print("running {}".format(directory))
    for log_t in tqdm(log_ts):
        pfname = 'xx_t{:.03f}'.format(log_t)
        ival = run_sirocco_model(log_t, version, macro_atoms=macro_atoms, data_string=data_string,
                                 data_folder=data_folder, banding=banding, nbands=nbands, nphotons=nphotons)
        if ival == 0:
            os.remove(pfname+'.wind_save')
            shutil.rmtree('diag_%s' % pfname)
            files = glob('%s*spec*' % pfname)
            # print(files)
            for one_file in files:
                os.remove(one_file)
            # print('Finished %s successfully (%d/%d)' %
            #      (pfname, i, len(log_ts)))
        else:
            print('Failed on %s  (%d/%d)' % (pfname, i, len(log_ts)))
        i += 1

    os.system("mv xx*.txt {}/".format(directory))
    os.system("rm -f *.out.pf *.sig *.old logfile")
    os.system("mv xx*.pf {}/".format(directory))

    return (directory)


if __name__ == "__main__":
    npts = 5
    version = sys.argv[1]
    run_multi_sirocco(log_tmin=3.5, log_tmax=6, npts=npts, version=version, macro_atoms=False,
                      data_string="standard80_noce", data_folder="data/", banding="logarithmic", nbands=10, nphotons=100000)
