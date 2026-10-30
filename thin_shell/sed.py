"""
SED shapes and conversions between source normalisation and ionization parameter.

Each SED is a dict with a "type" key, matching the Sirocco rad_type options:
    {"type": "bb", "t": 1e5}                                  blackbody, temperature in K
    {"type": "power", "alpha": -0.9}                          power law, L_nu ~ nu^alpha
    {"type": "cloudy", "alpha": -0.9, "low": 1, "hi": 1e4}    Cloudy "table power law": nu^2.5 below
                                                              low (eV) and nu^-2 above hi (eV)
    {"type": "brem", "t": 1.16e8, "alpha": -0.2}              L_nu ~ nu^alpha exp(-h nu / kT)
    {"type": "model", "file": "mean_qso.dat"}                 tabulated SED in Sirocco's model format:
                                                              wavelength (A) and F_lambda, any units

The source normalisation ("norm") is the central object radius in cm for a blackbody,
and the 2-10 keV luminosity in erg/s otherwise, as in Sirocco.

Definitions follow Sirocco's estimators:
    U  = Q(> 13.6 eV) / (4 pi r^2 nH c)       photons up to the high energy limit of the run
    xi = L(13.6 eV - 13.6 keV) / (nH r^2)
"""
import os
from functools import lru_cache

import numpy as np

# cgs constants, as in Sirocco's constants.h
PLANCK, BOLTZMANN, VLIGHT = 6.6262e-27, 1.38062e-16, 2.997925e10
STEFAN_BOLTZMANN, MPROT, MSOL, YR = 5.6696e-5, 1.672661e-24, 1.989e33, 3.1556925e7
HEV = 4.13620e-15  # h in eV s
XFREQMIN, XFREQMAX = 4.84e17, 2.42e18  # Sirocco's 2-10 keV band, in Hz
SED_TYPES = ("bb", "power", "cloudy", "brem", "model")


@lru_cache()
def read_model(fname):
    """wavelength (A) and F_lambda columns of a model SED file, sorted by wavelength"""
    lam, flam = np.loadtxt(fname, usecols=(0, 1), unpack=True)
    order = np.argsort(lam)
    return lam[order], flam[order]


def model_energy_range(fname):
    """the photon energy range (eV) covered by a model SED file"""
    lam = read_model(fname)[0]
    return HEV * VLIGHT * 1e8 / lam[-1], HEV * VLIGHT * 1e8 / lam[0]


def _lnu(sed, nu):
    """un-normalised L_nu"""
    kind, a = sed["type"], sed.get("alpha")
    if kind == "bb":
        # 4 pi R^2 * pi B_nu with R = 1 cm
        x = np.minimum(PLANCK * nu / (BOLTZMANN * sed["t"]), 700.0)
        return 8 * np.pi**2 * PLANCK * nu**3 / VLIGHT**2 / np.expm1(x)
    if kind == "power":
        return nu**a
    if kind == "cloudy":
        lo, hi = sed["low"] / HEV, sed["hi"] / HEV
        return np.where(nu < lo, lo**a * (nu / lo) ** 2.5,
                        np.where(nu > hi, hi**a * (nu / hi) ** -2.0, nu**a))
    if kind == "brem":
        return nu**a * np.exp(-PLANCK * nu / (BOLTZMANN * sed["t"]))
    if kind == "model":
        # interpolate F_lambda in log space as Sirocco does, zero outside the table, then L_nu ~ F_lambda lambda^2
        lam, flam = read_model(sed["file"])
        x = VLIGHT * 1e8 / nu
        f = np.exp(np.interp(np.log(x), np.log(lam), np.log(np.maximum(flam, 1e-300)), left=-np.inf, right=-np.inf))
        return f * x**2
    raise ValueError("unknown SED type {}".format(kind))


def integrate(fnu, e1, e2, npts=4000):
    """integrate fnu(nu) d nu between e1 and e2 in eV on a log grid"""
    nu = np.logspace(np.log10(e1 / HEV), np.log10(e2 / HEV), npts)
    return np.trapezoid(fnu(nu), nu)


def lnu(sed, nu, norm=1.0):
    """L_nu in erg/s/Hz: norm is R_star (cm) for a blackbody, L(2-10 keV) otherwise"""
    if sed["type"] == "bb":
        return norm**2 * _lnu(sed, nu)
    return norm * _lnu(sed, nu) / integrate(lambda f: _lnu(sed, f), XFREQMIN * HEV, XFREQMAX * HEV)


def ionization_parameters(sed, norm, r, nh, emax=1e5):
    """return (U, xi) at radius r (cm) and H density nh (cm^-3) for a source normalisation"""
    q = integrate(lambda f: lnu(sed, f, norm) / (PLANCK * f), 13.6, emax)
    lum = integrate(lambda f: lnu(sed, f, norm), 13.6, 13600.0)
    return q / (4 * np.pi * r**2 * nh * VLIGHT), lum / (nh * r**2)


def norm_from_ip(sed, value, kind, r, nh, emax=1e5):
    """the source normalisation needed to give U or xi (kind = "U" or "xi") equal to value"""
    u1, xi1 = ionization_parameters(sed, 1.0, r, nh, emax)
    ratio = value / {"U": u1, "xi": xi1}[kind]
    # the ionization parameter scales as R^2 for a blackbody and linearly with L(2-10 keV) otherwise
    return np.sqrt(ratio) if sed["type"] == "bb" else ratio


def elements(masterfile, data_dir="data"):
    """{symbol: (atomic number, abundance by number relative to H, atomic weight)} from the element file
    of a masterfile. The masterfile names the element file relative to data/, which is found in data_dir"""
    with open(masterfile) as f:
        elemfile = next(l.split()[0] for l in f if l.strip() and not l.startswith("#"))
    elemfile = os.path.join(data_dir, os.path.relpath(elemfile, "data"))
    elem = [l.split() for l in open(elemfile) if l.startswith("Element")]
    return {e[2]: (int(e[1]), 10 ** (float(e[3]) - 12), float(e[4])) for e in elem}


def rho2nh(masterfile, data_dir="data"):
    """Sirocco's conversion factor from mass density to nH, using the masterfile's element abundances"""
    elem = elements(masterfile, data_dir)
    return elem["H"][2] / (MPROT * sum(abun * weight for z, abun, weight in elem.values()))


def shell_mdot(nh, r, v, r2nh):
    """the shell mass loss rate in msol/yr giving density nh at radius r with velocity v"""
    return 4 * np.pi * r**2 * v * nh / r2nh * YR / MSOL
