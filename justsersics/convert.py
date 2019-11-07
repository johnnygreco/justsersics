import numpy as np
from scipy.special import gammaincinv, gamma
from astropy import units as u


__all__ = [
    'mu_0_to_mu_e_ave',
    'sersic_surface_brightness',
    'mu0_to_mtot',
    'mtot_to_mu0'
]


def _check_r_e_unit(r_e):
    if type(r_e) != u.Quantity:
        r_e = r_e * u.arcsec
    return r_e


def b_n(n):
    return gammaincinv(2.*n, 0.5)


def f_n(n):
    _bn = b_n(n)
    return gamma(2*n)*n*np.exp(_bn)/_bn**(2*n)


def mu_0_to_mu_e_ave(mu_0, n):
    mu_e = mu_0 + 2.5*b_n(n)/np.log(10)
    mu_e_ave = mu_e - 2.5*np.log10(f_n(n))
    return mu_e_ave


def sersic_surface_brightness(m_tot, r_e, n):
    r_e = _check_r_e_unit(r_e)
    area = np.pi * r_e.to('arcsec').value**2
    mu_e_ave = m_tot + 2.5 * np.log10(2 * area)
    mu_e = mu_e_ave + 2.5*np.log10(f_n(n))
    mu_0 = mu_e - 2.5*b_n(n)/np.log(10)
    return mu_0, mu_e, mu_e_ave


def mu0_to_mtot(mu_0, r_e, n):
    r_e = _check_r_e_unit(r_e)
    area = np.pi * r_e.to('arcsec').value**2
    mu_e = mu_0 + 2.5*b_n(n)/np.log(10) 	
    mu_e_ave = mu_e - 2.5*np.log10(f_n(n))
    m_tot = mu_e_ave - 2.5 * np.log10(2 * area)
    return m_tot


def mtot_to_mu0(m_tot, r_e, n):
    r_e = _check_r_e_unit(r_e)
    area = np.pi * r_e.to('arcsec').value**2
    mu_e_ave = m_tot + 2.5 * np.log10(2 * area)
    mu_e = mu_e_ave + 2.5*np.log10(f_n(n))
    mu_0 = mu_e - 2.5*b_n(n)/np.log(10)
    return mu_0
