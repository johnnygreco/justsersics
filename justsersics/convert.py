import numpy as np
from scipy.special import gammaincinv, gamma
from astropy import units as u


__all__ = [
    'mu_0_to_mu_e_ave',
    'sersic_surface_brightness',
    'mu0_to_mtot',
    'mtot_to_mu0', 
    'astro_units'
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


def r_eff_to_circ(r_e, ellip):
    q = 1 - ellip
    r_e = _check_r_e_unit(r_e)
    r_circ = r_e * np.sqrt(q)
    return r_circ


def mu_0_to_mu_e_ave(mu_0, n):
    mu_e = mu_0 + 2.5 * b_n(n) / np.log(10)
    mu_e_ave = mu_e - 2.5 * np.log10(f_n(n))
    return mu_e_ave


def sersic_surface_brightness(m_tot, r_e, n, ellip):
    r_circ = r_eff_to_circ(r_e, ellip)
    area = np.pi * r_circ.to('arcsec').value**2
    mu_e_ave = m_tot + 2.5 * np.log10(2 * area)
    mu_e = mu_e_ave + 2.5 * np.log10(f_n(n))
    mu_0 = mu_e - 2.5 * b_n(n) / np.log(10)
    return mu_0, mu_e, mu_e_ave


def mu0_to_mtot(mu_0, r_e, n, ellip):
    r_circ = r_eff_to_circ(r_e, ellip)
    area = np.pi * r_circ.to('arcsec').value**2
    mu_e = mu_0 + 2.5*b_n(n)/np.log(10) 	
    mu_e_ave = mu_e - 2.5*np.log10(f_n(n))
    m_tot = mu_e_ave - 2.5 * np.log10(2 * area)
    return m_tot


def mtot_to_mu0(m_tot, r_e, n, ellip):
    r_circ = r_eff_to_circ(r_e, ellip)
    area = np.pi * r_circ.to('arcsec').value**2
    mu_e_ave = m_tot + 2.5 * np.log10(2 * area)
    mu_e = mu_e_ave + 2.5*np.log10(f_n(n))
    mu_0 = mu_e - 2.5*b_n(n)/np.log(10)
    return mu_0


def astro_units(sersic_params, pixscale, zpt):
    I_e = sersic_params['amplitude']
    q = 1 - sersic_params['ellip']
    n = sersic_params['n']
    r_eff = sersic_params['r_eff'] * pixscale
    r_circ = r_eff * np.sqrt(q) 
    b_n = gammaincinv(2.0 * n, 0.5)
    mu_e = zpt - 2.5 * np.log10(I_e / pixscale**2)
    mu_0 = mu_e - 2.5 * b_n / np.log(10)
    f_n = gamma(2 * n) * n * np.exp(b_n) / b_n**(2 * n)
    mu_e_ave = mu_e - 2.5 * np.log10(f_n)
    A_eff = np.pi * (r_circ * pixscale)**2
    m_tot = mu_e_ave - 2.5 * np.log10(2 * A_eff)
    params = dict(
        m_tot = m_tot, 
        r_eff = r_eff, 
        mu_0 = mu_0, 
        mu_e = mu_e, 
        mu_e_ave = mu_e_ave
    )
    return params
