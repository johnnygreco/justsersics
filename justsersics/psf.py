import numpy as np
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel, Moffat2DKernel


__all__ = ['moffat_psf', 'gaussian_psf']


def moffat_psf(fwhm, pixscale, alpha=4.765, side=41, **kwargs):
    assert side % 2 != 0, 'side must be odd'
    width = fwhm / pixscale
    gamma = width / (2 * np.sqrt(2**(1/alpha) - 1))
    model = Moffat2DKernel(gamma=gamma, alpha=alpha, x_size=side, **kwargs)
    model.normalize()
    return model.array


def gaussian_psf(fwhm, pixscale, side=41, **kwargs):
    assert side % 2 != 0, 'side must be odd'
    sigma = gaussian_fwhm_to_sigma * fwhm / pixscale
    model = Gaussian2DKernel(x_stddev=sigma, y_stddev=sigma, 
                             x_size=side, **kwargs)
    model.normalize()
    return model.array
