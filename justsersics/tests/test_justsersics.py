import numpy as np
from astropy import units as u
from astropy.modeling.functional_models import Planar2D
from ..psf import moffat_psf, gaussian_psf
from ..justsersics import SersicFit
from ..justsersics import PSFConvolvedSersic2D


yy, xx = np.mgrid[:201, :201]
theta = 45 * u.deg
sersic_params = dict(r_eff=30, x_0=100, y_0=100, n=1, ellip=0.3, theta=theta)


def _check_fit(truth, meas):
    for p, v in truth.items():
        assert np.allclose(v, meas[p])


def test_sersic_fit_with_moffat_psf():
    psf = moffat_psf(0.7, 0.2)
    sersic = PSFConvolvedSersic2D(psf, **sersic_params)
    fitter = SersicFit(sersic(xx, yy), psf=psf)
    fitter.fit()
    _check_fit(sersic_params, fitter.fit_params)


def test_sersic_fit_with_gaussian_psf():
    psf = gaussian_psf(0.7, 0.2)
    sersic = PSFConvolvedSersic2D(psf, **sersic_params)
    fitter = SersicFit(sersic(xx, yy), psf=psf)
    fitter.fit()
    _check_fit(sersic_params, fitter.fit_params)


def test_sersic_with_tilted_plane():
    psf = moffat_psf(0.7, 0.2)
    params = dict(slope_x=0.01, slope_y=0.01, intercept=0.01)
    plane = Planar2D(**params)
    model = PSFConvolvedSersic2D(psf, **sersic_params) + plane
    fitter = SersicFit(model(xx, yy), psf=psf)
    fitter.fit(tilted_plane=True)
    params.update(sersic_params)
    _check_fit(params, fitter.fit_params)
