from __future__ import print_function, division

import numpy as np
import skimage.measure
from astropy.utils import lazyproperty
from astropy.modeling.models import Sersic2D
from astropy.modeling.functional_models import Planar2D
from astropy.modeling import fitting
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.convolution import convolve_fft
from photutils import EllipticalAperture
from .log import logger
DEFAULT_BOUNDS = dict(r_eff=(0, np.inf), n=(0.001, np.inf))


__all__ = ['PSFConvolvedSersic2D', 'SersicFit']


class PSFConvolvedSersic2D(Sersic2D):
    """
    PSF convolved astropy 2D Sersic model.
    Parameters
    ----------
    psf : ndarray
        Image psf
    """
    
    def __init__(self, psf, amplitude=1, r_eff=1, n=4, x_0=0, y_0=0, 
                 ellip=0, theta=0, **kwargs):
        super().__init__(amplitude, r_eff, n, x_0, y_0, ellip, theta, **kwargs)
        self.sersic_deconvolved = Sersic2D(amplitude, r_eff, n, x_0, y_0, 
                                           ellip, theta, **kwargs)
        psf /= psf.sum()
        self.psf = psf
    
    def evaluate(self, x, y, amplitude, r_eff, n, x_0, y_0, ellip, theta):
        sersic = self.sersic_deconvolved.evaluate(
            x, y, amplitude, r_eff, n, x_0, y_0, ellip, theta)
        sersic_conv = convolve_fft(sersic, self.psf, boundary='wrap', 
                                   normalize_kernel=True)
        return sersic_conv


class SersicFit(object):
    """
    Class for calculating morphological properties of a single source.
    """

    sersic_params = ['amplitude', 'ellip', 'n', 'r_eff', 'theta', 'x_0', 'y_0']
    tilted_plane_params = ['slope_x', 'slope_y', 'intercept']
    
    def __init__(self, image, psf=None, mask=None):
        
        self.image = image.astype(np.float64)
        self.psf = psf
        if mask is None:
            mask = np.zeros_like(image, dtype=bool)
        self.mask = mask
        self.masked_image = np.ma.masked_array(self.image, mask)
        self._sersic = None

    def set_psf(self, psf):
        self.psf = psf

    @lazyproperty
    def diagonal_length(self):
        dy, dx = self.image.shape
        return np.sqrt(dx**2 + dy**2)
    
    @lazyproperty
    def centroid(self):
        M = skimage.measure.moments(self.masked_image.filled(0), order=1)
        xc = M[0, 1] / M[0, 0]
        yc = M[1, 0] / M[0, 0]
        return np.array([xc, yc])

    @lazyproperty
    def covariance(self):
        xc, yc = self.centroid
        M = skimage.measure.moments_central(
            self.masked_image.filled(0), center=(yc, xc), order=2)
        cov = np.array([
            [M[0, 2], M[1, 1]],
            [M[1, 1], M[2, 0]]]
        ) / M[0, 0]
        return cov

    @lazyproperty
    def eigvals(self):
        eigvals = np.linalg.eigvals(self.covariance)
        eigvals = np.sort(np.abs(eigvals))[::-1] 
        return eigvals
    
    @lazyproperty
    def semimajor_axis(self):
        return np.sqrt(self.eigvals[0])
    
    @lazyproperty
    def semiminor_axis(self):
        return np.sqrt(self.eigvals[1])
    
    @lazyproperty
    def ellipticity(self):
        return 1 - self.semiminor_axis / self.semimajor_axis
    
    @lazyproperty
    def theta(self):
        """
        Angle with respect to the positive x-axis in degrees. 
        It is wrapped to be in [0 deg, 180 deg].
        """
        cov = self.covariance
        # find principal components and rotation angle of ellipse
        sigma_x2 = cov[0, 0]
        sigma_y2 = cov[1, 1]
        sigma_xy = cov[0, 1]
        theta = 0.5 * np.arctan2(2 * sigma_xy, (sigma_x2 - sigma_y2))
        wrapped = np.rad2deg(theta) % 360.0
        wrapped = wrapped - 180 * (wrapped > 180)
        return wrapped

    @property
    def residual(self):
        return self.image - self.model

    @property
    def image_init(self):
        x_0, y_0 = self.centroid
        s_ = np.s_[int(y_0) - 5:int(y_0) + 5, 
                   int(x_0) - 5:int(x_0) + 5]
        image_init = dict(
            x_0 = self.centroid[0],
            y_0 = self.centroid[1],
            r_eff = 1.5 * self.semimajor_axis,
            theta = np.deg2rad(self.theta),
            ellip = self.ellipticity,
            amplitude = 0.3 * self.image[s_].mean(),
            n = 1.0
        )
        return image_init

    @property
    def sersic(self):
        if self._sersic is None:
            image = self.model
        else:
            image = self._sersic
        return image

    def elliptical_mask(self, scale=2, use_sersic_model=False):
        if use_sersic_model:
            pars = self.fit_params
            ell = EllipticalAperture(
                [pars['x_0'], pars['y_0']], scale * pars['r_eff'],
                scale * pars['r_eff'] * (1 - pars['ellip']),  
                np.deg2rad(pars['theta'])
            )
        else:
            ell = EllipticalAperture(
                self.centroid, scale * self.semimajor_axis,
                scale * self.semiminor_axis,  np.deg2rad(self.theta)
            )
        ell_mask = ell.to_mask().to_image(self.image.shape).astype(bool)
        return ell_mask

    def fit(self, bounds={}, fixed={}, mask=None, psf=None, 
            fitter=LevMarLSQFitter, fitter_kw={}, tilted_plane=False, 
            **init_params): 
        """
        Fit 2D PSF-convoled Sersic function to image.
        """

        if self.psf is not None:
            psf = self.psf
        elif psf is None:
            raise Exception('Must provide PSF!')

        _fixed = fixed.copy()
        _bounds = DEFAULT_BOUNDS.copy()
        for k, v in bounds.items():
            _bounds[k] = v

        tilted_plane_fixed = {}
        tilted_plane_bounds = {}
        for p in self.tilted_plane_params:
            if p in _bounds.keys():
                tilted_plane_bounds[p] = _bounds.pop(p)
            if p in _fixed.keys():
                tilted_plane_fixed[p] = _fixed.pop(p)

        tilted_plane_init = {}
        init = self.image_init.copy()
        for k, v in init_params.items():
            if k in self.tilted_plane_params:
                tilted_plane_init[k] = v 
            elif k not in self.image_init.keys():
                raise Exception(f'{k} is not a valid initial parameter')
            else:
                init[k] = v

        model_init = PSFConvolvedSersic2D(psf, bounds=_bounds, 
                                          fixed=_fixed, **init)

        if tilted_plane:
            intercept = self.image_init['amplitude'] / 10
            kw = dict(slope_x=0, slope_y=0, intercept=intercept, 
                      bounds=tilted_plane_bounds, fixed=tilted_plane_fixed)
            for k, v in tilted_plane_init.items():
                kw[k] = v
            self.plane_init = Planar2D(**kw)
            self.model_init = model_init + self.plane_init
        else:
            self.model_init = model_init

        ny, nx = self.image.shape
        yy, xx = np.mgrid[:ny, :nx]

        if mask is not None:
            masked_image = self.masked_image.copy() 
            masked_image.mask |= mask.astype(bool)
        else:
            masked_image = self.masked_image

        if type(fitter) == fitting._FitterMeta:
            fitter = fitter()
        sersic_fit = fitter(self.model_init, xx, yy, 
                            masked_image, **fitter_kw)
        self.sersic_fitter = sersic_fit

        params = {}
        for par, val in zip(sersic_fit.param_names, sersic_fit.parameters):
            if tilted_plane:
                par = par[:-2]
            params[par] = val
        params['theta'] = np.rad2deg(params['theta'])

        if params['ellip'] < 0:
            logger.debug('ell < 0: flipping so that ell > 0')
            a = (1.0 - params['ellip']) * params['r_eff']
            b = params['r_eff']
            params['ellip'] = 1.0 - b/a
            params['r_eff'] = a
            params['theta'] -= 90.0

        if (params['theta'] > 180) or (params['theta'] < 0):
            msg = 'wrapping theta = {:.2f} within 0 < theta < 180'.\
                  format(params['theta'])
            logger.debug(msg)
            wrapped = params['theta'] % 360.0
            wrapped = wrapped - 180 * (wrapped > 180)
            params['theta'] = wrapped

        model_image = sersic_fit(xx, yy)
        if tilted_plane:
            plane = Planar2D(params['slope_x'], params['slope_y'], 
                             params['intercept'])
            self.plane = plane(xx, yy)
            self._sersic = model_image - self.plane
        else:
            self._sersic = None
        self.fit_params = params
        self.model = model_image
        self.fit_info = fitter.fit_info
