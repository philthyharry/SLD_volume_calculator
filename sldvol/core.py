import logging

import numpy as np

from sldvol.tools import find_volumes


log = logging.getLogger(__name__)


class SLDToolSetup:
    """Class for setting up and running sld->vol calculations"""

    def __init__(self, sample_x=None, sample_y=None, sample_err=None):

        self.solvent_y = [-0.56, 2.07, 6.35]
        self.solvent_x = [0.0, 38.0, 100.0]
        self.protein_y = [1.68, 2.3678, 3.49]
        self.protein_x = [0.0, 38.0, 100.0]
        self.lipid_y = [-0.39, -0.39, -0.39]
        self.lipid_x = [0.0, 38.0, 100.0]
        self.solvent_fit = np.array([0.0691, -0.56])
        self.protein_fit = np.array([0.0181, 1.68])
        self.lipid_fit = np.array([0.0, -0.390])
        self._sample_x = sample_x or [0, 0, 0]
        self._sample_y = sample_y or [0, 0, 0]
        self._sample_err = sample_err or [0, 0, 0]
        self.vs = None
        self.vp = None
        self.vl = None
        self.vl_samples = None
        self.vp_samples = None
        self.vs_samples = None

    @property
    def sample_x(self):
        return self._sample_x

    @sample_x.setter
    def sample_x(self, newx):
        self._sample_x = newx

    @property
    def sample_y(self):
        return self._sample_y

    @sample_y.setter
    def sample_y(self, newy):
        self._sample_y = newy

    @property
    def sample_err(self):
        return self._sample_err

    @sample_err.setter
    def sample_err(self, err):
        self._sample_err = err

    def find_volumes(self):
        """Calculates percentage fractions for solvent, protein and lipids in a single layer from
        the standard setup curves and sample measurements (slds at different d2o concentrations.
        """
        log.info('Calculating volumes ...')
        result, volume_samples = find_volumes(
            sample_x=self._sample_x,
            sample_y=self._sample_y,
            sample_err=self._sample_err,
            Sx = self.solvent_x,
            Sy=self.solvent_y,
            Px=self.protein_x,
            Py=self.protein_y,
            Lx=self.lipid_x,
            Ly=self.lipid_y,
            N=10000,
            return_samples=True)

        log.info('Finished calculating volumes (values are stored as vs, vp and vl for '\
                 'solvent, protein and lipids, respectively).')
        self.vs, self.vp, self.vl = result
        self.vs_samples, self.vp_samples, self.vl_samples = volume_samples

        log.info(
            'Result: \n\tSOLVENT = {:0.1f} +/- {:0.1f}%, \n\tPROTEIN = {:0.1f} +/- {:0.1f}%, '
            '\n\t  LIPID = {:0.1f} +/- {:0.1f}%'.format(
            *[round(v, 3) for v in [*result[0], *result[1], *result[2]]])
        )

    def show_plot(self):
        pass

    def show_bootstrapped_volumes(self):
        pass
