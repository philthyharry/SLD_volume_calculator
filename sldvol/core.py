import logging

import matplotlib.pyplot as plt
import numpy as np

from sldvol.tools import find_volumes


log = logging.getLogger(__name__)


class SLDToolSetup:
    """Class for setting up and running sld->vol calculations
    :param list sample_x: list of d2o concentrationsat which the experimental sample
    of interest was measured SLD values (given via sample_y; both lists should be equal
    in length)
    :param list sample_y: list of SLD values determined at different D2O concentrations
    (sld values should be provided via sample_x; both lists should be equal in length)
    """

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
        self._sample_x = sample_x
        self._sample_y = sample_y
        self._sample_err = sample_err
        self.vs = None
        self.vp = None
        self.vl = None
        self.vl_samples = None
        self.vp_samples = None
        self.vs_samples = None

        try:
            match = len(self._sample_y) == len(self._sample_x)
            if not match:
                raise ValueError(
                    'Lendths of sample_x and sample_y should be equal; '
                    'instead got {} and {} values (x and y, respectively)'.format(
                        len(self._sample_x), len(self._sample_y)))
        except TypeError:
            pass

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
        """Calculates percentage fractions for solvent, protein and lipids in a single
        layer from the standard setup curves and sample measurements (slds at different
        d2o concentrations."""
        log.info('Calculating volumes ...')
        result, volume_samples = find_volumes(
            sample_x=self._sample_x,
            sample_y=self._sample_y,
            sample_err=self._sample_err,
            Sx=self.solvent_x,
            Sy=self.solvent_y,
            Px=self.protein_x,
            Py=self.protein_y,
            Lx=self.lipid_x,
            Ly=self.lipid_y,
            N=10000,
            return_samples=True)

        log.info('Finished calculating volumes: vs, vp and vl ('
                 'solvent, protein and lipids, respectively).')
        self.vs, self.vp, self.vl = result
        self.vs_samples, self.vp_samples, self.vl_samples = volume_samples

        log.info(
            'Result: \n\n\tSOLVENT = {:0.1f} +/- {:0.1f}%, '
            '\n\tPROTEIN = {:0.1f} +/- {:0.1f}%, '
            '\n\t  LIPID = {:0.1f} +/- {:0.1f}%'.format(
                *[round(v, 3) for v in [*result[0], *result[1], *result[2]]])
        )

    def show_plot(self, **kwargs):
        """Plot the sld vs d2o concentration"""

        fig_params = {'figsize': (8, 4)}
        if kwargs:
            fig_params.update(kwargs)

        plt.figure(**fig_params)
        plt.plot(self.solvent_x, self.solvent_y, color='blue', label='Solvent'), \
        plt.plot(self.protein_x, self.protein_y, color='green', label='Protein'), \
        plt.plot(self.lipid_x, self.lipid_y, color='orange', label='Lipid')

        if self._sample_x:
            plt.plot(self._sample_x, self._sample_y, '--', color='black', label='Sample')
            plt.plot(self._sample_x, self._sample_y, 'o', color='red')
            plt.errorbar(
                self._sample_x, self._sample_y, yerr=self._sample_err,
                linestyle="None", ecolor="black", marker="None", elinewidth=1.5,
                capsize=4.5, capthick=1.5)

        plt.xlabel("D$_{2}$O concentration [%]"), \
        plt.ylabel("SLD x10$^{-6}$ [$\AA^{-2}$]"), \
        plt.title("Scattering Length Density vs D$_{2}$O concentration")
        plt.legend(loc='upper left', prop={'size': 12}, fancybox=True)

    def show_bootstrapped_volumes(self, return_fig=False, **kwargs):
        """Plot the histograms with bootstrapped distributions of sampled volumes
        (used to estimate the average values and the std. deviations)."""

        if self.vs_samples is None and self.vp_samples is None and self.vl_samples is None:
            log.warning('There is nothing to plot. Run find_volumes method first!')
            return

        fig_params = {'figsize': (13, 2.5)}
        if kwargs:
            fig_params.update(kwargs)

        fig, axes = plt.subplots(ncols=3, nrows=1, **fig_params)

        axes[0].hist(self.vs_samples, bins=80, density=True, color='grey', alpha=0.5)
        axes[0].axvline(self.vs[0], color='red', linestyle='--', alpha=0.6)
        axes[0].set_ylabel('probability')
        axes[0].set_xlabel('Volume %'), axes[0].set_title('Solvent')

        axes[1].hist(self.vp_samples, bins=80, density=True, color='grey', alpha=0.5)
        axes[1].axvline(self.vp[0], color='red', linestyle='--', alpha=0.6)
        axes[1].set_xlabel('Volume %'), axes[1].set_title('Protein')

        axes[2].hist(self.vl_samples, bins=80, density=True, color='grey', alpha=0.5)
        axes[2].axvline(self.vl[0], color='red', linestyle='--', alpha=0.6)
        axes[2].set_xlabel('Volume %'), axes[2].set_title('Lipid')

        fig.suptitle('Bootstrapped distributions of the estimated volume fractions',
                     fontsize=14, y=1.1)
        if return_fig:
            return fig
