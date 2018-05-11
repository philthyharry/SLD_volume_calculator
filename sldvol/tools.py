import logging

import numpy as np


log = logging.getLogger(__name__)


def calculate_volume_fractions(target_fit, solvent_fit, protein_fit, lipid_fit):
    """Calculates the volume fractions of a solvent, protein and lipid components in a
    sample from the experimental measurement at 3 concentrations and using the standard
    solvent, protein and lipid SLD/D20 concentration plots.

    :param float target_fit: list with slope and y-intercept from the experimental
    D2O-conc/SLD line fit
    :param list solvent_fit: list with solvent fitted slope and y-intercept values,
    eg. [0.0691, -0.56]
    :param list protein_fit: list with protein fitted slope and y-intercept values
    :param list lipid_fit: list with lipid fitted slope and y-intercept values
    """

    coeff_M = np.array([
        [solvent_fit[0], protein_fit[0], lipid_fit[0]],  # slopes
        [solvent_fit[1], protein_fit[1], lipid_fit[1]],  # y-intercepts
        [1, 1, 1]
    ])

    soltn_M = np.array([target_fit[0], target_fit[1], 1.0])

    # Solve the system of linear Equations:
    # vol[0] = solv, vol[1] = prot, vol[2] = lip
    volumes = np.linalg.solve(coeff_M, soltn_M)

    vs = volumes[0]*100
    vp = volumes[1]*100
    vl = volumes[2]*100

    return vs, vp, vl


def bootstrap_volumes_from_uncertainty(target_x, target_y, err_list,
                                       solvent_fit, protein_fit, lipid_fit, N=10000):
    """Bootstraps the sample volumes calculated from relative slopes of the curve by
    sampling different SLD values from within the uncertainty ranges.

    :param float target_x: D2O concentration values for which the target sample slds
    were measured
    :param float target_y: target sample SLD values determined experimentally
    :param list err_list: list of SLD uncertainties of the target sample
    :param list solvent_fit: list with solvent fitted values: [slope, y-intercept],
    eg. [0.0691, -0.56]
    :param list protein_fit: list with protein fitted values: [slope, y-intercept]
    :param list lipid_fit: list with lipid fitted values: [slope, y-intercept]
    :return: 2d array with Nx3 values: volume fractions of form
    [[solvent, protein, lipid], ...]
    :rtype <numpy.array(Nx3)>
    """
    samples_tmp = np.array(
        [np.random.normal(target_y[i], err_list[i]/3, N) for i in range(len(target_y))])
    samples = samples_tmp.T
    volumes = np.zeros([N, len(target_y)])

    for sn, sample in enumerate(samples):
        sample_fit = np.polyfit(target_x, sample, 1)
        volumes[sn] = calculate_volume_fractions(
            sample_fit, solvent_fit, protein_fit, lipid_fit)

    return volumes.T


def find_volumes(sample_x, sample_y, sample_err, Sx, Sy, Px, Py, Lx, Ly,
                 N=10000, return_samples=False):
    """Combines all the input data (slds and d2o concentrations) for sample, solvent,
    protein and lipid and calculates the volume fractions of the latter three in the
    sample with estimated errors."""

    solvent_fit = np.polyfit(Sx, Sy, 1)
    protein_fit = np.polyfit(Px, Py, 1)
    lipid_fit = np.polyfit(Lx, Ly, 1)

    vs, vp, vl = bootstrap_volumes_from_uncertainty(
        sample_x, sample_y, sample_err, solvent_fit, protein_fit, lipid_fit, N=N)

    vsm, vpm, vlm = vs.mean(), vp.mean(), vl.mean()

    # Check if the volumes are > zero (ie make sense):
    if vsm < 0.0 or vsm > 100.0 or vpm < 0.0 or vpm > 100.0 or vlm < 0.0 or vlm > 100.0:
        raise ValueError(
            "ERROR! Some of the calculated volume values make no sense "
            "(eg. have negative values or are larger than 100%): \n"
            "V_sol = {0:2}, V_prot = {1:2}, V_lip = {2:2} \nCheck your "
            "experimental data.".format(round(vs, 4), round(vp, 4), round(vl, 4)))

    result = list(zip([vsm, vpm, vlm], [vs.std(), vp.std(), vl.std()]))

    if return_samples:
        return result, [vs, vp, vl]
    return result
