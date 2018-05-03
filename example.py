import logging

import numpy as np
import matplotlib.pyplot as plt

from helpers.core import find_volumes

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)


# SETUP data points for Solvent, Protein and lipid for SLD/D2O conc curves:

# Lipids:
lipid_x = [0, 38, 100]
lipid_y = [-0.39, -0.39, -0.39]
lipid_fit = np.polyfit(lipid_x, lipid_y, 1)   # L_fit = [ 0.0,  -0.390]  # tuple of (slope, intercept)

# Proteins:
protein_x = [0, 38, 100]
protein_y = [1.68, 2.3678, 3.49]
protein_fit = np.polyfit(protein_x, protein_y, 1)    # P_fit = [ 0.0181,  1.68  ]

# Solvent:
solvent_x = [0, 38, 100]
solvent_y = [-0.56, 2.07, 6.35]
solvent_fit = np.polyfit(solvent_x, solvent_y, 1)    #S_fit = [0.0691, -0.56]


if __name__ == '__main__':

    logging.info('Starting logging...')

    # User provides the SLD/D2O concentrations and sld uncertainty data:
    sample_d2o_conc = [0, 38, 100]
    sample_slds = [0.017, 1.8985, 4.961]
    sample_sld_errors = [0.08, 0.13, 0.012]

    sample_fit_parameters = np.polyfit(sample_d2o_conc, sample_slds, 1)

    # Calculate the volume fractions with uncertainties
    # and return bootstrapped volume samples for plotting:
    res, bootstrapped_vols = find_volumes(
        sample_d2o_conc, sample_slds, sample_sld_errors,
        solvent_x, solvent_y, protein_x, protein_y, lipid_x, lipid_y,
        return_samples=True)

    print("\nVolume fractions:\n\n\tSolvent: {:0.1f} +/- {:0.1f} [%]" \
          "\n\tProtein: {:0.1f} +/- {:0.1f} [%]" \
          "\n\tLipid: {:0.1f} +/- {:0.1f} [%]\n\n".format(
            res[0][0], res[0][1], res[1][0],
            res[1][1], res[2][0], res[2][1]
        )
    )

    # ----- PLOT the SLD vs [D2O]: -----

    x = np.linspace(0, 100, 3)

    plt.figure(figsize=(8, 4))
    plt.plot(sample_d2o_conc, solvent_y, color='blue', label='Solvent')
    plt.plot(sample_d2o_conc, protein_y, color='green', label='Protein')
    plt.plot(sample_d2o_conc, lipid_y, color='orange', label='Lipid')
    plt.errorbar(
        sample_d2o_conc, sample_slds, yerr=sample_sld_errors,
        linestyle="None", ecolor="black", marker="None", elinewidth=1.5,
        capsize=4.5, capthick=1.5)
    plt.plot(sample_d2o_conc, sample_slds, '--', color='black', label='Sample')
    plt.plot(sample_d2o_conc, sample_slds, 'o', color='red')
    legend = plt.legend(loc='upper left', prop={'size': 12}, fancybox=True)
    plt.xlabel("D$_{2}$O concentration [%]")
    plt.ylabel("SLD x10$^{-6}$ [$\AA^{-2}$]")
    plt.title("Scattering Length Density vs D$_{2}$O concentration")
    plt.show()
