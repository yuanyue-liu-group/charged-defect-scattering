"""Step 3 -- Interpolate Delta V_H(q_par, z) onto the dense transport k-grid.

The rescreened potential is computed on the FFT grid of the large supercell
(720x720 in-plane).  The Boltzmann transport step needs it at arbitrary
k - k' on a much denser mesh (300x300 and beyond), so each z plane is fitted
with a linear radial-basis-function interpolator and pickled.

Real and imaginary parts are fitted separately because RBFInterpolator is
real-valued.  Points with |q_par| > 3 Bohr^-1 are dropped -- the potential is
negligible there and they only slow the fit down.

Run as::

    python 04_pot_to_rbf.py

Adapted from the original ``potk2rbf.py`` + ``job-sh``.
"""

import pickle

import numpy as np
from scipy.interpolate import RBFInterpolator

from potcorr import PotCorr, cells

POT_DIR = '/anvil/projects/x-che190065/rjguo/mos2/potential/'
OUT_DIR = POT_DIR + 'rbf_model/'
NZ = 60                      # z planes to fit
PERT_NAMES = ['0.03']        # defect-charge / perturbation labels

pert = PotCorr(cells.mos2_12to48)
pert.fft_init()

X = pert.fft_kxx[:, :, pert.fft_nz // 2].flatten()
Y = (pert.fft_kyy[:, :, pert.fft_nz // 2] * 2 * np.sqrt(3) / 3
     + np.sqrt(3) / 3 * pert.fft_kxx[:, :, pert.fft_nz // 2]).flatten()
S = np.sqrt(X ** 2 + Y ** 2)

mask = S < 3
X, Y = X[mask], Y[mask]
print("Total points for interpolation:", len(X))
points = np.vstack((X.ravel(), Y.ravel())).T

for pert_name in PERT_NAMES:
    print("Reading perturbation potential in reciprocal space")
    pertpot_k = np.load(f'{POT_DIR}2d_new/pot_k_{pert_name}.npy')

    for part in ('imag', 'real'):
        print(f"Starting RBF interpolation for the {part} part")
        models = []
        for i in range(NZ):
            print(i, part)
            values = getattr(pertpot_k[:, :, i].flatten(), part)[mask]
            models.append(RBFInterpolator(points, values, kernel='linear',
                                          epsilon=1, neighbors=50))

        out = f'{OUT_DIR}pertpot_k_{pert_name}_rbf_model_{part}.pkl'
        print("Saving RBF list ->", out)
        with open(out, 'wb') as fh:
            pickle.dump(models, fh)
