"""Step 2A -- Rescreening, fully 3D-periodic route.   Eq. (5) / (6) / (13)

Screens the bare defect charge with the dielectric function of a *large*
supercell that includes the free carriers:

    W = v eps^-1 = v (1 - v chi0)^-1              Eq. (6)
    Delta V_H(q+G) = sum_G' W(q+G, q+G') Delta rho_b(q+G')   Eq. (13)

This is the route used for the FIG. 2(d) 9x9 validation.  It keeps periodic
boundary conditions in all three directions, so the vacuum spacing must be
large enough (or the 2D-truncated kernel used) to suppress image interaction.
For the mobility numbers in the main text see ``03_rescreen_quasi2d.py``.

Reproduces ``notebooks/mos2_test.ipynb`` cells 6-12.
"""

import numpy as np

from potcorr import PotCorr, cells

cell = cells.mos2_12to24

pc = PotCorr(cell)
pc.fft_init()
pc.get_vcoul()

# bare charge produced by 01_descreen.py, on the large grid
pc.rho_bare_r = np.load(cell['folder'] + 'rho_bare_12-24.npy')
pc.rho_bare_k = np.fft.fftn(pc.rho_bare_r)
pc.rho2pot_bare(kernel='3d')

# chi0 of the large, *doped* cell (chi0_b + chi0_c)
pc.read_chi(cell['folder'] + '../chi_24x24_fermi_-0.15/chimat.h5',
            cell['folder'] + '../chi_24x24_fermi_-0.15/chi0mat.h5')
pc.epsmat_init()
pc.get_epsmat_from_chi(G_ind_cut=600, kernel='3d')
pc.epsmat_inv(G_ind_cut=600)
del pc.Chi0, pc.Chi1

# Delta V_bare -> Delta V_H, i.e. apply W
pc.pot_bare2tot(kernel='3d', G_ind_cut=100)

pc.write_xsf(ftype='pot_tot', filedir=pc.folder + 'pot_rescreened.xsf')


# ---------------------------------------------------------------------------
# Symmetry-accelerated variant
# ---------------------------------------------------------------------------
# For production runs the eps^-1 matrices are computed only in the irreducible
# BZ and unfolded.  The equivalent chain is:
#
#   k_symmetry_map = pc.get_k_symmetry_map()
#   epsym_dict     = pc.get_epsym_dict(k_symmetry_map, ecut=10)
#   pc.epsmat_irrbz_init(k_symmetry_map, epsym_dict)
#   pc.get_epsinv_mat_irrbz(k_symmetry_map, epsym_dict, G_ind_cut=1000)
#   pc.get_wcoul_irrbz(k_symmetry_map, epsym_dict, G_ind_cut=1000, kernel='2d')
#   pot_k = pc.rho_ext2pot_tot_irrbz(epsym_dict, pc.rho_bare_k, G_ind_cut=1300)
#
# or, via the cached dictionary form,
#   gen_wcoul_mat_dict() -> gen_phi_G_dict() -> map_phi_G()
# with read_epsym_dict / read_epsmat_dict / read_wcoul_mat_dict reloading the
# pickles instead of recomputing them.  See notebooks/check_doped.ipynb.
