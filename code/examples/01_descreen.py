"""Step 1 -- Descreening: Delta rho_tot -> Delta rho_b.   Eq. (11) / (12)

Strips the bound-electron screening of the small defect supercell off the DFT
density difference, recovering the bare excess charge of the defect.

The code takes the equivalent potential-space route rather than the literal
rho-space form of Eq. (12):

    Delta rho_tot --v--> Delta V_H,WJ --eps=1-v chi0_b--> Delta V_bare
                   --v^-1, +Q/Omega--> Delta rho_b

Inputs
------
  drho.xsf                 Delta rho_tot from QE (defect cell minus pristine)
  ../chi_12x12/chimat.h5   chi0_b of the primitive cell (BerkeleyGW)
  ../chi_12x12/chi0mat.h5

Output
------
  rho_bare.xsf / rho_bare_12-24.npy

Reproduces ``notebooks/mos2_test.ipynb`` cell 1.  Note that the notebook calls
``ttkw.get_epsmat(...)``, which no longer exists -- the method is now
``get_epsmat_from_chi(...)``, used below.
"""

import numpy as np

from potcorr import PotCorr, cells

cell = cells.mos2_12x12

pc = PotCorr(cell)
pc.fft_init()
pc.get_vcoul()

# Delta rho_tot from the DFT supercell
pc.read_rho_tot(pc.folder + 'drho.xsf')

# chi0_b of the primitive cell
pc.read_chi(cell['folder'] + '../chi_12x12/chimat.h5',
            cell['folder'] + '../chi_12x12/chi0mat.h5')

# map every FFT point onto a (q index, G index) pair of the BGW matrices
pc.epsmat_init()

# eps = 1 - v chi0_b
pc.get_epsmat_from_chi(G_ind_cut=1000, kernel='3d')
del pc.Chi0, pc.Chi1

# Delta V_H,WJ = v (Delta rho_tot - Q/Omega); the Q/Omega term is implicit in
# get_vcoul's v_coul[0,0,0] = 0
pc.rho2pot_tot(kernel='3d')

# Delta V_bare(G) = sum_G' eps(q,G,G') Delta V(q,G')
pc.pot_tot2bare(kernel='3d', G_ind_cut=1000)

# Delta rho_b = Delta V_bare / v, with G=0 set to Q/Omega
pc.pot2rho_bare(kernel='3d', ncharge=1)

pc.write_xsf(ftype='rho_bare', filedir=pc.folder + 'rho_bare.xsf')

# --- hand off to step 2: embed the 360x360x225 result in the 720x720 grid ---
rho0 = pc.rho_bare_r.real
rho_new = np.zeros((720, 720, 225))
rho_new[180:540, 180:540, :] = rho0
np.save(cell['folder'] + 'rho_bare_12-24.npy', rho_new)
