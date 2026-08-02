"""Step 2B -- Rescreening, quasi-2D open-boundary route.  Eq. (S7)/(S11)-(S13)

This is the screening actually used for the mobility results in the main text.
The out-of-plane direction is kept in real space with open boundaries, so the
defect potential has no periodic images along z.

    v(q_par; z,z')      = 2 pi exp(-q_par|z-z'|)/q_par           Eq. (S7)
    chi0 = chi0_b(q_par; z,z') + chi(q_par) xi(z) xi(z')         Eq. (14)/(S24)
    eps                 = 1 - v chi0                             Eq. (S11)
    W                   = eps^-1 v                               Eq. (S12)
    Delta V_H(q_par,z)  = int dz' W Delta rho_b(q_par,z')        Eq. (S13)

Prerequisites, all produced elsewhere in the workflow:

  chi_300k_<n>_rbf.pkl  Adler-Wiser chi0_c(q_par) per T and carrier density,
                        from notebooks/read_wan.ipynb            Eq. (14)/(S27)
  xi_Gz.npy             CBM envelope xi(z) in G_z space,
                        from notebooks/wfc_read.ipynb            Eq. (S22)
  chimat.h5/chi0mat.h5  chi0_b of the primitive cell (BerkeleyGW)
  rho_ext.npy           Delta rho_b from 01_descreen.py

Reproduces ``notebooks/2d_screening_chi.ipynb``.

Read ``docs/METHOD.md`` (Open questions) before using the absolute magnitude
of the output: two prefactor conventions in this route are unresolved.
"""

import numpy as np

from potcorr import PotCorr, cells, const, quasi2d

# --- 1. primitive-cell setup and chi0_b ------------------------------------
DIEL = '/anvil/projects/x-che190065/rjguo/mos2/dielectric/'

pc = PotCorr(cells.mos2_unit_tt)
pc.fft_init()
pc.read_epsinv(eps1=DIEL + 'chi_for_zz/12x12_cut30/chimat.h5',
               eps0=DIEL + 'chi_for_zz/12x12_cut30/chi0mat.h5')

# --- 2. z grid and the separable free-carrier kernel ------------------------
NGZ = 61            # G_z points taken from BerkeleyGW
NGZ_L = 450         # refined z grid
LZ_ANG = 24.0       # out-of-plane period of the BGW cell, Angstrom
LZ = LZ_ANG / const.Bohr_R

z, dist_zz = quasi2d.z_grid(half_width=LZ_ANG / 2, nz=NGZ_L)
xi, xizz = quasi2d.load_xi(DIEL + 'xi/xi_Gz.npy', nGz_l=NGZ_L, Lz=LZ)
chi_rbf = quasi2d.load_chi_carrier_rbf(DIEL + 'adler-wiser/chi_300k_3e13_rbf.pkl')

# --- 3. W(q_par; z, z') on the BGW q-mesh ----------------------------------
qabs_list, wlist = quasi2d.build_wcoul_over_q(
    Chi_list=[(pc.Eps0, range(1)),        # q -> 0 matrix: the Gamma point
              (pc.Eps1, range(1, 12))],   # finite-q matrix
    dist_zz=dist_zz,
    xizz=xizz,
    chi_rbf=chi_rbf,
    lattpara_a=pc.lattpara_unit[0],
    Lz=LZ,
    nGz=NGZ,
    nGz_l=NGZ_L,
    q_cutoff=3.6,
)

# --- 4. continuous in q_par ------------------------------------------------
NZ_POT = 225        # z planes of the large-supercell density
w_interp = quasi2d.interpolate_wcoul(qabs_list, wlist, nz_pot=NZ_POT)

# --- 5. apply to the bare density on the large supercell -------------------
pc2 = PotCorr(cells.mos2_12to48)
pc2.fft_init()

rho_bare_l = np.load('/anvil/projects/x-che190065/rjguo/mos2/'
                     'potential/dft/6x6/rho_bare_large.npy')

fft_q_abs = np.sqrt(
    pc2.fft_kxx[:, :, pc2.fft_nz // 2] ** 2
    + (pc2.fft_kyy[:, :, pc2.fft_nz // 2] * 2 * np.sqrt(3) / 3
       + np.sqrt(3) / 3 * pc2.fft_kxx[:, :, pc2.fft_nz // 2]) ** 2)
fft_q_abs[0, 0] = 0.00027      # q=0 is outside the interpolation range

pot_k_z = quasi2d.pot_from_rho_bare(
    rho_bare_l, w_interp, fft_q_abs,
    lattpara=pc2.lattpara, nz_pot=NZ_POT, z_window=30, progress=True)

np.save('pot_k_quasi2d.npy', pot_k_z)
