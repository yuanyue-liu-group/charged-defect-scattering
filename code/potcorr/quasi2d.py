"""Quasi-2D open-boundary rescreening: W(q_par, z, z') and Delta V_H(q_par, z).

This is the screening route actually used for the mobility results in the main
text of R. Guo et al., Phys. Rev. Lett. 135, 126302 (2025) -- Eqs. (S7),
(S11), (S12), (S13) of the Supplemental Material.  It differs from the fully
3D-periodic route in :mod:`potcorr.core` in that the out-of-plane direction is
kept in real space with open boundary conditions, so no artificial periodic
images of the defect potential appear along z:

    v(q_par; z,z')      = 2 pi exp(-q_par |z-z'|) / q_par            Eq. (S7)
    chi0(q_par; z,z')   = chi0_b(q_par; z,z') + chi0_c(q_par; z,z')  Eq. (14)
    eps(q_par; z,z')    = delta_zz' - v chi0                         Eq. (S11)
    W(q_par; z,z')      = eps^-1 v                                   Eq. (S12)
    Delta V_H(q_par, z) = int dz' W(q_par; z,z') Delta rho_b(q_par, z')
                                                                     Eq. (S13)

chi0_b is the bound-electron polarizability from BerkeleyGW, back-transformed
from (G_z, G_z') to (z, z').  chi0_c is the free-carrier part, built as a
separable product chi(q_par) xi(z) xi(z') of an Adler-Wiser in-plane
polarizability and the CBM envelope xi(z) (Eqs. S22, S24, S27).

Everything here was previously inline in ``2d_screening_chi.ipynb``; the
numerics are reproduced verbatim.

Units: Rydberg atomic units (lengths in Bohr, energies in Ry), except where a
function argument is explicitly documented as angstrom.
"""

import pickle

import numpy as np

from . import const

__all__ = [
    "bgw_chi_GzGz",
    "pad_chi_GzGz",
    "chi_bound_zz",
    "chi_carrier_zz",
    "chi_2DEG_analytic",
    "fermi_wavevector",
    "load_chi_carrier_rbf",
    "load_xi",
    "z_grid",
    "vcoul_2d",
    "epsilon_2d",
    "wcoul_2d",
    "build_wcoul_over_q",
    "interpolate_wcoul",
    "pot_from_rho_bare",
]


# ---------------------------------------------------------------------------
# chi0_b : BerkeleyGW (G_z, G_z') -> (z, z')
# ---------------------------------------------------------------------------

def bgw_chi_GzGz(Chi, qind, Gy=0, nGz=61):
    """Pull the chi0_b(q; G_z, G_z') sub-block out of a BerkeleyGW chi matrix.

    Selects the G-vectors with in-plane components (0, `Gy`) and all G_z in
    the FFT-ordered window of `nGz` points, and returns the dense complex
    ``(nGz, nGz)`` block.

    Parameters
    ----------
    Chi : potcorr.bgw.Epsmat
        Already read via ``read_epsh5()`` (chimat.h5 or chi0mat.h5).
    qind : int
        Index into ``Chi.qpts``.
    Gy : int
        In-plane G component to fix (the paper uses the G_par = 0 block plus
        a few umklapp shells, selected by the caller).
    nGz : int
        Number of G_z points; must not exceed what the BGW cutoff contains.
    """
    nmtx_ = Chi.nmtx
    qind_ = qind

    Gz_list = np.fft.fftfreq(nGz, d=1 / nGz)
    Gz_list2 = -np.fft.fftfreq(nGz, d=1 / nGz)

    Gind_list = []
    Gind_list2 = []
    for i in Gz_list:
        G_vec_ = tuple((0, Gy, i))
        G_ind_rho = Chi.G_vec2ind[G_vec_]
        G_ind_eps = Chi.gind_rho2eps[qind_][G_ind_rho]
        Gind_list.append(G_ind_eps)

    for i in Gz_list2:
        G_vec_ = tuple((0, Gy, i))
        G_ind_rho = Chi.G_vec2ind[G_vec_]
        G_ind_eps = Chi.gind_rho2eps[qind_][G_ind_rho]
        Gind_list2.append(G_ind_eps)

    chi_GzGz = np.zeros([len(Gind_list), len(Gind_list)], dtype=complex)
    for i, g0 in enumerate(Gind_list):
        for j, g1 in enumerate(Gind_list2):
            chi_real = Chi.mat[qind_, 0, 0, g1 - 1, g0 - 1, 0]
            chi_imag = Chi.mat[qind_, 0, 0, g1 - 1, g0 - 1, 1]
            chi_GzGz[i, j] = chi_real + 1j * chi_imag
    return chi_GzGz


def pad_chi_GzGz(chi_GzGz, nGz, nGz_l):
    """Zero-pad a ``(nGz, nGz)`` G_z block onto a finer ``(nGz_l, nGz_l)`` grid.

    FFT-ordered quadrant copy: low-|G_z| corners are preserved and the new
    high-|G_z| components are set to zero, i.e. the z-resolution is refined
    without inventing content.
    """
    chi_GzGz_l = np.zeros([nGz_l, nGz_l], dtype=complex)
    lo = nGz // 2 + 1
    hi = -nGz // 2 + 1
    chi_GzGz_l[:lo, :lo] = chi_GzGz[:lo, :lo]
    chi_GzGz_l[:lo, hi:] = chi_GzGz[:lo, hi:]
    chi_GzGz_l[hi:, :lo] = chi_GzGz[hi:, :lo]
    chi_GzGz_l[hi:, hi:] = chi_GzGz[hi:, hi:]
    return chi_GzGz_l


def chi_bound_zz(chi_GzGz_l, Lz):
    """chi0_b(q; G_z, G_z') -> chi0_b(q; z, z') by a 2D inverse FFT.

    Parameters
    ----------
    chi_GzGz_l : (nGz_l, nGz_l) complex ndarray
        Output of :func:`pad_chi_GzGz`.
    Lz : float
        Out-of-plane period of the *primitive-cell* BerkeleyGW calculation,
        in Bohr (24 Angstrom / Bohr_R for the published MoS2 runs).

    The ``nGz_l**2 / Lz**2`` factor converts the discrete inverse DFT back to
    the continuum double Fourier transform in z and z'.
    """
    nGz_l = chi_GzGz_l.shape[0]
    return np.fft.ifftn(chi_GzGz_l) * nGz_l * nGz_l / Lz ** 2


# ---------------------------------------------------------------------------
# chi0_c : free-carrier polarizability, separable in z
# ---------------------------------------------------------------------------

def fermi_wavevector(E_f, eff_m):
    """k_F in Bohr^-1 from a Fermi energy in eV and a band mass in m_e."""
    return np.sqrt(2 * eff_m * E_f / 2 / const.Ry2eV)


def chi_2DEG_analytic(q, kf, eff_m, Ns=1, Nv=1):
    """Analytic zero-temperature 2DEG (Stern) polarizability, in Ry units.

    Used only as the FIG. S2(b) reference against the Lindhard/Adler-Wiser
    result; the production runs use :func:`load_chi_carrier_rbf`.
    """
    return np.where(
        q > 2 * kf,
        -Ns * Nv * eff_m / 4 / np.pi * (1 - np.sqrt(1 - 4 * kf ** 2 / q ** 2)),
        -Ns * Nv * eff_m / 4 / np.pi,
    )


def load_chi_carrier_rbf(path):
    """Load a pickled ``scipy.interpolate.RBFInterpolator`` for chi0_c(q_par).

    These are produced by the Adler-Wiser sum in ``notebooks/read_wan.ipynb``
    (Eq. 14 / S27), one file per temperature and carrier concentration, named
    e.g. ``chi_300k_3e13_rbf.pkl``.  The interpolator takes q_par in
    Angstrom^-1 and returns chi in Ry units.
    """
    with open(path, "rb") as fh:
        return pickle.load(fh)


def load_xi(xi_Gz_file, nGz_l, Lz, nGz_src=225):
    """Load the CBM envelope xi(z) and form the separable kernel xi(z) xi(z').

    Parameters
    ----------
    xi_Gz_file : str
        ``xi_Gz.npy``, the G_z-space envelope written by
        ``notebooks/wfc_read.ipynb`` (Eq. S22, xi(z) = sum_{r_par} |u_CBM|^2).
    nGz_l : int
        Target number of z points.
    Lz : float
        Out-of-plane period in Bohr.
    nGz_src : int
        Number of G_z points stored in `xi_Gz_file`.

    Returns
    -------
    xi : (nGz_l,) real ndarray
    xizz : (nGz_l, nGz_l) real ndarray
        ``np.outer(xi, xi)``, the z-dependence of chi0_c in Eq. (S24).
    """
    xi_Gz = np.load(xi_Gz_file)

    xi_Gz_l = np.zeros(nGz_l, dtype=complex)
    xi_Gz_l[: nGz_src // 2] = xi_Gz[: nGz_src // 2]
    xi_Gz_l[-nGz_src // 2:] = xi_Gz[-nGz_src // 2:]

    xi = (np.fft.ifftn(xi_Gz_l) * nGz_l / Lz).real
    return xi, np.outer(xi, xi)


def chi_carrier_zz(q_abs, chi_rbf, xizz, spin_factor=0.5):
    """chi0_c(q_par; z, z') = spin_factor * chi(q_par) * xi(z) xi(z').

    `q_abs` is in Bohr^-1; the RBF interpolator is queried in Angstrom^-1.
    `spin_factor` reproduces the 0.5 used in the published runs (the
    Adler-Wiser sum in ``read_wan.ipynb`` is spin-summed).
    """
    return float(spin_factor * chi_rbf([[q_abs * const.Bohr_R]])) * xizz


# ---------------------------------------------------------------------------
# v, eps, W in the (q_par; z, z') representation
# ---------------------------------------------------------------------------

def z_grid(half_width=12.0, nz=450):
    """Open-boundary z grid and the |z - z'| distance matrix.

    Parameters
    ----------
    half_width : float
        Half the z window, in **Angstrom**.
    nz : int
        Number of z points.

    Returns
    -------
    z : (nz,) ndarray in Bohr
    dist_zz : (nz, nz) ndarray in Bohr
    """
    z = np.linspace(-half_width, half_width, nz) / const.Bohr_R
    dist_zz = np.abs(z[:, np.newaxis] - z[np.newaxis, :])
    return z, dist_zz


def vcoul_2d(q_abs, dist_zz):
    """v(q_par; z, z') = 2 pi exp(-q_par |z-z'|) / q_par   --  Eq. (S7).

    The in-plane Fourier transform of 1/|r| for a single q_par, in Ry units
    (i.e. the e^2 = 2 convention already folded in).  Open in z: no image
    charges, unlike the ``v_coul2d`` truncated kernel in :mod:`potcorr.core`.
    """
    return 2 * np.pi * np.exp(-q_abs * dist_zz) / q_abs


def epsilon_2d(v_2d, chi_zz_total, prefactor=None, dz=None):
    """eps(q_par; z, z') = delta_zz' - prefactor * v chi0   --  Eq. (S11).

    Parameters
    ----------
    v_2d, chi_zz_total : (nz, nz) ndarray
    prefactor : float or None
        Scalar multiplying ``v @ chi``.  ``None`` means 2 pi, the value used
        for the published runs.
    dz : float or None
        If given, multiplies `prefactor` -- pass ``Lz / nz`` in Bohr to use
        the discretised integral measure of Eq. (S11) instead.
    """
    coeff = 2 * np.pi if prefactor is None else prefactor
    if dz is not None:
        coeff = coeff * dz
    return np.eye(v_2d.shape[0]) - coeff * (v_2d @ chi_zz_total)


def wcoul_2d(epsilon, v_2d):
    """W(q_par; z, z') = eps^-1 v   --  Eq. (S12)."""
    return np.linalg.inv(epsilon) @ v_2d


def build_wcoul_over_q(
    Chi_list,
    dist_zz,
    xizz,
    chi_rbf,
    lattpara_a,
    Lz,
    nGz=61,
    nGz_l=450,
    q_cutoff=3.6,
    Gy_list=(0,),
    spin_factor=0.5,
):
    """Assemble W(q_par; z, z') on the BerkeleyGW q-mesh.

    Loops over the supplied (Epsmat, q-index range) pairs, builds
    chi0 = chi0_b + chi0_c, forms eps and inverts it.  Above `q_cutoff`
    (Bohr^-1) screening is switched off (eps = I), so W falls back to the
    bare 2D Coulomb kernel -- the BGW chi is not converged out there and its
    contribution to the scattering integral is negligible.

    Parameters
    ----------
    Chi_list : sequence of (Epsmat, iterable_of_q_indices)
        Typically ``[(Eps0, range(1)), (Eps1, range(1, 12))]`` -- the q -> 0
        matrix supplies the single Gamma point, the finite-q matrix the rest.
    lattpara_a : float
        In-plane lattice constant a of the *primitive* cell, in Angstrom.
    Lz : float
        Out-of-plane period of the BGW cell, in Bohr.

    Returns
    -------
    qabs_list : list of float   -- |q_par| in Bohr^-1
    wlist : list of (nGz_l, nGz_l) complex ndarray
    """
    wlist = []
    qabs_list = []
    for Chi, qinds in Chi_list:
        qpts = Chi.qpts[:]
        for Gy in Gy_list:
            for qind in qinds:
                q_abs = np.sqrt(
                    (qpts[qind, 0]) ** 2
                    + (
                        np.sqrt(3) / 3 * (qpts[qind, 0])
                        + 2 * np.sqrt(3) / 3 * (qpts[qind, 1] + Gy)
                    ) ** 2
                ) * 2 * np.pi / lattpara_a * const.Bohr_R

                v_2d = vcoul_2d(q_abs, dist_zz)

                if q_abs < q_cutoff:
                    chi_GzGz = bgw_chi_GzGz(Chi, qind=qind, Gy=Gy, nGz=nGz)
                    chi_GzGz_l = pad_chi_GzGz(chi_GzGz, nGz, nGz_l)
                    chi_zz_l = chi_bound_zz(chi_GzGz_l, Lz)
                    chi_zz_c = chi_carrier_zz(q_abs, chi_rbf, xizz, spin_factor)
                    epsilon = epsilon_2d(v_2d, chi_zz_c + chi_zz_l)
                else:
                    epsilon = np.eye(nGz_l)

                wlist.append(wcoul_2d(epsilon, v_2d))
                qabs_list.append(q_abs)
    return qabs_list, wlist


def interpolate_wcoul(qabs_list, wlist, nz_pot):
    """Make W(q_par) continuous in q_par, one 1D interpolant per (z, z') pair.

    The BGW q-mesh is far coarser than the 720x720 in-plane FFT grid of the
    large supercell, so each matrix element is linearly interpolated in
    |q_par| before being evaluated on that grid.

    `nz_pot` (225 in the published runs) may be smaller than ``nGz_l`` (450),
    in which case only the leading ``nz_pot`` rows and columns of each W are
    kept.

    Returns
    -------
    list of list of ``scipy.interpolate.interp1d``, indexed ``[iz][jz]``.
    """
    from scipy.interpolate import interp1d

    q = np.asarray(qabs_list)
    W = np.asarray(wlist)
    return [
        [interp1d(q, W[:, i, j]) for j in range(nz_pot)]
        for i in range(nz_pot)
    ]


def pot_from_rho_bare(
    rho_bare_l,
    w_interp,
    fft_q_abs,
    lattpara,
    nz_pot=225,
    z_window=30,
    progress=False,
):
    """Delta V_H(q_par, z) = int dz' W(q_par; z,z') Delta rho_b(q_par, z')  --  Eq. (S13).

    For each output plane z, the in-plane FFT of the bare defect density at
    z' is multiplied by W(q_par; z, z') and summed over z' with the trapezoidal
    weight ``dz' = lattpara[2] / nz_pot``.

    Parameters
    ----------
    rho_bare_l : (nx, ny, nz) complex ndarray
        Bare defect charge density on the *large* supercell grid, already
        re-centred on the origin and zero-padded (see
        :func:`potcorr.grid.rho_expand`).
    w_interp : nested list
        Output of :func:`interpolate_wcoul`.
    fft_q_abs : (nx, ny) ndarray
        |q_par| for each in-plane FFT point, in Bohr^-1.  Set the (0,0)
        element to a small positive number rather than 0 -- ``vcoul_2d``
        diverges there and the interpolants are not defined at exactly 0.
    lattpara : sequence of 3 floats
        Large-supercell lattice parameters in Bohr.
    z_window : int
        Only z' planes within +/- `z_window` of the slab centre contribute;
        outside the slab Delta rho_b is numerically zero and skipping those
        planes is what makes the double loop affordable.

    Notes
    -----
    The in-plane measure ``lattpara[0]*lattpara[1]/nx/ny*sqrt(3)/2`` is the
    real-space cell area per FFT point for the hexagonal (ibrav=4) lattice;
    it converts the unnormalised ``np.fft.fftn`` into the continuum
    Delta rho_b(q_par, z').
    """
    nx, ny, _ = rho_bare_l.shape
    area_elem = lattpara[0] * lattpara[1] / nx / ny * np.sqrt(3) / 2
    dz = lattpara[2] / nz_pot

    pot_k_z = np.zeros([nx, ny, nz_pot], dtype=complex)
    for i in range(nz_pot):
        if progress:
            print(i, end=" ", flush=True)
        for j in range(nz_pot // 2 - z_window, nz_pot // 2 + z_window):
            rho_k = np.fft.fftn(rho_bare_l[:, :, j]) * area_elem
            pot_k_z[:, :, i] += w_interp[i][j](fft_q_abs) * rho_k * dz
    return pot_k_z
