"""potcorr -- descreening and rescreening of charged-defect potentials.

Reference implementation of the potential-correction scheme in

    R. Guo, ... and Y. Liu, "Charged defect scattering ...",
    Phys. Rev. Lett. 135, 126302 (2025).

The problem: a charged defect in a DFT supercell carries the screening
response of *that* supercell, at *that* carrier concentration.  What the
transport calculation needs is the defect potential screened by the real,
dilute, doped crystal.  The fix is two steps:

  1. **Descreening** (Eq. 11/12) -- divide out the bound-electron response of
     the small supercell to recover the bare excess charge Delta rho_b.
  2. **Rescreening** (Eq. 5/13) -- screen Delta rho_b again with the
     dielectric function of a large cell that includes the free carriers,
     W = (1 - v chi0)^-1 v.

Two rescreening routes are provided:

  * :mod:`potcorr.core` -- fully 3D-periodic, in a plane-wave (q+G) basis.
    Used for the FIG. 2(d) validation.
  * :mod:`potcorr.quasi2d` -- quasi-2D with open boundaries along z.  This is
    the route behind the mobility results in the main text.

Quick start::

    from potcorr import PotCorr, cells

    pc = PotCorr(cells.mos2_12x12)
    pc.fft_init()
    pc.get_vcoul()
    pc.read_rho_tot(pc.folder + 'drho.xsf')
    pc.read_chi(chi1='chimat.h5', chi0='chi0mat.h5')
    pc.epsmat_init()
    pc.get_epsmat_from_chi(G_ind_cut=1000, kernel='3d')
    pc.rho2pot_tot(kernel='3d')
    pc.pot_tot2bare(kernel='3d', G_ind_cut=1000)
    pc.pot2rho_bare(kernel='3d', ncharge=1)   # -> pc.rho_bare_r

See ``docs/METHOD.md`` for the equation-by-equation map into the code, and
for the two open numerical questions in the quasi-2D route.
"""

from . import const
from . import io_xsf
from . import bgw
from . import grid
from . import cells
from . import quasi2d
from .core import PotCorr, filter_G, round_tuple, progress_bar

__version__ = "0.1.0"

__all__ = [
    "PotCorr",
    "filter_G",
    "round_tuple",
    "progress_bar",
    "const",
    "io_xsf",
    "bgw",
    "grid",
    "cells",
    "quasi2d",
]
