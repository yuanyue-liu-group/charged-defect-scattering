# potcorr

Descreening and rescreening of charged-defect potentials from DFT supercells.

Reference implementation for

> R. Guo, … and Y. Liu,
> *Charged defect scattering*, **Phys. Rev. Lett. 135, 126302 (2025)**.

## The problem

A charged defect in a DFT supercell carries the screening response of *that*
supercell at *that* carrier concentration. What a transport calculation needs
is the defect potential as screened by the real, dilute, doped crystal. Making
the supercell big enough is not an option — the screening length in a lightly
doped 2D semiconductor is tens of nanometres.

So the potential is corrected in two steps:

1. **Descreening** — divide out the bound-electron response of the small
   supercell to recover the bare excess charge Δρ_b.  *Eq. (11)/(12)*
2. **Rescreening** — screen Δρ_b again with the dielectric function of a large
   cell that includes the free carriers, W = (1 − vχ⁰)⁻¹ v.  *Eq. (5)/(13)*

```
        small supercell                    large supercell
   ┌────────────────────────┐        ┌────────────────────────┐
   │  Δρ_tot  (DFT)         │        │                        │
   │      │                 │        │                        │
   │   descreen  Eq. 11     │        │    rescreen  Eq. 13    │
   │      ▼                 │  ───▶  │           ▼            │
   │  Δρ_b   (bare charge)  │  zero  │  ΔV_H  (transport)     │
   └────────────────────────┘   pad  └────────────────────────┘
```

Two rescreening routes are implemented:

- **3D-periodic** (`potcorr.core`) — plane-wave (q+G) basis, periodic in all
  three directions. Used for the FIG. 2(d) validation.
- **Quasi-2D** (`potcorr.quasi2d`) — open boundaries along z, so the defect
  potential has no periodic images out of plane. **This is the route behind
  the mobility results in the main text.**

## Layout

```
potcorr/
  core.py      PotCorr — descreening + 3D-periodic rescreening (q+G basis)
  quasi2d.py   quasi-2D open-boundary rescreening, W(q_par; z, z')
  bgw.py       BerkeleyGW chi / eps^-1 HDF5 readers
  grid.py      zero-pad / coarsen between small and large supercells
  io_xsf.py    XSF readers and writers for QE densities and potentials
  cells.py     cell definitions (lattice, supercell, FFT grid, data paths)
  const.py     Ry / Bohr units
examples/      runnable drivers, one per workflow stage
notebooks/     the notebooks the paper was actually produced with
docs/METHOD.md equation-by-equation map into the code
```

## Install

```bash
pip install -e .
```

Requires `numpy`, `scipy`, `h5py`. `mpi4py` is optional — `core` picks up MPI
automatically when launched under `mpirun`.

## Use

Descreening:

```python
from potcorr import PotCorr, cells

pc = PotCorr(cells.mos2_12x12)
pc.fft_init()
pc.get_vcoul()
pc.read_rho_tot(pc.folder + 'drho.xsf')
pc.read_chi(chi1='chimat.h5', chi0='chi0mat.h5')
pc.epsmat_init()
pc.get_epsmat_from_chi(G_ind_cut=1000, kernel='3d')

pc.rho2pot_tot(kernel='3d')                   # Δρ_tot → ΔV_H,WJ
pc.pot_tot2bare(kernel='3d', G_ind_cut=1000)  # ΔV_H,WJ → ΔV_bare
pc.pot2rho_bare(kernel='3d', ncharge=1)       # ΔV_bare → Δρ_b
```

Quasi-2D rescreening:

```python
from potcorr import quasi2d

z, dist_zz = quasi2d.z_grid(half_width=12.0, nz=450)
xi, xizz   = quasi2d.load_xi('xi_Gz.npy', nGz_l=450, Lz=24.0/0.529177)
chi_rbf    = quasi2d.load_chi_carrier_rbf('chi_300k_3e13_rbf.pkl')

qs, ws  = quasi2d.build_wcoul_over_q([(Eps0, range(1)), (Eps1, range(1, 12))],
                                     dist_zz, xizz, chi_rbf,
                                     lattpara_a=3.18, Lz=24.0/0.529177)
w_interp = quasi2d.interpolate_wcoul(qs, ws, nz_pot=225)
pot_k_z  = quasi2d.pot_from_rho_bare(rho_bare_l, w_interp, fft_q_abs,
                                     lattpara=pc2.lattpara)
```

Full drivers for each stage are in `examples/`:

| Script | Stage |
|---|---|
| `01_descreen.py` | Δρ_tot → Δρ_b |
| `02_rescreen_3d.py` | Δρ_b → ΔV_H, 3D-periodic |
| `03_rescreen_quasi2d.py` | Δρ_b → ΔV_H, quasi-2D open boundary |
| `04_pot_to_rbf.py` | ΔV_H onto the dense transport k-grid |

Units are Rydberg atomic units throughout: energies in Ry, lengths in Bohr.
The one exception is that a few function arguments take Å, always documented.

## Inputs this repository does not contain

- QE ground-state runs and the `drho.xsf` density differences
- BerkeleyGW `chimat.h5` / `chi0mat.h5` (χ⁰_b) and `epsmat.h5` / `eps0mat.h5`
- the Adler–Wiser χ⁰_c pickles `chi_{T}_{n}_rbf.pkl`, and the CBM envelope
  `xi_Gz.npy` — produced by `notebooks/read_wan.ipynb` and
  `notebooks/wfc_read.ipynb`
- the scattering-matrix step, Eq. (S16)–(S20), which consumes ΔV_H

`cells.py` and the examples hard-code absolute paths into the author's Anvil
scratch tree. Edit `folder` in the relevant cell entry for your own data.

## See also

- [`../README.md`](../README.md) — the method write-up: derivation, validation
  and the full nine-stage workflow this code covers stages 3–7 of
- [EDI](https://github.com/yuanyue-liu-group/EDI) — the electron–defect
  interaction and Boltzmann transport solver that consumes ΔV_H

## License

MIT. See [LICENSE](LICENSE).
