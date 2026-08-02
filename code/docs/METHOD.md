# Equation-to-code map

Equation numbers refer to

> R. Guo *et al.*, *Phys. Rev. Lett.* **135**, 126302 (2025),
> and its Supplemental Material (`Eq. S…`).

Line numbers are for the files in this repository.

---

## 1. Descreening — Eq. (11) / (12)

> Δρ_b = Δρ_tot − χ⁰_b · v · (Δρ_tot − Q/Ω)

The code does **not** transcribe the ρ-space form of Eq. (12). It takes the
algebraically equivalent potential-space route:

```
Δρ_tot ──v──▶ ΔV_H,WJ ──ε = 1 − vχ⁰_b──▶ ΔV_bare ──v⁻¹, +Q/Ω──▶ Δρ_b
```

| Code | Role |
|---|---|
| `const.py` | Ry / Bohr units |
| `cells.py` | small-supercell entries, e.g. `mos2_6x6_z24`, `mos2_12x12` (SM "Small supercell calculations") |
| `bgw.py:Epsmat.read_epsh5` | reads χ⁰_b from the primitive cell's `chimat.h5` / `chi0mat.h5` |
| `grid.py:read_dat`, `io_xsf.py:Rho` | reads Δρ_tot from QE |
| `core.py:99` `__init__`, `core.py:139` `fft_init` | builds the q_s+G grid |
| `core.py:225` `get_vcoul` | v = 8π/\|q+G\|² (Ry). **`v_coul[0,0,0] = 0` is the implicit δ_{q_s+G′,0}·Q/Ω of Eq. (12)** |
| `core.py:256` `read_chi` | loads χ⁰_b |
| `core.py:606` `epsmat_init` | splits every FFT point into q_s+G and maps it onto BerkeleyGW's (q index, G index) — the whole index bookkeeping of the Eq. (12) sum |
| `core.py:762` `get_epsmat_from_chi(kernel='3d')` | ε = 1 − v·χ⁰_b |
| `core.py:993` `rho2pot_tot` | ΔV_H,WJ = v(Δρ_tot − Q/Ω) |
| **`core.py:1173` `pot_tot2bare`** | **ΔV_bare(G) = Σ_G′ ε(q,G,G′) ΔV(q,G′)** |
| **`core.py:1004` `pot2rho_bare(ncharge=Q)`** | **Δρ_b = ΔV_bare/v, with G=0 set to Q/Ω** |
| `core.py:1207` `write_xsf(ftype='rho_bare')` | output |

**Driver:** `examples/01_descreen.py` (from `notebooks/mos2_test.ipynb` cell 1,
`ncharge=1`). `notebooks/Na_test.ipynb` c33–c36 runs the same flow for a 3D
system with `ncharge=9`.

**Validation:** `notebooks/rho_bare.ipynb` (FIG. 2a–c),
`notebooks/rho_bare2.ipynb` (FIG. S1, z = 20–36 Å convergence).

Descreening does **not** need χ⁰_c, W, ε⁻¹, or the symmetry maps.

---

## 2. Bridging the small and large supercells

The SM's four-step transfer (FFT → real space → zero the outer region → FFT):

| Code | Role |
|---|---|
| `grid.py:rho_expand` | zero-pad (`mode='constant'`) |
| `grid.py:downsample_3d_array` | coarsen |

The notebooks also contain hand-written quadrant-shuffle versions:
`mos2_test.ipynb` c2–c4 (12×12 → 24×24) and `2d_screening_chi.ipynb` c24–c27
(6×6 → 24×24, 360 → 720 in-plane points).

---

## 3A. Rescreening — fully 3D-periodic — Eq. (5) / (6) / (13)

Used for the FIG. 2(d) 9×9 validation.

| Code | Role |
|---|---|
| `core.py:262` `read_epsinv` | reads the large cell's ε⁻¹ on its q-grid (doped, i.e. including χ⁰_c) |
| `core.py:225` `get_vcoul` (`is2d=True` → `v_coul2d`) | 2D-truncated Coulomb kernel |
| `core.py:514` `vcoul2d0modify`, `core.py:521` `vcoul0modify` | NNS correction to v as q → 0 |
| `core.py:805` `epsmat_inv`, `core.py:816` `get_epsinv_mat` | ε⁻¹ |
| **`core.py:932` `get_wcoul`** | **W = v·ε⁻¹ ← Eq. (6)** |
| **`core.py:1066` `rho_ext2pot_tot`** | **ΔV_H(q+G) = Σ_G′ W(q+G,q+G′) Δρ_b(q+G′) ← Eq. (13)** |

Symmetry acceleration (enabled in all production runs; unfolds ε⁻¹ from the
irreducible BZ to the full BZ):

`filter_G` (`core.py:46`), `get_k_symmetry_map` (`:292`), `get_epsym_dict`
(`:380`), `gen_mat_q` / `gen_mat_dict` (`:465`, `:492`),
`epsmat_irrbz_init` (`:665`), `get_epsinv_mat_irrbz` (`:840`),
`get_wcoul_irrbz` (`:865`), `rho_ext2pot_tot_irrbz` (`:1101`).

Equivalent dictionary-based chain: `gen_wcoul_mat_dict` (`:552`) →
`gen_phi_G_dict` (`:583`) → `map_phi_G` (`:592`), with
`read_epsym_dict` / `read_epsmat_dict` / `read_wcoul_mat_dict`
(`:443`, `:448`, `:458`) reloading the cached pickles.

**Drivers:** `examples/02_rescreen_3d.py`; `notebooks/check_doped.ipynb`
c1–c5, c13–c14, c20; `notebooks/mos2_test.ipynb` c6–c12.

---

## 3B. Rescreening — quasi-2D, open boundaries — Eq. (S7) / (S11)–(S13)

**This is the route behind the mobility results in the main text.** It used to
live entirely inside `2d_screening_chi.ipynb`; it is now `potcorr/quasi2d.py`.

| `quasi2d.py` | Notebook cell | Role |
|---|---|---|
| `bgw_chi_GzGz` | c4 `BGW2obz` | pulls χ(q; G_z,G_z′) out of BerkeleyGW |
| `z_grid` | c6 | z grid and \|z − z′\| matrix |
| `load_xi` | c7 | ξ(z), `xizz = outer(ξ,ξ)` — Eq. (S22)/(S24) |
| `vcoul_2d` | c11 | v = 2π e^{−q\|z−z′\|}/q — **Eq. (S7)** |
| `pad_chi_GzGz`, `chi_bound_zz` | c11 | χ⁰_b(q∥, z, z′) |
| `chi_carrier_zz` | c11 | χ⁰_c = ½ χ(q) ξ(z)ξ(z′) — **Eq. (14)/(S27)** |
| `epsilon_2d` | c11 | ε = 1 − v χ⁰ — **Eq. (S11)** |
| `wcoul_2d` | c11 | W = ε⁻¹ v — **Eq. (S12)** |
| `build_wcoul_over_q` | c11 | the whole q-loop, with the q > 3.6 Bohr⁻¹ cutoff |
| `interpolate_wcoul` | c16 | W(q∥) made continuous, one `interp1d` per (z,z′) |
| **`pot_from_rho_bare`** | **c36** | **ΔV_H(q∥,z) = ∫dz′ W Δρ_b — Eq. (S13)** |
| `chi_2DEG_analytic` | c4 `get_chi_2DEG` | analytic 2DEG χ, the FIG. S2(b) Lindhard reference |

**Driver:** `examples/03_rescreen_quasi2d.py`.

**Where χ⁰_c comes from** (Adler–Wiser, Eq. 14 / S27):

- `notebooks/read_wan.ipynb` c17 — the AW sum (f_i − f_f)/(ε_i − ε_f)·φ_overlap,
  with the Fermi level solved from the carrier concentration (c14);
  c22–c24 fit and pickle `chi_{T}_{n}_rbf.pkl`
- `notebooks/wfc_read.ipynb` c26–c31 — ξ(z) = Σ_{r∥}|u_CBM|², saved as
  `xi_z.npy` / `xi_Gz.npy` (Eq. S22)
- the wavefunction overlaps φ_{k,k+q∥+G∥} come from `get_wfc_overlap_kp_ks_npy`
  in the scattering-matrix code (`calcmat.py`), which is **not** part of this
  repository

**Final step to the transport grid:** `examples/04_pot_to_rbf.py` (RBF
interpolation of δV(k), one model per z plane).

---

## What was deliberately left out

Present in the original `potcorr.py` / working directory, dropped here because
it is not part of descreening or rescreening:

| Dropped | Why |
|---|---|
| `get_greenfunc`, and the `kernel='gf'` branches | real-space Green-function kernel, unused in the paper |
| `get_vcoul_lr`, `get_vcoul_sr` | Ewald long/short-range splitting, unused |
| `get_pointc_pot_bare` | the point-charge baseline the paper argues against |
| `epsmat_init_interp`, `get_chimat_interp` | an abandoned interpolation variant |
| `read_chi_doped` | referenced `self.Chi0_doped_filename`, which `__init__` no longer sets |
| `read_eps.py:eps_model` | spline fit of the ε diagonal, used only for plots |
| `calcmat.py` | the scattering matrix elements, Eq. (S16)–(S20) — a separate step |
| `dat2xsf.py`, `get_eps.py`, `new_0.py`, `mpi_test.py`, `test.py`, `*_test.py` | scratch and superseded code |
| `output/`, `*.npy`, `myjob.*`, `kpt_test.dat` | run artefacts and data |
