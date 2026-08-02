# Notebooks

Preserved as the record of what was actually run for the paper. Outputs have
been stripped (18.8 MB → 0.3 MB); the code cells are otherwise untouched.

For the cleaned-up, importable versions of the same calculations see
`../examples/`.

| Notebook | What it does | Figures |
|---|---|---|
| `mos2_test.ipynb` | descreening driver (c1), 12×12 → 24×24 transfer (c2–c4), 3D-periodic rescreening (c6–c12) | — |
| `Na_test.ipynb` | the same descreening flow for a 3D system, `ncharge=9` (c33–c36) | — |
| `check_doped.ipynb` | 3D-periodic rescreening with the doped ε⁻¹, symmetry-accelerated | FIG. 2(d) |
| `2d_screening_chi.ipynb` | the quasi-2D open-boundary route, now `potcorr/quasi2d.py` | FIG. S2 |
| `rho_bare.ipynb` | descreening validation | FIG. 2(a–c) |
| `rho_bare2.ipynb` | vacuum-thickness convergence, z = 20–36 Å | FIG. S1 |
| `read_wan.ipynb` | Adler–Wiser χ⁰_c, fits and pickles `chi_{T}_{n}_rbf.pkl` | — |
| `wfc_read.ipynb` | CBM envelope ξ(z) → `xi_z.npy` / `xi_Gz.npy` | — |

## Running them

They import `lab`, `read_eps`, `const`, `io_xsf` as top-level modules, which is
the pre-packaging layout. Add a first cell:

```python
import _legacy_shim  # noqa: F401
```

and the old names resolve to the packaged modules. One rename is deliberately
not shimmed because it is a real API change:

```
ttkw.get_epsmat(...)  ->  ttkw.get_epsmat_from_chi(...)
```

`mos2_test.ipynb` c1 and c8 still use the old name and will raise
`AttributeError` as written.

## Data paths

Every notebook hard-codes absolute paths under
`/anvil/projects/x-che190065/rjguo/`. None of that data is in this repository.
