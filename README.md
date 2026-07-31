# Accurate *Ab Initio* Method for Charged Defect Scattering

Project page and method reference for:

> R. Guo, K. Kim, Z. Xiao, and Y. Liu,
> *Accurate Ab Initio Method for Charged Defect Scattering*,
> **Phys. Rev. Lett. 135, 126302 (2025)** — [doi:10.1103/42pc-t2bx](https://doi.org/10.1103/42pc-t2bx)

📄 **Project page:** https://yuanyue-liu-group.github.io/charged-defect-scattering/

## What the method does

Charged defect scattering matrix elements require the perturbation potential ΔV, whose Hartree
component is long ranged — too long ranged to converge in any affordable DFT supercell. Existing
approaches sidestep this by *assuming* a distribution for the bare (unscreened) excess charge:
a point charge, a Gaussian, or a Kohn–Sham defect orbital. Those assumptions limit accuracy.

This work determines the bare charge from first principles instead. Two operations do the work:

**Descreening** — recover the bare charge from the total density change of a *small* supercell,
using the bound-electron polarizability of the primitive cell:

```
Δρ_b = Δρ_tot − χ⁰_b · v · (Δρ_tot − Q/Ω)
```

Δρ_b turns out to be localized around the defect, which is why a small supercell suffices.

**Rescreening** — transfer Δρ_b to an arbitrarily large supercell and screen it there:

```
ΔV_H^(NS→CS) = W · Δρ_b,    W = (I − vχ⁰)⁻¹ v,    χ⁰ = χ⁰_b + χ⁰_c
```

Free carriers enter through χ⁰_c, evaluated with a modified Adler–Wiser formula rather than a
Thomas–Fermi or Lindhard model. For 2D materials the out-of-plane direction is kept in real space
so the open boundary and finite thickness are handled correctly.

## Main result

Applied to 2D MoS₂ with −1 charged sulfur vacancies (V<sub>S</sub><sup>−1</sup>) and neutral oxygen
substitutes of sulfur (O<sub>S</sub>), the charged-defect-limited mobility shows a strong
temperature dependence at low temperature and low carrier concentration — contrary to the
conventional expectation that defect-limited mobility is temperature insensitive. The neutral
defect shows no such behavior. The origin is free-carrier screening of the long-range scattering
potential.

## Pipeline

| Stage | Step | Tool |
|---|---|---|
| 1 | Ground-state DFT: primitive cell + PS / NS / CS supercells | Quantum ESPRESSO |
| 2 | Bound-electron polarizability χ⁰_b | BerkeleyGW |
| 3 | Descreening → Δρ_b | `potcorr` |
| 4 | Transfer Δρ_b to the large supercell | `potcorr` |
| 5 | Free-carrier polarizability χ⁰_c, per (T, n_e) | `potcorr` |
| 6 | Rescreening → ΔV_H | `potcorr` |
| 7 | Interpolate onto the transport mesh | `potcorr` |
| 8 | Scattering matrix elements | `potcorr` |
| 9 | Relaxation times + Boltzmann transport | [EDI](https://github.com/yuanyue-liu-group/EDI) |

Stage 5 is an independent branch and must be repeated for every temperature and carrier
concentration.

## Related

- [EDI](https://github.com/yuanyue-liu-group/EDI) — electron–defect interaction and transport solver

## Citation

```bibtex
@article{Guo2025ChargedDefect,
  title   = {Accurate Ab Initio Method for Charged Defect Scattering},
  author  = {Guo, Rongjing and Kim, Kwangrae and Xiao, Zhongcan and Liu, Yuanyue},
  journal = {Physical Review Letters},
  volume  = {135},
  number  = {12},
  pages   = {126302},
  year    = {2025},
  doi     = {10.1103/42pc-t2bx}
}
```

## Acknowledgments

Office of Naval Research (N00014-23-1-2400), NSF (2425545), Welch Foundation (F-1959), and
DOE EERE (DE-EE0052781). Calculations performed on the clusters of ACCESS, TACC, and NREL.
