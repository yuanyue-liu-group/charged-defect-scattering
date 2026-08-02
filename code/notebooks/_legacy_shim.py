"""Compatibility shim for the notebooks in this directory.

The notebooks were written against the old flat layout, where `lab`,
`read_eps`, `const`, `io_xsf` and `utli` were top-level modules sitting next
to `potcorr.py`.  Run this first and their imports resolve against the
packaged versions:

    >>> import _legacy_shim   # noqa: F401
    >>> import potcorr, lab, read_eps, const, io_xsf

One rename is *not* shimmed, because it is a real API change:

    ttkw.get_epsmat(...)   ->   ttkw.get_epsmat_from_chi(...)

`mos2_test.ipynb` c1 and c8 still use the old name and will raise
AttributeError as written.
"""

import sys

import potcorr
from potcorr import bgw, cells, const, grid, io_xsf

for _name, _mod in [
    ("lab", cells),
    ("read_eps", bgw),
    ("const", const),
    ("io_xsf", io_xsf),
    ("utli", grid),
]:
    sys.modules.setdefault(_name, _mod)
