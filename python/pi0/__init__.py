"""rgd-pi0: the Python extraction stage of the RG-D pi0 nuclear multiplicity analysis.

Reads Stage B's binned ROOT file (uproot only -- no ROOT, no dictionary) and
produces the four observables:

===============  ==========================================================
:mod:`pi0.ratio`       ``R_A`` per 4D bin
:mod:`pi0.broadening`  ``Delta<pT2>`` per 3D bin, **sideband-subtracted**
:mod:`pi0.bsa`         ``A_LU``, dilution-corrected, polarization mandatory
:mod:`pi0.qa`          diagnostics, including the abscissa-vs-centre plot
===============  ==========================================================

with :mod:`pi0.config` owning every cut value and every index formula,
:mod:`pi0.io` owning the file reading and the provenance gate, and
:mod:`pi0.extract` owning the yield extraction and the count-weighted abscissa.

Three invariants hold across the package:

* **No cut value is hard-coded.** Everything comes from ``config/cuts.json``
  and the two grid JSONs, and a missing key raises rather than defaults.
* **No result is ever reported at a geometric bin centre.** Every abscissa is
  the count-weighted, sideband-subtracted mean; a bin whose abscissa cannot be
  computed is dropped and counted.
* **Provenance is enforced, not decorated.** A file whose provenance says its
  photons came from a fallback GBT model, or whose grids are placeholders, will
  not silently produce a physics number.
"""

from __future__ import annotations

__version__ = "0.1.0"

__all__ = ["config", "io", "extract", "ratio", "broadening", "bsa", "qa", "__version__"]

# Submodules are NOT imported eagerly: `python -m pi0.qa` would then find pi0.qa
# already in sys.modules before runpy executes it, which runpy warns about and
# which can execute a module twice. `from pi0 import qa` imports it on demand.
