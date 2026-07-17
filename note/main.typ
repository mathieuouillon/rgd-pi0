#import "template/lib.typ": *

#show: make-glossary
#register-glossary(acronyms)

#show: note.with(
  title: "Neutral Pion Production off Nuclear Targets",
  subtitle: "Multiplicity ratios, transverse-momentum broadening and beam-spin asymmetries in RG-D",
  author: "Mathieu Ouillon",
  affiliation: "Mississippi State University",
  email: "ouillon@jlab.org",
  doc-name: "RG-D π⁰ Analysis Note",
  version: "Draft — internal circulation",
)

#include "sections/01_introduction.typ"
#include "sections/02_dataset.typ"
#include "sections/03_event_selection.typ"
#include "sections/04_pi0_reconstruction.typ"
#include "sections/05_binning.typ"
#include "sections/06_background.typ"
#include "sections/07_multiplicity_ratio.typ"
#include "sections/08_pt_broadening.typ"
#include "sections/09_bsa.typ"
#include "sections/10_systematics.typ"
#include "sections/11_results.typ"
#include "sections/12_summary.typ"
#include "sections/99_appendix.typ"
