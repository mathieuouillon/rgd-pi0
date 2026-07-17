// Acronyms used throughout the note. First use auto-expands to long form,
// subsequent uses show the short form only. Reference an entry with `@KEY`
// (lowercase or uppercase as defined). For sentence-initial caps use
// `#Gls("KEY")`; force long form with `@KEY:long`.

#import "@preview/glossarium:0.5.10": make-glossary, register-glossary, print-glossary, gls, glspl

#let acronyms = (
  (key: "HTCC",     short: "HTCC",     long: "High Threshold Cherenkov Counter"),
  (key: "ECAL",     short: "ECAL",     long: "Electromagnetic Calorimeter"),
  (key: "PCAL",     short: "PCAL",     long: "Preshower Calorimeter"),
  (key: "ECIN",     short: "ECIN",     long: "Inner Electromagnetic Calorimeter"),
  (key: "ECOUT",    short: "ECOUT",    long: "Outer Electromagnetic Calorimeter"),
  (key: "FTOF",     short: "FTOF",     long: "Forward Time-of-Flight"),
  (key: "DC",       short: "DC",       long: "Drift Chamber"),
  (key: "MIP",      short: "MIP",      long: "minimum ionizing particle"),
  (key: "FD",       short: "FD",       long: "Forward Detector"),
  (key: "CD",       short: "CD",       long: "Central Detector"),
  (key: "FT",       short: "FT",       long: "Forward Tagger"),
  (key: "EB",       short: "EB",       long: "Event Builder"),
  (key: "CLAS12",   short: "CLAS12",   long: "CEBAF Large Acceptance Spectrometer at 12 GeV"),
  (key: "CEBAF",    short: "CEBAF",    long: "Continuous Electron Beam Accelerator Facility"),
  (key: "COATJAVA", short: "COATJAVA", long: "CLAS12 Offline Analysis Toolkit (Java)"),
  (key: "RGA",      short: "RG-A",     long: "Run Group A"),
  (key: "RGC",      short: "RG-C",     long: "Run Group C"),
  (key: "RGD",      short: "RG-D",     long: "Run Group D"),
  (key: "SF",       short: "SF",       long: "sampling fraction"),
  (key: "DSCB",     short: "DSCB",     long: "double-sided Crystal Ball"),
  (key: "CCDB",     short: "CCDB",     long: "Calibration Constants Database"),
  (key: "RCDB",     short: "RCDB",     long: "Run Conditions Database"),
  (key: "TMD",      short: "TMD",      long: "Transverse Momentum Distribution"),
  // ---- Added for the pi0 analysis note ----
  (key: "DIS",      short: "DIS",      long: "deep inelastic scattering"),
  (key: "SIDIS",    short: "SIDIS",    long: "semi-inclusive deep inelastic scattering"),
  (key: "QCD",      short: "QCD",      long: "Quantum Chromodynamics"),
  (key: "BSA",      short: "BSA",      long: "beam-spin asymmetry"),
  (key: "SSA",      short: "SSA",      long: "single-spin asymmetry"),
  (key: "MLM",      short: "MLM",      long: "maximum likelihood method"),
  (key: "FIFO",     short: "FIFO",     long: "first-in, first-out"),
  (key: "HWP",      short: "HWP",      long: "half-wave plate"),
  (key: "MC",       short: "MC",       long: "Monte Carlo"),
  (key: "GBT",      short: "GBT",      long: "gradient boosted decision tree"),
  (key: "SWIF2",    short: "SWIF2",    long: "Scientific Workflow Indefatigable Factotum"),
  (key: "HIPO",     short: "HIPO",     long: "High Performance Output"),
  (key: "FF",       short: "FF",       long: "fragmentation function"),
  (key: "nFF",      short: "nFF",      long: "nuclear-modified fragmentation function"),
  (key: "PDF",      short: "PDF",      long: "parton distribution function"),
)

#let glossary-section() = [
  = Glossary
  #print-glossary(acronyms, show-all: true, disable-back-references: true)
]
