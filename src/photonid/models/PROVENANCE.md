# Photon GBT models — provenance

**Do not edit these files.** They are generated CatBoost model data, vendored
verbatim. An edit — even a reformat — silently changes the classifier's output,
and nothing in the test suite would catch a one-digit change among 203,339
numeric tokens.

## What they are

Five gradient-boosted-tree classifiers that separate real photons from the
neutral background (principally neutrons and split-off calorimeter clusters) in
the CLAS12 electromagnetic calorimeter. Each takes a 45-component feature vector
(`FloatFeatureCount = 45`) built by `src/photonid/Features.cpp`, and returns a
score `x`; the photon is accepted when `sigmoid(x) > threshold`
(`photon.gbt_threshold` in `config/cuts.json`, currently 0.78).

## Where they come from

They are CLAS12 training products, distributed by **Iguana** (LGPL v3):

    https://github.com/JeffersonLab/iguana
    src/iguana/algorithms/clas12/PhotonGBTFilter/models/*.cpp

Verified byte-for-byte identical in payload to Iguana's copies — 203,339
numeric tokens each, compared token-by-token. Only the extension differs
(`.cpp` → `.hpp`), because here they are `#include`d rather than compiled as
separate translation units.

See `NOTICE` at the repository root for the attribution, and `COPYING.LESSER` +
`COPYING` for the licence.

## Run coverage — read this before trusting a result

The model map (`src/photonid/RunRangeModelMap.cpp`) is a faithful port of
Iguana's, and covers:

| runs          | model                  |
|---------------|------------------------|
| 5032–5666     | RGA in/outbending pass 1, 2 |
| 6156–6783     | RGA inbending pass 1, 2     |
| 11093–11571   | RGA in/outbending pass 1, 2 |
| 16042–16772   | RGC Summer 2022 pass 1      |

**It stops at run 16772. There is no trained model for RG-D (runs
18305–19131), in Iguana or anywhere else.**

Iguana's response is to log a warning and fall back to `RGA_inbending_pass1`.
The superseded clas-framework port dropped even the warning and fell back
silently — which is where this analysis's largest unquantified photon
systematic came from.

This project instead **refuses**: `select_model()` throws unless
`photon.allow_rga_fallback` is explicitly set to `true` in `config/cuts.json`,
and when it is, the fallback is stamped into the output provenance
(`gbt.fallback_used`). Applying an RG-A *inbending* model to RG-D *outbending*
nuclear-target data is a real assumption, not a formality: the feature vector
includes neighbour multiplicity, and calorimeter occupancy differs between
LD₂ and Sn — which is precisely the direction in which the efficiency would
fail to cancel in a target ratio.

Getting a real RG-D model is the single highest-value fix available to the
photon selection.
