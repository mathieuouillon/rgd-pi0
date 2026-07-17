# Unit tests for src/core -- pure physics, no ROOT and no HIPO, no data file.
#
# ASSUMPTION: the parent exposes `pi0_core_dep`, a declare_dependency() that
# carries the pi0_core library *and* `src/` as an include directory (so that
# `#include "core/Kinematics.hpp"` resolves). If the parent instead exposes the
# raw library plus the include dir, replace the `dependencies:` line below with
#
#     link_with: pi0_core,
#     include_directories: inc_src,
#
# This file is reached via `subdir('tests')` from the top-level meson.build,
# which must be added there (it is not, as of writing).
#
# Catch2 is not installed system-wide on the dev machine, so this falls back to
# the wrap at external/catch2.wrap (external/ is this project's subproject_dir).
# `catch2-with-main` supplies Catch2's own main(); test_core.cpp defines none.

catch2_dep = dependency(
  'catch2-with-main',
  version: '>=3.0.0',
  required: true,
  fallback: ['catch2', 'catch2_with_main_dep'],
  # Catch2's own SelfTest suite is built by default and does not compile
  # cleanly under Apple clang 21 (TEMPLATE_PRODUCT_TEST_CASE / get_wrapper).
  # We only want the library, and we do not want to install it.
  default_options: ['tests=false', 'install=false'],
)

test_core_exe = executable(
  'test_core',
  sources: files('test_core.cpp'),
  dependencies: [pi0_core_dep, catch2_dep],
  # Tests are built, never shipped.
  install: false,
  # The core is meant to be unit-testable standalone; keep the test binary
  # honest about warnings even though the project default is already level 2.
  override_options: ['warning_level=3'],
)

test(
  'core',
  test_core_exe,
  # The Trento random-configuration scan runs 2000 configurations; everything
  # else is closed-form. This is generous for a pure-arithmetic suite.
  timeout: 60,
  # TAP lets meson count individual assertions rather than one pass/fail.
  # CAVEAT: Catch2's TAP reporter emits WARN(...) as a line meson's TAP parser
  # rejects ("Unknown TAP output lines"), which reports a fully passing run as
  # an ERROR. Use INFO(...) + CHECK(...) in test_core.cpp instead of WARN.
  protocol: 'tap',
  args: ['--reporter', 'tap'],
  suite: 'unit',
)

# ---------------------------------------------------------------------------
# test_vertex -- the vz correction (src/vertex).
#
# Unlike test_core this one DOES need a data file: the correction is defined by
# data/Vz/vz_corrector_params.txt, so the loader and the reference-value tests
# read it. The parent exposes its absolute path as `vz_params_file` and has
# already checked that it exists, so a missing file is a configure-time error,
# not a mid-test throw.
#
# The path is injected as a macro rather than passed as an argv argument because
# Catch2 supplies main() and owns the command line; a stray positional argument
# would be parsed as a test-name filter. test_vertex.cpp #errors out if the
# macro is not defined -- it has no fallback path, so it cannot silently test
# the wrong file (or silently skip).
#
# Depends on pi0_vertex_dep, NOT pi0_core_dep: src/vertex includes nothing from
# src/core and needs no ROOT. Keeping the link line honest keeps the "this is
# plain C++17" claim checkable.
# ---------------------------------------------------------------------------
test_vertex_exe = executable(
  'test_vertex',
  sources: files('test_vertex.cpp'),
  dependencies: [pi0_vertex_dep, catch2_dep],
  cpp_args: ['-DPI0_VZ_PARAMS_FILE="' + vz_params_file + '"'],
  install: false,
  override_options: ['warning_level=3'],
)

test(
  'vertex',
  test_vertex_exe,
  # Reads a 460K params file a handful of times; everything else is arithmetic.
  timeout: 60,
  protocol: 'tap',
  args: ['--reporter', 'tap'],
  suite: 'unit',
)

# ---------------------------------------------------------------------------
# test_photonid -- the ported CLAS12 GBT photon classifier (src/photonid).
#
# No data file and no HIPO: the REC::Particle / REC::Calorimeter RVec columns
# are hand-built in the test. That is the point of build_features() taking RVecs
# rather than a hipo::bank -- the port is checkable without a skim.
#
# warning_level=3 despite pi0_photonid itself building at level 0: the level-0
# override there is scoped to the target that includes the verbatim CatBoost
# codegen in photonid/models/. Nothing generated is included from here (the
# models are reached through pi0::photonid::models::*(), which hands back a
# function pointer), so this binary is held to the strict level and the
# hand-written headers get checked at it.
# ---------------------------------------------------------------------------
test_photonid_exe = executable(
  'test_photonid',
  sources: files('test_photonid.cpp'),
  dependencies: [pi0_photonid_dep, catch2_dep],
  install: false,
  override_options: ['warning_level=3'],
)

test(
  'photonid',
  test_photonid_exe,
  # Pure arithmetic over hand-built events; no model is ever evaluated (the
  # scoring tests use stub margins), so nothing here is slow.
  timeout: 60,
  protocol: 'tap',
  args: ['--reporter', 'tap'],
  suite: 'unit',
)

# ---------------------------------------------------------------------------
# test_selection -- the cut config (src/config) and electron/photon selection
# (src/selection).
#
# Reads the REAL config/cuts.json, whose absolute path the parent exposes as
# `cuts_json_file` and has already checked exists. Testing against a hand-written
# copy of the config would prove only that the copy parses; the whole point of
# these tests is that the SHIPPED file says what the analysis claims it says.
#
# Same macro-not-argv reasoning as test_vertex above: Catch2 owns the command
# line, so a positional path would be read as a test-name filter.
#
# PI0_TEST_TMPDIR is this test's own build directory. The "throws on a missing
# key" tests work by deleting one key from a parsed copy of the real config and
# writing the result somewhere; they need a writable directory that is not the
# source tree. Each test removes its file.
#
# nlohmann_json_dep is needed HERE, directly: pi0_config_dep deliberately does
# not re-export it (it is private to Cuts.cpp), but this test mutates a parsed
# config document to build the negative cases, so it is a real user of the
# library in its own right.
# ---------------------------------------------------------------------------
test_selection_exe = executable(
  'test_selection',
  sources: files('test_selection.cpp'),
  dependencies: [pi0_selection_dep, nlohmann_json_dep, catch2_dep],
  cpp_args: [
    '-DPI0_CUTS_JSON="' + cuts_json_file + '"',
    '-DPI0_TEST_TMPDIR="' + meson.current_build_dir() + '"',
  ],
  install: false,
  override_options: ['warning_level=3'],
)

test(
  'selection',
  test_selection_exe,
  # Parses a ~15K JSON file once per negative case (about 20 of them) and does
  # closed-form arithmetic otherwise.
  timeout: 60,
  protocol: 'tap',
  args: ['--reporter', 'tap'],
  suite: 'unit',
)

# ---------------------------------------------------------------------------
# test_binning -- the frozen factorized binning (src/core/Binning).
#
# Reads the REAL config/binning/*.json, whose absolute paths the parent exposes
# as `grid_a_json_file` / `grid_b_json_file` and has already checked exist. Same
# reasoning as test_selection's use of cuts_json_file: testing against a
# hand-written copy of the grid would prove only that the copy parses, and the
# whole point is that the SHIPPED grids say what the analysis claims they say.
#
# Same macro-not-argv reasoning too: Catch2 supplies main() and owns the command
# line, so a positional path would be parsed as a test-name filter.
#
# PI0_TEST_TMPDIR is this test's own build directory, shared with test_selection
# and used for the same job: the negative cases work by mutating a parsed copy
# of a real grid (deleting a key, reversing an axis) and writing the result
# somewhere outside the source tree. Each test removes its file.
#
# nlohmann_json_dep is needed HERE, directly, and for the same reason
# test_selection needs it: pi0_core_dep deliberately does not re-export it (it
# is private to Binning.cpp), but this test builds its negative cases by
# mutating a parsed grid document, so it is a real user of the library in its
# own right.
# ---------------------------------------------------------------------------
test_binning_exe = executable(
  'test_binning',
  sources: files('test_binning.cpp'),
  dependencies: [pi0_core_dep, nlohmann_json_dep, catch2_dep],
  cpp_args: [
    '-DPI0_GRID_A_JSON="' + grid_a_json_file + '"',
    '-DPI0_GRID_B_JSON="' + grid_b_json_file + '"',
    '-DPI0_TEST_TMPDIR="' + meson.current_build_dir() + '"',
  ],
  install: false,
  override_options: ['warning_level=3'],
)

test(
  'binning',
  test_binning_exe,
  # The round-trip tests walk all 1400 4D bins and all 280 3D bins, which is
  # pure integer arithmetic; the negative cases parse a ~5K JSON file about 20
  # times.
  timeout: 60,
  protocol: 'tap',
  args: ['--reporter', 'tap'],
  suite: 'unit',
)

# ---------------------------------------------------------------------------
# test_donorpool -- the frozen event-mixing donor pool (src/stageB_bin).
#
# Reads the REAL config/cuts.json, for the same reason test_selection and
# test_binning read their real inputs: the note's eq:pool-key claims 224 pool
# bins, and that is a claim about the SHIPPED file, not about a copy of it that a
# test happens to agree with. Same macro-not-argv reasoning as those two, too --
# Catch2 owns the command line, so a positional path would be read as a
# test-name filter.
#
# Depends on pi0_config_dep rather than pi0_stageB_bin_dep alone: the pool grid
# arrives through Cuts::mixing, so this test is a user of the config as well as
# of the pool. pi0_stageB_bin_dep comes along inside pi0_config_dep (Cuts.hpp
# includes stageB_bin/PoolGrid.hpp for MixingCuts), so it is not listed twice.
#
# No PI0_TEST_TMPDIR: unlike test_selection and test_binning this suite writes no
# files. Its negative cases mutate a PoolGridCuts in memory, because a donor pool
# is built from a struct rather than parsed from a document -- which is the point
# of MixingCuts living in src/stageB_bin and not behind a loader.
# ---------------------------------------------------------------------------
test_donorpool_exe = executable(
  'test_donorpool',
  sources: files('test_donorpool.cpp'),
  dependencies: [pi0_config_dep, catch2_dep],
  cpp_args: ['-DPI0_CUTS_JSON="' + cuts_json_file + '"'],
  install: false,
  override_options: ['warning_level=3'],
)

test(
  'donorpool',
  test_donorpool_exe,
  # The uniformity test builds 400 pools of 100 offers each -- 40k draws through
  # mt19937_64, which is microseconds. Everything else is a handful of offers.
  timeout: 60,
  protocol: 'tap',
  args: ['--reporter', 'tap'],
  suite: 'unit',
)

# ---------------------------------------------------------------------------
# test_inputchain -- the N slim files of one Stage B run
# (src/stageB_bin/InputChain.hpp): their canonical order, the cross-checks that
# must hold across them, and the seed they derive.
#
# No config file and no PI0_TEST_TMPDIR, and unlike every other suite here no
# real input of any kind: a Provenance is an ordered vector of string pairs, so
# a chain -- including the corrupt ones that must be refused -- is built in
# code. That is the point of InputChain.hpp taking a ChainInput rather than a
# TFile. The chains these tests must reject (an LD2 file chained with an Sn one,
# a truncated skim hidden among full ones) are exactly the chains it would be
# most expensive to keep real fixtures of.
#
# Depends on pi0_core_dep alone, which is the whole dependency: the header
# includes util/Provenance.hpp and util/Sha256.hpp -- both header-only, both
# reached through the `src/` include directory that every pi0_*_dep carries --
# and ROOT, for TFile, which Provenance.hpp needs and pi0_core_dep supplies. It
# does NOT need pi0_config_dep or pi0_stageB_bin_dep: the header parses no cut
# and constructs no DonorPool, and the absence on this line is what keeps
# "deciding which files are one run is separable from reading them" checkable.
# ---------------------------------------------------------------------------
test_inputchain_exe = executable(
  'test_inputchain',
  sources: files('test_inputchain.cpp'),
  dependencies: [pi0_core_dep, catch2_dep],
  install: false,
  override_options: ['warning_level=3'],
)

test(
  'inputchain',
  test_inputchain_exe,
  # Sorts and refuses chains of at most three hand-built inputs, and hashes a
  # handful of 64-character strings. There is nothing here that is not immediate.
  timeout: 60,
  protocol: 'tap',
  args: ['--reporter', 'tap'],
  suite: 'unit',
)
