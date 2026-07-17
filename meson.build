project(
  'rgd_pi0_analysis',
  'cpp',
  version: '0.1.0',
  meson_version: '>=1.2',
  # hipo-cpp is vendored at external/hipo-cpp, so point meson's subproject
  # machinery at external/ instead of the default subprojects/. This is also
  # where catch2.wrap lives.
  subproject_dir: 'external',
  default_options: {
    'cpp_std':       'c++17',  # must match ROOT: `root-config --cflags` -> -std=c++17
    'buildtype':     'release',
    'warning_level': '2',
  },
)

fs = import('fs')

# ---------------------------------------------------------------------------
# ROOT
#
# Found via CMake config mode, exactly as hipo-cpp finds it (see
# external/hipo-cpp/meson.build). Using the same method and module list keeps
# us and the subproject on one consistent ROOT.
#
# Note we deliberately do NOT use `root-config --cflags`: ROOT's CMake config
# pulls in the VDT::VDT imported target via find_dependency(Vdt), so the -I for
# <vdt/vdtMath.h> (included by ROOT 6.40's <ROOT/RVec.hxx>) comes along for
# free. root-config does not emit that -I. See the ROOT_INCLUDE_PATH note below
# for the rootcling side of the same problem.
# ---------------------------------------------------------------------------
ROOT_dep = dependency(
  'ROOT',
  method: 'cmake',
  version: '>=6.36',
  modules: [
    'ROOT::Core',
    'ROOT::Hist',
    'ROOT::Physics',
    'ROOT::RIO',
    'ROOT::ROOTDataFrame',
    'ROOT::ROOTVecOps',
    'ROOT::Tree',
  ],
)

# ---------------------------------------------------------------------------
# rootcling / ROOT_INCLUDE_PATH preflight
#
# THE TRAP: the C++ compile is fine (ROOT_dep carries VDT::VDT's include dir),
# but rootcling is not the C++ compiler and does not see meson's -I flags.
# hipo-cpp's dictionary custom_target hardcodes its rootcling flags and never
# forwards ROOT_dep's compile args, so generating libHipoDataFrame's dictionary
# dies with:
#
#     fatal error: 'vdt/vdtMath.h' file not found
#
# The fix is the environment variable ROOT_INCLUDE_PATH, which is the mechanism
# rootcling intends for exactly this (upstream hipo-cpp uses it too, see
# external/hipo-cpp/meson/this_hipo.sh.in). rootcling reads it at BUILD time, so
# it must be exported in the shell that runs ninja -- meson machine files have
# no [env] section, and meson.add_devenv() only affects `meson devenv`, not the
# build. Hence: check here and fail loudly at configure time rather than let the
# user hit a confusing rootcling error ten minutes into the build.
# ---------------------------------------------------------------------------
# WHERE vdt LIVES IS PLATFORM-DEPENDENT, so ask rather than guess. On the JLab
# farm it ships inside ROOT's own include directory; under Homebrew it is a
# separate formula in /opt/homebrew/include. ROOT's incdir is tried FIRST, and it
# is the case that needs no environment at all -- see below.
root_incdir = run_command('root-config', '--incdir', check: true).stdout().strip()

vdt_incdir = ''
foreach candidate : [root_incdir, '/opt/homebrew/include', '/usr/local/include', '/usr/include']
  if vdt_incdir == '' and fs.exists(candidate / 'vdt' / 'vdtMath.h')
    vdt_incdir = candidate
  endif
endforeach

if vdt_incdir == ''
  error(
    'vdt/vdtMath.h not found, but ROOT\'s <ROOT/RVec.hxx> includes it.\n' +
    '\n' +
    'Looked in ROOT\'s own include dir (' + root_incdir + '), then\n' +
    '/opt/homebrew/include, /usr/local/include, /usr/include.\n' +
    '\n' +
    'On macOS:  brew install vdt\n' +
    'Elsewhere: install vdt, or use a ROOT build that ships it.',
  )
endif

# ROOT_INCLUDE_PATH is only needed when vdt sits somewhere rootcling does not
# already look. rootcling always searches ROOT's own incdir -- that is where
# ROOT's headers are, so it reads it by construction -- and demanding the
# variable in that case would be demanding the user point it at a directory it
# cannot miss. That is the JLab farm, where this builds with no environment set.
if vdt_incdir != root_incdir
  root_include_path = run_command('printenv', 'ROOT_INCLUDE_PATH', check: false).stdout().strip()

  if not root_include_path.split(':').contains(vdt_incdir)
    error(
      'ROOT_INCLUDE_PATH does not contain "' + vdt_incdir + '".\n' +
      '\n' +
      'vdt is at ' + vdt_incdir + ', which is NOT inside ROOT\'s include dir\n' +
      '(' + root_incdir + '), so rootcling will not find it by itself. It needs\n' +
      '<vdt/vdtMath.h> to generate the libHipoDataFrame dictionary; without it\n' +
      'the build fails with\n' +
      '  fatal error: \'vdt/vdtMath.h\' file not found\n' +
      '\n' +
      'Export it in the shell you build from, then re-run meson setup:\n' +
      '\n' +
      '    export ROOT_INCLUDE_PATH=' + vdt_incdir + '\n' +
      '\n' +
      'It is needed for `ninja` too, not just `meson setup`. See README.md.',
    )
  endif
endif

# Make ROOT_INCLUDE_PATH available to `meson devenv -C build`, so that running
# the built programs (and any interactive ROOT session that loads our
# dictionaries) resolves headers the same way the build did.
devenv = environment()
devenv.prepend('ROOT_INCLUDE_PATH', vdt_incdir)
meson.add_devenv(devenv)

# ---------------------------------------------------------------------------
# hipo-cpp subproject -> libhipo4 + libHipoDataFrame
#
# DO NOT MODIFY external/hipo-cpp. It is a pristine v4.4.1 checkout; everything
# we need is pulled out of it read-only via get_variable() below.
# ---------------------------------------------------------------------------
hipo_proj = subproject(
  'hipo-cpp',
  default_options: {
    'dataframes':     true,   # build the RDataFrame extension -> libHipoDataFrame
    'build_examples': false,
    'build_tests':    false,  # hipo-cpp's own option, not ours (see meson.options)
  },
)

# hipo4/meson.build ends with a declare_dependency() named `hipo_dep`; it
# already carries include_directories('..') (so we #include "hipo4/reader.h"),
# link_with: hipo_lib, and fmt_dep. Note -D__LZ4__ is deliberately NOT needed by
# consumers: it only guards code in hipo4/*.cpp, never in the public headers.
hipo_dep = hipo_proj.get_variable('hipo_dep')

# The dataframes extension (extensions/dataframes/meson.build) defines
#   dataframe_lib     : library() object -- libHipoDataFrame (+ .rootmap, _rdict.pcm)
#   dataframe_headers : list of ABSOLUTE header paths (RHipoDS.hxx, RIguanaDS.hxx,
#                       ThreadedAlgo.hxx) -- absolute because rootcling needs that
#                       for its generated #include directives
# but -- unlike hipo4/ -- it never calls declare_dependency(). So we assemble one
# here from the pieces it does expose.
#
# The headers are installed flat and are included as `#include "RHipoDS.hxx"`,
# so the -I is their containing directory. We pass it as a raw compile_arg
# rather than include_directories(): the latter would point from the master
# project into a subproject's source tree and trip meson's subproject sandbox.
dataframe_incdir = fs.parent(hipo_proj.get_variable('dataframe_headers')[0])

HipoDataFrame_dep = declare_dependency(
  link_with:    hipo_proj.get_variable('dataframe_lib'),
  compile_args: ['-I' + dataframe_incdir],
  dependencies: [hipo_dep, ROOT_dep],
)

# ---------------------------------------------------------------------------
# nlohmann/json -- the JSON parser, used by exactly two loaders.
#
# Declared HERE, above pi0_core, rather than down beside pi0_config where it
# used to live: pi0_core now carries src/core/Binning.cpp, whose Binning::load()
# reads config/binning/*.json, and meson resolves these in order.
#
# It is a PRIVATE dependency of both libraries that use it -- neither Cuts.hpp
# nor Binning.hpp includes it, so it appears in the static_library() calls and
# in NEITHER declare_dependency(). Nothing downstream inherits it. The standing
# rule is unchanged and unweakened: Cuts::load() is the only reader of
# config/cuts.json, and Binning::load() is the only reader of
# config/binning/*.json. Neither file may be parsed anywhere else.
#
# The system pkg-config package is used when present; external/nlohmann_json.wrap
# is the fallback so a fresh checkout configures without brew. Header-only, so
# the wrap costs a download and essentially no build time.
# ---------------------------------------------------------------------------
nlohmann_json_dep = dependency(
  'nlohmann_json',
  version:  '>=3.9.0',  # json::type_name() + the parse_error hierarchy
  fallback: ['nlohmann_json', 'nlohmann_json_dep'],
)

# ---------------------------------------------------------------------------
# pi0_core -- pure analysis code, no HIPO, no ROOT beyond ROOT::VecOps::RVec.
# Everything here must be unit-testable without a data file.
#
# Meson has no source globbing by design, so src/core/*.cpp is spelled out.
# ADD NEW src/core/*.cpp FILES TO THIS LIST or they will not be compiled.
#
# Constants.hpp and Types.hpp are intentionally absent: both are header-only
# (inline constexpr / inline member functions), so they have no .cpp to compile.
#
# Binning.cpp is the one file here that can read a file, via Binning::load().
# That does not break the "unit-testable without a data file" contract above:
# the geometry -- Grid1D::find(), the flat index formulae, provenance_hash() --
# is pure and is tested on hand-built grids, and load() is tested on documents
# the test writes itself. The grids are NOT data in the sense vz_corrector_params
# is; they are config, they live in the repo, and they are small.
# ---------------------------------------------------------------------------
pi0_core_sources = files(
  'src/core/Binning.cpp',
  'src/core/Kinematics.cpp',
  'src/core/Pairing.cpp',
)

# Project-relative includes: #include "core/Kinematics.hpp"
pi0_inc = include_directories('src')

pi0_core = static_library(
  'pi0_core',
  pi0_core_sources,
  include_directories: pi0_inc,
  # ROOT for ROOT::VecOps::RVec only; nlohmann is PRIVATE to Binning.cpp and is
  # deliberately absent from pi0_core_dep below.
  dependencies:        [ROOT_dep, nlohmann_json_dep],
  install:             false,
)

pi0_core_dep = declare_dependency(
  link_with:            pi0_core,
  include_directories:  pi0_inc,
  dependencies:         [ROOT_dep],
)

# ---------------------------------------------------------------------------
# pi0_vertex -- the RG-D vz correction (src/vertex).
#
# A library of its own rather than more files in pi0_core, for two reasons:
#   * pi0_core is documented above as "unit-testable without a data file", and
#     the vz correction is *defined* by data/Vz/vz_corrector_params.txt -- its
#     loader test needs that file. Keeping it separate preserves that contract.
#   * pi0_core pulls in ROOT for RVec; this needs nothing but the C++ stdlib.
#
# It takes pi0_inc only for the project-relative include style
# (#include "vertex/VzCorrector.hpp"); it includes nothing from core/.
#
# Declared BEFORE pi0_config because Cuts.hpp includes vertex/VzCorrector.hpp
# for VertexCuts (the /vertex block's runtime image, defined next to the code
# that consumes it -- same arrangement as PairingCuts in core/Pairing.hpp), and
# meson resolves these in order.
#
# ADD NEW src/vertex/*.cpp FILES TO THIS LIST or they will not be compiled.
# ---------------------------------------------------------------------------
pi0_vertex = static_library(
  'pi0_vertex',
  files('src/vertex/VzCorrector.cpp'),
  include_directories: pi0_inc,
  install:             false,
)

pi0_vertex_dep = declare_dependency(
  link_with:           pi0_vertex,
  include_directories: pi0_inc,
)

# Absolute path to the vz correction parameters. The file is data, not a build
# input, so it is not copied into the build dir; callers open it where it lives.
# Checked here so a missing/renamed file is a configure-time error rather than a
# runtime throw inside a test.
vz_params_file = meson.project_source_root() / 'data' / 'Vz' / 'vz_corrector_params.txt'

if not fs.exists(vz_params_file)
  error('vz corrector parameters not found at ' + vz_params_file)
endif

# ---------------------------------------------------------------------------
# pi0_stageB_bin -- Stage B: the frozen event-mixing donor pool (src/stageB_bin).
#
# Declared BEFORE pi0_config, like pi0_vertex above and for the same reason:
# Cuts.hpp includes stageB_bin/PoolGrid.hpp for MixingCuts (the /mixing block's
# runtime image, defined next to the code that consumes it -- the same
# arrangement as PairingCuts in core/Pairing.hpp and VertexCuts in
# vertex/VzCorrector.hpp), and meson resolves these in order.
#
# THE ARROW POINTS THIS WAY ROUND ON PURPOSE. It would seem more natural for a
# "Stage B" library to depend on the config rather than the config on it -- but
# DonorPool takes a PoolGridCuts, not a Cuts, so it never includes Cuts.hpp and
# there is no cycle. The payoff is that the donor pool is unit-testable with a
# grid built in code, with no JSON parsed and no config file on disk, which is
# the same contract src/core and src/vertex hold to.
#
# It takes pi0_core_dep for Grid1D (core/Binning.hpp): the pool grid's axes reuse
# the project's ONE axis type rather than introducing a second edge-vector-with-
# a-lookup. Sharing the type is emphatically NOT sharing the edges -- the pool
# grid's Q^2/x_B edges are hand-authored config product edges and Grid A's are
# machine-fitted equal-statistics edges. See src/stageB_bin/PoolGrid.hpp.
#
# ADD NEW src/stageB_bin/*.cpp FILES TO THIS LIST or they will not be compiled.
# ---------------------------------------------------------------------------
pi0_stageB_bin = static_library(
  'pi0_stageB_bin',
  files(
    'src/stageB_bin/DonorPool.cpp',
    'src/stageB_bin/PoolGrid.cpp',
  ),
  include_directories: pi0_inc,
  dependencies:        [pi0_core_dep],  # Grid1D, from core/Binning.hpp
  install:             false,
)

pi0_stageB_bin_dep = declare_dependency(
  link_with:           pi0_stageB_bin,
  include_directories: pi0_inc,
  # PUBLIC: PoolGrid.hpp includes core/Binning.hpp, so anyone who includes it
  # needs Grid1D's definition and its out-of-line nbins()/find().
  dependencies:        [pi0_core_dep],
)

# ---------------------------------------------------------------------------
# pi0_config -- the runtime image of config/cuts.json (the Cuts struct).
#
# Separate from pi0_selection because Cuts is what EVERY other module needs
# (selection, vertex, photonid, the skim), and none of them should have to link
# the electron/photon selection just to read a threshold.
#
# It depends on pi0_core for PairingCuts, on pi0_vertex for VertexCuts, and on
# pi0_stageB_bin for MixingCuts: all three of those blocks are DEFINED by the
# library that consumes them, and aggregated here. The arrows all point into
# pi0_config, never out of it, so there is no cycle -- pi0_vertex still depends
# on nothing but the C++ stdlib, and every one of those libraries stays usable
# (and testable) without ever parsing a JSON file.
#
# nlohmann/json is a PRIVATE dependency: Cuts.hpp does not include it, so it
# appears in the static_library() below and NOT in the declare_dependency().
# Nothing downstream inherits it, and nothing else in the project may use it to
# read a cut -- Cuts::load() is the only reader of cuts.json. It is declared
# near the top of this file, above pi0_core, because src/core/Binning.cpp needs
# it too (for config/binning/*.json, which is not a cut file).
# ---------------------------------------------------------------------------
pi0_config = static_library(
  'pi0_config',
  files('src/config/Cuts.cpp'),
  include_directories: pi0_inc,
  # core for PairingCuts, vertex for VertexCuts, stageB_bin for MixingCuts.
  dependencies:        [pi0_core_dep, pi0_vertex_dep, pi0_stageB_bin_dep, nlohmann_json_dep],
  install:             false,
)

pi0_config_dep = declare_dependency(
  link_with:           pi0_config,
  include_directories: pi0_inc,
  # All three are PUBLIC: Cuts.hpp includes core/Pairing.hpp,
  # vertex/VzCorrector.hpp and stageB_bin/PoolGrid.hpp, so anyone who includes
  # Cuts.hpp needs them too.
  dependencies:        [pi0_core_dep, pi0_vertex_dep, pi0_stageB_bin_dep],
)

# ---------------------------------------------------------------------------
# pi0_selection -- electron and photon selection (src/selection).
#
# Deliberately does NOT contain the vertex cut: that is target-dependent, needs
# the correction parameterisation, and for the solid foils it is what decides
# which nucleus an event belongs to. It lives in pi0_vertex, and pass_electron()
# takes its verdict as a bool -- so nothing here calls into pi0_vertex, and
# there is no cycle. (pi0_vertex does arrive on the link line anyway, via
# pi0_config_dep, because Cuts carries the /vertex block. That is the aggregate
# config header, not a dependency of the selection logic: the bool argument is
# still the only way a vertex verdict reaches pass_electron.)
#
# The GBT score is likewise absent: pi0_photonid owns the Catboost models, and
# pass_photon() takes the score (or a callable producing it). Same reason, same
# non-cycle.
#
# ADD NEW src/selection/*.cpp FILES TO THIS LIST or they will not be compiled.
# ---------------------------------------------------------------------------
pi0_selection = static_library(
  'pi0_selection',
  files(
    'src/selection/ElectronSelection.cpp',
    'src/selection/PhotonSelection.cpp',
    'src/selection/SamplingFraction.cpp',
  ),
  include_directories: pi0_inc,
  dependencies:        [pi0_config_dep, ROOT_dep],  # ROOT for RVec only
  install:             false,
)

pi0_selection_dep = declare_dependency(
  link_with:           pi0_selection,
  include_directories: pi0_inc,
  dependencies:        [pi0_config_dep, ROOT_dep],
)

# Absolute path to the cut configuration. Like vz_params_file below, it is data
# rather than a build input: it is not copied into the build dir, and callers
# open it where it lives. Checked at configure time so a missing or renamed file
# is a configure error rather than a throw inside a test.
cuts_json_file = meson.project_source_root() / 'config' / 'cuts.json'

if not fs.exists(cuts_json_file)
  error('cut configuration not found at ' + cuts_json_file)
endif

# Absolute paths to the two frozen binning grids (src/core/Binning.hpp). Same
# treatment as cuts_json_file: config, not a build input, opened where it lives,
# and checked at configure time so a missing or renamed file is a configure
# error rather than a throw inside a test.
#
# BOTH ARE PLACEHOLDERS TODAY. They carry cuts.json's product edges and
# provenance.source == "placeholder"; the real ones are written by the make_grid
# tool from real data. The tests read these SHIPPED files rather than a
# hand-written copy -- a test that parses its own copy of the config proves only
# that the copy parses.
grid_a_json_file = meson.project_source_root() / 'config' / 'binning' / 'grid_A_q2_xb.json'
grid_b_json_file = meson.project_source_root() / 'config' / 'binning' / 'grid_B_z_pt2.json'

foreach grid_file : [grid_a_json_file, grid_b_json_file]
  if not fs.exists(grid_file)
    error('binning grid not found at ' + grid_file)
  endif
endforeach

# ---------------------------------------------------------------------------
# pi0_photonid -- the CLAS12 GBT photon classifier, ported from clas-framework.
#
# Uses core/Constants.hpp and ROOT's RVec, but NOT HIPO: it consumes RDataFrame
# RVec columns, so it stays unit-testable without a data file -- same rule as
# src/core.
#
# ADD NEW src/photonid/*.cpp FILES TO THIS LIST or they will not be compiled.
# src/photonid/models/*.hpp are absent by design: they are header-only (`inline`
# applicators over static tables) and are included by RunRangeModelMap.cpp
# alone. That one TU is ~25k lines of generated code and dominates this
# target's build time.
#
# warning_level=0 applies to THIS TARGET ONLY, and is not laziness:
# RunRangeModelMap.cpp includes the five models/*.hpp, which are verbatim
# CatBoost codegen under a standing "never edit" rule. Warnings nobody is
# allowed to act on are noise that hides the ones we can. Features.cpp and
# PhotonGBT.cpp are still held to warning_level=3 -- the test binary in tests/
# compiles against their headers at that level, and they carry no generated
# code of their own.
# ---------------------------------------------------------------------------
pi0_photonid_sources = files(
  'src/photonid/Features.cpp',
  'src/photonid/PhotonGBT.cpp',
  'src/photonid/RunRangeModelMap.cpp',
)

pi0_photonid = static_library(
  'pi0_photonid',
  pi0_photonid_sources,
  include_directories: pi0_inc,
  dependencies:        [ROOT_dep],  # for ROOT::VecOps::RVec only
  override_options:    ['warning_level=0'],
  install:             false,
)

pi0_photonid_dep = declare_dependency(
  link_with:            pi0_photonid,
  include_directories:  pi0_inc,
  dependencies:         [pi0_core_dep, ROOT_dep],
)

# ---------------------------------------------------------------------------
# The code version string, stamped into every output's provenance.
#
# Defined HERE rather than inside the program that needed it first. It used to
# live in src/stageA_skim/meson.build; make_grid now stamps it too, and
# make_grid's grid JSON is cross-checked against the slim files stageA_skim
# wrote. Two subdirs computing "the same" version string independently would be
# two copies of a constant -- the defect this project spends more words warning
# about than any other, and one that would drift exactly once, silently, years
# from now. Defined above both subdir() calls because meson resolves in order.
#
# The source tree is not a git checkout, so there is no commit to stamp and
# vcs_tag() would have nothing to read. The project version is what is
# available; the string says so explicitly rather than leaving a reader to
# wonder whether "0.1.0" is a tag. If this tree is ever put under version
# control, replace this with vcs_tag() and every provenance record in the
# project improves with no change to any .cxx.
# ---------------------------------------------------------------------------
pi0_code_version = meson.project_version() + ' (meson project version; source tree is not under VCS, so no commit is recorded)'

# ---------------------------------------------------------------------------
# programs and tests
#
# src/tools/meson.build may use: pi0_core_dep, pi0_config_dep, nlohmann_json_dep,
#                                hipo_dep, HipoDataFrame_dep, ROOT_dep,
#                                pi0_code_version.
#                                Note the asymmetry inside it: dump_columns is a
#                                HIPO diagnostic, make_grid reads the Stage A
#                                slim schema and takes NO hipo dependency at all.
#                                That absence is what makes "Stage A's output is
#                                a real boundary" a claim the link line checks.
# src/stageA_skim/meson.build may use: every pi0_*_dep, hipo_dep,
#                                HipoDataFrame_dep, ROOT_dep, vz_params_file,
#                                pi0_code_version.
#                                It is the one target that links them all -- the
#                                libraries are deliberately ignorant of each
#                                other, so the composition happens there.
# src/stageB_bin/meson.build may use: pi0_config_dep, pi0_core_dep,
#                                pi0_stageB_bin_dep, ROOT_dep, pi0_code_version.
#                                It builds the stageB_bin EXECUTABLE only; the
#                                pi0_stageB_bin LIBRARY's sources are listed
#                                above, because pi0_config depends on it and must
#                                therefore be declared before this subdir() runs.
#                                Like make_grid, it takes no hipo dependency and
#                                no nlohmann dependency: it reads the Stage A slim
#                                schema through ROOT, and its config through
#                                Cuts::load/Binning::load. Both absences are
#                                checkable claims rather than comments.
# tests/meson.build      may use: pi0_core_dep, pi0_config_dep, pi0_selection_dep,
#                                pi0_vertex_dep, pi0_photonid_dep, pi0_stageB_bin_dep,
#                                nlohmann_json_dep, vz_params_file,
#                                cuts_json_file, grid_a_json_file,
#                                grid_b_json_file, catch2_dep
#
# src/util/ has no target and appears in no list: Sha256.hpp is header-only, and
# every consumer already carries src/ as an include directory through one of the
# pi0_*_dep dependencies.
# ---------------------------------------------------------------------------
subdir('src/tools')
subdir('src/stageA_skim')
subdir('src/stageB_bin')

if get_option('build_tests')
  # Catch2 v3, fetched via external/catch2.wrap (same version hipo-cpp pins).
  # `catch2-with-main` provides Catch2::Catch2WithMain, so test binaries need
  # no main() of their own.
  #
  # tests=false is REQUIRED, not tidiness: Catch2's own option defaults to true,
  # and `ninja` builds every target in the build dir including subprojects', so
  # Catch2 3.14.0's SelfTest suite gets compiled as part of our build -- and it
  # does not compile under Apple clang 21 ("unknown type name 'TestInit'" out of
  # catch_template_test_registry.hpp). We want the library, not its test suite.
  # install=false keeps a vendored test framework out of our install tree.
  catch2_dep = dependency(
    'catch2-with-main',
    version: '>=3.0.0',
    default_options: {
      'tests':   false,
      'install': false,
    },
  )
  subdir('tests')
endif

# ---------------------------------------------------------------------------
summary(
  {
    'ROOT':             ROOT_dep.version(),
    'vdt include dir':  vdt_incdir,
    'ROOT_INCLUDE_PATH': root_include_path,
    'build tests':      get_option('build_tests'),
  },
  section: 'rgd_pi0_analysis',
)
