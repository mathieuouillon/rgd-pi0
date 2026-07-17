# ---------------------------------------------------------------------------
# stageA_skim -- the HIPO -> slim ROOT TTree pass.
#
# This is the only target in the project that links EVERYTHING: it is where the
# five independently-testable libraries are finally composed against real data.
# Each of them is deliberately ignorant of the others (selection does not know
# about the GBT or the vertex; photonid does not know about Cuts), so the wiring
# has to live somewhere, and this is it.
#
# Inherited from the parent scope via subdir(): pi0_core_dep, pi0_config_dep,
# pi0_selection_dep, pi0_photonid_dep, pi0_vertex_dep, hipo_dep,
# HipoDataFrame_dep, ROOT_dep, vz_params_file.
#
# install: false, matching src/tools. Nothing in this project installs yet; the
# skim is run from the build dir. Flip this (and add an install_dir) when there
# is a production deployment to install into.
# ---------------------------------------------------------------------------

# The vz parameters are DATA, not a build input: the file is not copied into the
# build dir, and the program opens it where it lives. Injected as a macro rather
# than defaulted in C++ so that a moved or renamed params file is a link-line
# change here, not a silent runtime fallback. main.cxx #errors if it is absent,
# so there is no path on which the skim runs without the correction. Same
# mechanism, and the same reasoning, as tests/meson.build.
#
# The parent has already checked that vz_params_file exists (configure-time
# error if not), so this cannot inject a path to nothing.

# pi0_code_version is defined in the TOP-LEVEL meson.build, not here. It moved
# there when make_grid started stamping it too: two subdirs deriving "the same"
# version string independently is two copies of a constant. Inherited via
# subdir(), like the pi0_*_dep variables above.

executable(
  'stageA_skim',
  sources: files('main.cxx'),
  dependencies: [
    pi0_config_dep,     # Cuts::load -- every threshold
    pi0_selection_dep,  # find_trigger_electron, pass_electron, pass_photon_scored
    pi0_photonid_dep,   # CaloMap, build_features, score, select_model
    pi0_vertex_dep,     # VzCorrector: correct() + pass_window()
    pi0_core_dep,       # compute_dis, angle_between_deg, Constants, NeumaierSum/partition_*
    hipo_dep,
    HipoDataFrame_dep,  # RHipoDS
    ROOT_dep,           # RDataFrame, TFile, TTree, TNamed, TObjString, EnableThreadSafety
    thread_dep,         # std::thread pool for the opt-in multi-threaded (--threads) skim; util/Progress.hpp
  ],
  cpp_args: [
    '-DPI0_VZ_PARAMS_FILE="' + vz_params_file + '"',
    '-DPI0_CODE_VERSION="' + pi0_code_version + '"',
  ],
  install: false,
  # Held to the strict level like the rest of the hand-written code. The
  # warning_level=0 override on pi0_photonid is scoped to the target that
  # compiles the verbatim CatBoost codegen; nothing generated is included here
  # (the models arrive as a function pointer from select_model()).
  override_options: ['warning_level=3'],
)
