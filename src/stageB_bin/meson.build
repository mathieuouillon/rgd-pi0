# ---------------------------------------------------------------------------
# stageB_bin -- the slim TTree -> binned spectra pass.
#
# Inherited from the parent scope via subdir() (meson evaluates subdir() in the
# caller's variable scope, so these are inherited, not re-declared):
# pi0_config_dep, pi0_core_dep, pi0_stageB_bin_dep, ROOT_dep, pi0_code_version.
#
# NOTE this directory is reached by TWO things in the top-level meson.build, and
# they are not the same thing:
#   * the `pi0_stageB_bin` STATIC LIBRARY, whose sources (DonorPool.cpp,
#     PoolGrid.cpp) are listed up there because pi0_config depends on it for
#     MixingCuts -- it must be declared before pi0_config, i.e. before this
#     subdir() could possibly run;
#   * this file, which builds the EXECUTABLE and must run after pi0_config_dep
#     exists.
# Hence the library's sources are named in the parent and only main.cxx is named
# here. Do not "tidy" the library sources into this file: it would move them
# after pi0_config and the configure would fail on an undefined variable.
#
# NO HIPO, and the absence is the point -- the same claim make_grid's link line
# makes. This program reads the Stage A slim ROOT schema and never a raw HIPO
# file, so Stage A's output being a real boundary is something the build checks
# rather than something a comment asserts.
#
# nlohmann_json_dep is likewise ABSENT. This program reads two grid JSONs and one
# cuts.json, and it parses neither: Binning::load() and Cuts::load() are the only
# readers of those files in the project, and they arrive through pi0_core_dep and
# pi0_config_dep with nlohmann private to them. A JSON dependency on this line
# would be the first sign someone had started parsing config here.
#
# PI0_CODE_VERSION is injected exactly as src/stageA_skim and src/tools/make_grid
# do, from the parent's pi0_code_version. Recomputing the same string here would
# be a second copy of a constant.
#
# install: false, matching the rest of the project. Run it from the build dir.
# ---------------------------------------------------------------------------

executable(
  'stageB_bin',
  sources: files('main.cxx'),
  dependencies: [
    pi0_config_dep,      # Cuts::load -- the cuts, the m_gg axis, the phi_h axis, the mixing block
    pi0_core_dep,        # Binning::load, find_gg_pairs, admissible_pair, compute_sidis, Grid1D, NeumaierSum
    pi0_stageB_bin_dep,  # DonorPool, DonorPhoton
    ROOT_dep,            # TFile, TTree, TChain, TTreeReader, TNamed, TObjString
    thread_dep,          # std::thread pool for the opt-in multi-threaded pass 1; util/Progress.hpp
  ],
  cpp_args: [
    '-DPI0_CODE_VERSION="' + pi0_code_version + '"',
  ],
  install: false,
  # Held to the strict level like every other piece of hand-written code here.
  override_options: ['warning_level=3'],
)
