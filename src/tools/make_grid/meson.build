# ---------------------------------------------------------------------------
# make_grid -- Stage A slim files -> config/binning/grid_{A,B}*.json.
#
# Inherited from the parent scope via subdir() (meson evaluates subdir() in the
# caller's variable scope, so these are inherited, not re-declared):
# pi0_config_dep, pi0_core_dep, nlohmann_json_dep, ROOT_dep.
#
# NO HIPO. This program reads the slim ROOT schema, never a raw HIPO file, so
# neither hipo_dep nor HipoDataFrame_dep belongs on this link line -- and their
# absence is the checkable statement that Stage A's output really is a boundary
# rather than a convention.
#
# nlohmann_json_dep DIRECTLY, and this does not break the rule the top-level
# meson.build states ("nothing else in the project may use it to read a cut --
# Cuts::load() is the only reader of cuts.json"). make_grid WRITES a grid file;
# it reads its cuts through Cuts::load like everything else. tests/meson.build
# takes the same dependency for the same kind of reason. The rule is about who
# reads cuts.json, not about who may link a JSON library.
#
# PI0_CODE_VERSION is injected exactly as src/stageA_skim does, and for the same
# reason: the emitted grid's provenance has to name the code that produced it.
# The value is defined there; recomputing the same string here would be a second
# copy of a constant, so it is read out of the parent scope instead. That is why
# this subdir() must come after src/stageA_skim's in the top-level meson.build.
#
# install: false, matching the rest of the project. Run it from the build dir.
# ---------------------------------------------------------------------------

executable(
  'make_grid',
  sources: files('main.cxx'),
  dependencies: [
    pi0_config_dep,     # Cuts::load -- the grid shape, the pairing cuts, the z window
    pi0_core_dep,       # find_gg_pairs, compute_sidis, Types
    nlohmann_json_dep,  # writing the grid JSON (NOT reading a cut)
    ROOT_dep,           # TFile, TTree, TTreeReader, TNamed, ROOT::EnableThreadSafety
    thread_dep,         # std::thread pool for opt-in --threads sample collection (util/Progress.hpp too)
  ],
  cpp_args: [
    '-DPI0_CODE_VERSION="' + pi0_code_version + '"',
  ],
  install: false,
  # Held to the strict level like every other piece of hand-written code here.
  override_options: ['warning_level=3'],
)
