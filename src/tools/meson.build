# ---------------------------------------------------------------------------
# Assumptions, checked against the top-level meson.build as it stood when this
# file was written (not assumed blind):
#
#   * The root meson.build defines `hipo_dep`, `HipoDataFrame_dep` and
#     `ROOT_dep`. They are visible here because subdir() evaluates in the
#     parent's variable scope -- these are inherited, not re-declared.
#   * There is NO `inc_src` include_directories in the parent today, so this
#     file does not reference one; naming it would be a hard error. It is not
#     needed either: dump_columns.cxx includes only "RHipoDS.hxx", whose path
#     HipoDataFrame_dep already carries as a -I compile arg. If src/ later
#     grows shared headers and the parent defines `inc_src`, add
#     `include_directories: inc_src` below.
#
# WIRING: this file is reached only if src/meson.build does `subdir('tools')`.
# That line does not exist yet -- without it dump_columns is never built.
#
# install: false -- dump_columns is a development diagnostic for pinning the
# column-name contract, not part of the analysis pipeline. Run it from the
# build dir. Flip to true if it ever needs to ship alongside rgd_pi0.
# ---------------------------------------------------------------------------

executable(
  'dump_columns',
  sources: files('dump_columns.cxx'),
  dependencies: [hipo_dep, HipoDataFrame_dep, ROOT_dep],
  install: false,
)

# make_grid -- the equal-statistics grid generator. Its own directory because,
# unlike dump_columns, it is not a single-file diagnostic: it is the tool that
# defines the production binning, and its output is a version-controlled config
# file that every downstream result is hashed against.
#
# NOTE the dependency asymmetry with dump_columns above: make_grid takes no HIPO
# dependency at all. It reads the Stage A slim schema, which is the point --
# Stage A's output being a real boundary is a claim its link line can check.
subdir('make_grid')
