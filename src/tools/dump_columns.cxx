// ---------------------------------------------------------------------------
// dump_columns -- pin down the RDataFrame column-name contract for HIPO files.
//
// Why this tool exists:
//   Everything downstream in this project reaches HIPO bank data through
//   RDataFrame *column names*, and those are NOT the names you see in
//   `hipo-utils -dump`. RHipoDS rewrites them: the HIPO column
//   `REC::Particle.px` is reachable from RDataFrame only as `REC_Particle_px`.
//   Getting that wrong yields a runtime "Unknown column" error, not a compile
//   error -- so the mapping is worth observing against a real file rather than
//   remembering. This tool makes the contract checkable.
//
// It prints (a) the entry count, (b) each column with its exact C++ type and
// whether it is a vector, and (c) the raw -> translated naming rule, derived by
// diffing the two name tables RHipoDS keeps rather than hardcoded here.
//
// Two things that look like bugs but are not:
//   * RHipoDS calls setTags(0) internally. Tag 0 is the PHYSICS tag, so the
//     entry count below counts ONLY tag-0 events. A file may hold many
//     non-tag-0 events (scaler, calibration, ...); they are excluded on
//     purpose. A count lower than `hipo-utils -info` reports is expected.
//   * Multi-row banks surface as ROOT::VecOps::RVec<T> (one entry per particle
//     per event), single-row banks as bare scalars. The SHAPE column says
//     which; read it before writing a Define() that assumes one or the other.
//
// The rest of the project builds its RDataFrame via either
//   auto ds = std::make_unique<RHipoDS>(file, n_inspect);
//   auto df = ROOT::RDataFrame(std::move(ds));
// or the MakeHipoDataFrame(file) convenience wrapper. This tool deliberately
// stops at the datasource: constructing the RDataFrame would std::move(ds) and
// leave nothing to interrogate, and the column contract lives on the
// datasource anyway.
//
// usage: dump_columns <file.hipo> [name-prefix]        (prefix default: REC_)
// ---------------------------------------------------------------------------

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include "RHipoDS.hxx"

namespace {
// C++17: no std::string::starts_with (that is C++20, and this project pins
// cpp_std=c++17 to match ROOT). rfind(p, 0) == 0 anchors the search at 0.
bool starts_with(const std::string& s, const std::string& prefix) {
    return s.rfind(prefix, 0) == 0;
}
}  // namespace

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " <file.hipo> [name-prefix]\n"
                  << "  name-prefix  filter for the column table; defaults to \"REC_\".\n"
                  << "               Pass \"\" to list every column.\n";
        return 1;
    }
    const std::string file = argv[1];
    // argc > 2 rather than argv[2] != nullptr: an explicit empty prefix ("")
    // is meaningful here -- it means "list everything".
    const std::string prefix = (argc > 2) ? argv[2] : "REC_";

    // n_inspect = how many events RHipoDS scans up front to learn the schema
    // and decide vector-vs-scalar per column. 10000 is the RHipoDS default.
    auto ds = std::make_unique<RHipoDS>(file, 10000);

    const std::vector<std::string>& cols = ds->GetColumnNames();

    std::cout << "file             : " << file << '\n'
              << "tag-0 entries    : " << ds->GetEntries()
              << "   (PHYSICS tag only -- non-tag-0 events excluded by design)\n"
              << "columns          : " << cols.size() << '\n'
              << "name translation : " << (ds->fColumnNameTranslation ? "ON" : "OFF") << '\n';

    // --- naming convention, derived from this file -------------------------
    // RHipoDS stores the raw HIPO names in fAllColumnsPreTranslated and the
    // RDataFrame-visible names in fAllColumns, index-for-index. Diffing them
    // shows the rule as it actually applied here, instead of asserting it.
    std::cout << "\n--- naming convention ---\n";
    if (!ds->fColumnNameTranslation) {
        std::cout << "Translation is OFF: RDataFrame column names are the raw HIPO names.\n";
    } else {
        std::cout << "HIPO \"::\" and \".\" each become \"_\" in the RDataFrame column name.\n"
                  << "Observed in this file:\n";
        int shown = 0;
        for (std::size_t i = 0; i < cols.size() && shown < 3; ++i) {
            if (ds->fAllColumnsPreTranslated[i] == cols[i]) continue;  // name needed no rewrite
            std::cout << "  HIPO " << std::left << std::setw(30) << ds->fAllColumnsPreTranslated[i]
                      << " ->  RDataFrame " << cols[i] << '\n';
            ++shown;
        }
        if (shown == 0) std::cout << "  (no column name in this file required rewriting)\n";
    }

    // --- column table ------------------------------------------------------
    std::cout << "\n--- columns matching prefix \"" << prefix << "\" ---\n"
              << std::left << std::setw(34) << "NAME" << std::setw(30) << "TYPE" << "SHAPE\n";

    std::size_t n_match = 0;
    for (std::size_t i = 0; i < cols.size(); ++i) {
        if (!starts_with(cols[i], prefix)) continue;
        ++n_match;
        // Look the type up BY NAME, not by index: that exercises the exact
        // name -> type path the rest of the project relies on, so a name that
        // prints here is a name that resolves.
        std::cout << std::left << std::setw(34) << cols[i] << std::setw(30) << ds->GetTypeName(cols[i])
                  << (ds->fColumnTypeIsVector[i] ? "VECTOR" : "SCALAR") << '\n';
    }

    std::cout << '\n' << n_match << " of " << cols.size() << " columns match \"" << prefix << "\"";
    if (n_match == 0) std::cout << "   (no match -- try another prefix, or \"\" to list all)";
    std::cout << '\n';

    return 0;
}
