#pragma once

/// \file Provenance.hpp
/// \brief The provenance block every stage writes into its output, and reads out
///        of its input. Header-only.
///
/// WHY THIS IS SHARED RATHER THAN COPIED INTO EACH PROGRAM
/// ------------------------------------------------------
/// It began as a local struct in stageA_skim/main.cxx, which was fine while
/// exactly one program wrote one. It is here now because stageB_bin does not
/// merely WRITE a provenance block -- it READS Stage A's out of the slim file and
/// PROPAGATES IT FORWARD. That makes the layout (a `provenance` TDirectory of
/// TNamed(key -> value), plus a TObjString of the whole block as text) a contract
/// BETWEEN two programs rather than a detail inside one, and a contract with two
/// independent implementations is a contract that holds until somebody edits one
/// of them. Same argument as util/Sha256.hpp, which was extracted for the same
/// reason on the same day: a fingerprint whose value is that everyone computes it
/// identically may not have two implementations.
///
/// The chain this exists to protect: the analysis this project replaces could not
/// reproduce its own production, because nothing recorded which cuts, which
/// model, or which code had made a given file -- so a plot could not be traced
/// back to the thing that made it, and by the time anyone noticed, the inputs
/// were gone. Every stage's output must carry everything needed to re-run it AND
/// everything its inputs said about themselves. A stage that records its own
/// parameters but drops its input's has merely moved the hole one step downstream.
///
/// Header-only (`inline`) so it needs no meson target: every consumer already
/// carries `src/` as an include directory and already links ROOT.

#include <array>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <TDirectory.h>
#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TNamed.h>
#include <TObjString.h>

namespace pi0::util {

/// An ordered list of (key, value) strings, written into a ROOT file.
///
/// ORDERED, and a vector rather than a map on purpose: the block is read by
/// humans as often as by programs, and the order it is written in is the order it
/// makes sense in. A std::map would sort "config.sha256" above "run" and scatter
/// the fields that belong together.
struct Provenance {
    std::vector<std::pair<std::string, std::string>> entries;

    void add(std::string key, std::string value) { entries.emplace_back(std::move(key), std::move(value)); }

    /// Write into `file` as a TDirectory named `dir_name` of TNamed(key -> value).
    ///
    /// \param dir_name  the directory to create. Defaults to "provenance", which
    ///                  is where a reader looks for a file's OWN provenance. A
    ///                  stage propagating its input's block forward writes it
    ///                  under a different name, so that an inherited value can
    ///                  never be mistaken for one this stage stands behind.
    /// \throws std::runtime_error if the directory cannot be created.
    void write(TFile& file, const char* dir_name = "provenance") const {
        TDirectory* dir = file.mkdir(dir_name);
        if (dir == nullptr) {
            throw std::runtime_error(std::string("cannot create the '") + dir_name +
                                     "' directory in the output file");
        }
        dir->cd();
        for (const auto& [k, v] : entries) {
            TNamed n(k.c_str(), v.c_str());
            n.Write();
        }
        file.cd();
    }

    /// Write the whole block as one TObjString, for humans.
    ///
    /// Both forms on purpose: the TNamed directory is what a program greps for a
    /// single field, this is what a person prints. uproot reads both.
    void write_text(TFile& file, const char* obj_name = "provenance_text") const {
        file.cd();
        TObjString text(as_text().c_str());
        text.Write(obj_name);
    }

    [[nodiscard]] std::string as_text() const {
        std::size_t w = 0;
        for (const auto& [k, v] : entries) {
            (void)v;
            w = std::max(w, k.size());
        }
        std::ostringstream os;
        for (const auto& [k, v] : entries) os << std::left << std::setw(static_cast<int>(w)) << k << " : " << v << '\n';
        return os.str();
    }

    /// Read a provenance directory back out of `file`, IN WRITE ORDER.
    ///
    /// Order comes from the directory's key list, which ROOT keeps in the order
    /// the objects were written. That matters for propagation: a block that
    /// arrives shuffled is still correct and still unreadable.
    ///
    /// \return an empty Provenance if `dir_name` does not exist. That is not an
    ///         error here -- deciding whether a file without provenance is
    ///         acceptable is the caller's judgement, not this reader's, and a
    ///         caller that cares can test `entries.empty()`. What this must never
    ///         do is invent a value.
    [[nodiscard]] static Provenance read(TFile& file, const char* dir_name = "provenance") {
        Provenance p;
        TDirectory* dir = file.GetDirectory(dir_name);
        if (dir == nullptr) return p;

        TList* keys = dir->GetListOfKeys();
        if (keys == nullptr) return p;

        TIter next(keys);
        while (TObject* o = next()) {
            auto* key = dynamic_cast<TKey*>(o);
            if (key == nullptr) continue;
            // Read through the key rather than dir->Get<TNamed>(name): a key list
            // preserves write order and can carry cycles, and Get() by name would
            // silently collapse duplicates onto the newest.
            auto* named = dynamic_cast<TNamed*>(key->ReadObj());
            if (named == nullptr) continue;
            p.add(named->GetName(), named->GetTitle());
        }
        return p;
    }

    /// \return the value for `key`, or "" if absent. "" is also a legal value, so
    ///         a caller that must distinguish the two should search `entries`.
    [[nodiscard]] std::string get(const std::string& key) const {
        for (const auto& [k, v] : entries) {
            if (k == key) return v;
        }
        return {};
    }
};

/// The current UTC time as an ISO-8601 string.
[[nodiscard]] inline std::string utc_now() {
    const std::time_t t = std::time(nullptr);
    std::tm tm{};
#if defined(_WIN32)
    gmtime_s(&tm, &t);
#else
    gmtime_r(&t, &tm);
#endif
    std::array<char, 32> buf{};
    std::strftime(buf.data(), buf.size(), "%Y-%m-%dT%H:%M:%SZ", &tm);
    return std::string(buf.data());
}

}  // namespace pi0::util
