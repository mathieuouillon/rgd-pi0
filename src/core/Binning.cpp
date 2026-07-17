#include "core/Binning.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

// nlohmann/json is PRIVATE to this translation unit, exactly as it is private to
// config/Cuts.cpp: Binning.hpp does not include it, so nothing downstream
// acquires the dependency by including this grid.
//
// The same operator[] discipline applies here and for the same reason. nlohmann's
// operator[] on a non-const object DEFAULT-INSERTS null for a missing key, which
// is the silent-default failure both this file and cuts.json exist to prevent.
// Every access below goes through GridReader.

namespace pi0 {
namespace {

using nlohmann::json;

/// Navigates the grid document and throws with the full key path and the file
/// name. Deliberately a near-twin of the Reader in config/Cuts.cpp rather than a
/// shared helper: sharing it would mean src/core depends on src/config, which
/// inverts the project's dependency arrows (everything points INTO pi0_config,
/// never out of it). Forty lines of duplication is the cheaper of the two.
class GridReader {
   public:
    GridReader(const json& root, std::string file, std::string grid)
        : m_root(root), m_file(std::move(file)), m_grid(std::move(grid)) {}

    [[nodiscard]] const json& node(const std::string& path) const {
        const json* cur = &m_root;
        std::string walked;
        std::istringstream parts(path);
        std::string part;
        while (std::getline(parts, part, '.')) {
            if (!cur->is_object()) {
                fail(path, "is unreachable: '" + walked + "' is not an object");
            }
            const auto it = cur->find(part);
            walked += (walked.empty() ? "" : ".") + part;
            if (it == cur->end()) {
                fail(path, "is MISSING (no '" + walked + "')");
            }
            cur = &(*it);
        }
        return *cur;
    }

    [[nodiscard]] const json& array(const std::string& path) const {
        const json& n = node(path);
        if (!n.is_array()) fail(path, wrong_type("an array", n));
        return n;
    }

    /// A string read from an arbitrary node rather than from a path off the
    /// root -- for array elements, which node() cannot address. `display_path`
    /// must spell the full path a reader would have walked, e.g.
    /// "axes[0].name": an error naming only the leaf is an error nobody can act
    /// on.
    [[nodiscard]] std::string str_in(const json& parent, const std::string& key,
                                     const std::string& display_path) const {
        const json& n = child(parent, key, display_path);
        if (!n.is_string()) fail(display_path, wrong_type("a string", n));
        return n.get<std::string>();
    }

    [[nodiscard]] const json& array_in(const json& parent, const std::string& key,
                                       const std::string& display_path) const {
        const json& n = child(parent, key, display_path);
        if (!n.is_array()) fail(display_path, wrong_type("an array", n));
        return n;
    }

    [[noreturn]] void fail(const std::string& path, const std::string& why) const {
        throw std::runtime_error("Binning::load: grid " + m_grid + ": " + m_file + ": key '" + path +
                                 "' " + why);
    }

    [[noreturn]] void fail_file(const std::string& why) const {
        throw std::runtime_error("Binning::load: grid " + m_grid + ": " + m_file + ": " + why);
    }

   private:
    [[nodiscard]] const json& child(const json& parent, const std::string& key,
                                    const std::string& display_path) const {
        if (!parent.is_object()) fail(display_path, "is unreachable: its parent is not an object");
        const auto it = parent.find(key);
        if (it == parent.end()) fail(display_path, "is MISSING");
        return *it;
    }

    [[nodiscard]] static std::string wrong_type(const std::string& wanted, const json& got) {
        return "is not " + wanted + " (found " + std::string(got.type_name()) + ")";
    }

    const json& m_root;
    std::string m_file;
    std::string m_grid;
};

/// Parse one entry of the document's `axes` array into a validated Grid1D.
///
/// `display` is that entry's path for error messages, e.g. "axes[0]".
Grid1D read_axis(const GridReader& r, const json& node, const std::string& display) {
    Grid1D g;
    g.name = r.str_in(node, "name", display + ".name");
    if (g.name.empty()) {
        r.fail(display + ".name",
               "is the empty string. The axis name is stamped into provenance_hash() and into "
               "every error message; an unnamed axis makes both useless.");
    }

    const json& edges = r.array_in(node, "edges", display + ".edges");

    // < 2 edges is 0 bins: a grid that accepts nothing. It cannot be a typo
    // worth guessing at, so it is an error rather than an empty axis.
    if (edges.size() < 2) {
        r.fail(display + ".edges", "has " + std::to_string(edges.size()) +
                                       " edge(s). An axis needs at least 2 (i.e. at least one "
                                       "bin); fewer means a grid that bins nothing at all.");
    }

    g.edges.reserve(edges.size());
    for (std::size_t i = 0; i < edges.size(); ++i) {
        const std::string at = display + ".edges[" + std::to_string(i) + "]";
        if (!edges[i].is_number()) {
            r.fail(at, "is not a number (found " + std::string(edges[i].type_name()) + ")");
        }
        const double e = edges[i].get<double>();

        // A NaN or an infinity would defeat the monotonicity check below (every
        // comparison against NaN is false, so a NaN edge would sail through as
        // "not <=  the previous one") and then silently mis-bin everything.
        if (!std::isfinite(e)) {
            r.fail(at, "is not finite. A non-finite edge defeats the monotonicity check and then "
                       "mis-files every event that reaches it.");
        }
        g.edges.push_back(e);
    }

    // STRICTLY increasing. Equal adjacent edges are rejected too: a zero-width
    // bin can never be filled -- find() is half-open, so [e, e) is empty -- and
    // it would show up downstream as a permanently empty bin with no
    // explanation. That is a config bug, not a physics choice.
    for (std::size_t i = 1; i < g.edges.size(); ++i) {
        if (!(g.edges[i] > g.edges[i - 1])) {
            std::ostringstream os;
            os.precision(17);
            os << "is not strictly increasing: edges[" << (i - 1) << "] = " << g.edges[i - 1]
               << " but edges[" << i << "] = " << g.edges[i]
               << ". Out-of-order edges do not crash -- the binary search over them is merely "
                  "meaningless -- so every event would be filed into the wrong bin and every "
                  "number downstream would be quietly wrong. Refusing to load.";
            r.fail(display + ".edges", os.str());
        }
    }
    return g;
}

/// Load one grid file into a Grid2D.
///
/// `grid_label` is "A" or "B", for error messages only. `expect_x` / `expect_y`
/// are the axis names this grid must carry.
///
/// The axis names are CHECKED, not merely read. The whole point of the flat
/// index formula in the header is that (Q^2, x_B) is the slow-then-fast order of
/// Grid A; handing this loader a file whose axes are swapped would produce a
/// perfectly valid Binning that bins x_B on the Q^2 axis and reports it under
/// the wrong name -- and would still pass every geometry check, because the two
/// axes have the same shape.
Grid2D load_grid(const std::string& path, const std::string& grid_label,
                 const std::string& expect_x, const std::string& expect_y) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Binning::load: grid " + grid_label + ": cannot open '" + path +
                                 "'");
    }

    json root;
    try {
        in >> root;
    } catch (const json::parse_error& e) {
        throw std::runtime_error("Binning::load: grid " + grid_label + ": '" + path +
                                 "' is not valid JSON: " + e.what());
    } catch (const json::exception& e) {
        // NOT redundant with the parse_error catch above, and not defensive
        // padding. nlohmann's parser throws out_of_range (406) -- which is a
        // json::exception but NOT a json::parse_error -- for a number literal
        // that overflows a double, e.g. an edge written as 1e999. Without this
        // arm that exception escapes load() as a json::exception, and json is
        // PRIVATE to this translation unit: Binning.hpp promises std::runtime_error
        // and does not include nlohmann, so a caller following the header has
        // no type to catch and the process dies on an unhandled exception at
        // the one moment it is trying to tell someone the grid is malformed.
        throw std::runtime_error("Binning::load: grid " + grid_label + ": '" + path +
                                 "' could not be read as JSON: " + e.what());
    }

    const GridReader r(root, path, grid_label);
    if (!root.is_object()) {
        r.fail_file("is valid JSON but not an object");
    }

    const json& axes = r.array("axes");
    if (axes.size() != 2) {
        r.fail("axes", "has " + std::to_string(axes.size()) +
                           " entries; a Grid2D has exactly 2 (the slow axis first).");
    }

    Grid2D g;
    g.x = read_axis(r, axes[0], "axes[0]");
    g.y = read_axis(r, axes[1], "axes[1]");

    const auto require_named = [&](const Grid1D& axis, const std::string& expected,
                                   const std::string& display) {
        if (axis.name != expected) {
            r.fail(display + ".name",
                   "is '" + axis.name + "', but grid " + grid_label + "'s axis there must be '" +
                       expected +
                       "'. Axis ORDER is load-bearing: it is the slow-then-fast order of the flat "
                       "index formula in core/Binning.hpp, which the downstream Python decodes "
                       "against. A swapped pair would bin and index cleanly while labelling every "
                       "result with the wrong variable.");
        }
    };
    require_named(g.x, expect_x, "axes[0]");
    require_named(g.y, expect_y, "axes[1]");

    // /provenance is required but not stored: it exists so that a grid file can
    // always say where its edges came from, and the placeholders say
    // "placeholder" rather than pretending. Reading it here is what stops it
    // from becoming another declared-set-and-never-read key -- the failure mode
    // the whole config/ discipline is a monument to. It is NOT hashed: it is
    // documentation, and re-wording it must not invalidate a result.
    const json& prov = r.node("provenance");
    if (!prov.is_object()) {
        r.fail("provenance", "is not an object (found " + std::string(prov.type_name()) + ")");
    }
    if (!prov.contains("source") || !prov.at("source").is_string()) {
        r.fail("provenance.source",
               "is missing or is not a string. Every grid must say where its edges came from: "
               "'placeholder' until make_grid has run, and the scan's identity after.");
    }
    if (!prov.contains("n_events") || !prov.at("n_events").is_number_integer()) {
        r.fail("provenance.n_events",
               "is missing or is not an integer. A placeholder grid says 0; a real one says how "
               "many events its equal-statistics edges were computed from.");
    }
    return g;
}

// FNV-1a, 64-bit. Chosen over anything cryptographic because the question is
// "are these the same edges", not "did someone forge these edges", and because
// it is fifteen lines with no dependency -- provenance stamping must never be
// the reason a build needs a crypto library.
constexpr std::uint64_t kFnvOffset = 1469598103934665603ULL;
constexpr std::uint64_t kFnvPrime = 1099511628211ULL;

void hash_bytes(std::uint64_t& h, const std::string& s) {
    for (const char c : s) {
        h ^= static_cast<std::uint64_t>(static_cast<unsigned char>(c));
        h *= kFnvPrime;
    }
}

/// A double as a canonical decimal string. %.17g is the shortest format that
/// round-trips every IEEE-754 double, so two doubles serialise identically iff
/// they are the same double -- which is exactly the equivalence
/// provenance_hash() must report. snprintf rather than a stringstream so the
/// caller's locale cannot reach it and turn '.' into ',' (a hash that depends on
/// LC_NUMERIC is not a hash).
std::string canonical(double v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    return std::string(buf);
}

void hash_axis(std::uint64_t& h, const Grid1D& g) {
    hash_bytes(h, g.name);
    hash_bytes(h, "|");
    // The count is hashed, not just the values, so that concatenating two axes
    // differently cannot collide.
    hash_bytes(h, std::to_string(g.edges.size()));
    for (const double e : g.edges) {
        hash_bytes(h, ",");
        hash_bytes(h, canonical(e));
    }
    hash_bytes(h, ";");
}

}  // namespace

// ---------------------------------------------------------------------------
// Grid1D
// ---------------------------------------------------------------------------

int Grid1D::nbins() const {
    return edges.size() < 2 ? 0 : static_cast<int>(edges.size()) - 1;
}

int Grid1D::find(double v) const {
    if (edges.size() < 2) return -1;

    // Written as !(v >= lo) rather than (v < lo) so that a NaN takes the -1
    // branch: every comparison against NaN is false, so !(NaN >= lo) is true.
    // A NaN reaching upper_bound below would be a value with no ordering
    // against the range, i.e. an arbitrary bin.
    if (!(v >= edges.front())) return -1;
    if (!(v <= edges.back())) return -1;

    // THE TOP EDGE. A value exactly on edges.back() is in range but has no
    // half-open bin: [lo, hi) never contains hi. The old find_1d_bin() let it
    // fall through and returned -1, which silently DROPPED every p_T^2 sitting
    // exactly on the top edge. Put it in the last bin instead -- the top edge is
    // the kinematic limit, and discarding the events that reach it is not a
    // defensible reading of "outside the grid".
    if (v == edges.back()) return nbins() - 1;

    // upper_bound gives the first edge strictly greater than v; the bin is the
    // one starting at the edge before it. Half-open [lo, hi) falls out of this
    // directly: v == some interior edge yields that edge's own bin.
    const auto it = std::upper_bound(edges.begin(), edges.end(), v);
    return static_cast<int>(it - edges.begin()) - 1;
}

// ---------------------------------------------------------------------------
// Grid2D
// ---------------------------------------------------------------------------

int Grid2D::find(double vx, double vy) const {
    const int ix = x.find(vx);
    if (ix < 0) return -1;
    const int iy = y.find(vy);
    if (iy < 0) return -1;
    return ix * y.nbins() + iy;  // row-major, y fast. See Binning's header doc.
}

int Grid2D::ncells() const { return x.nbins() * y.nbins(); }

// ---------------------------------------------------------------------------
// Binning
// ---------------------------------------------------------------------------

int Binning::find_4d(double q2, double xb, double z, double pt2) const {
    const int a = A.find(q2, xb);
    if (a < 0) return -1;
    const int b = B.find(z, pt2);
    if (b < 0) return -1;
    return a * B.ncells() + b;  // B cell fast.
}

int Binning::n4d() const { return A.ncells() * B.ncells(); }

int Binning::find_3d(double q2, double xb, double z) const {
    const int a = A.find(q2, xb);
    if (a < 0) return -1;
    const int iz = B.x.find(z);
    if (iz < 0) return -1;
    return a * B.x.nbins() + iz;  // z fast. p_T^2 is the observable, not a bin.
}

int Binning::n3d() const { return A.ncells() * B.x.nbins(); }

std::string Binning::provenance_hash() const {
    std::uint64_t h = kFnvOffset;
    // The grid tags are hashed too, so that swapping Grid A's file for Grid B's
    // is a different binning rather than the same multiset of edges.
    hash_bytes(h, "A:");
    hash_axis(h, A.x);
    hash_axis(h, A.y);
    hash_bytes(h, "B:");
    hash_axis(h, B.x);
    hash_axis(h, B.y);

    char buf[32];
    std::snprintf(buf, sizeof(buf), "%016llx", static_cast<unsigned long long>(h));
    return std::string(buf);
}

Binning Binning::load(const std::string& grid_a_json, const std::string& grid_b_json) {
    Binning b;
    b.A = load_grid(grid_a_json, "A", "q2", "xb");
    b.B = load_grid(grid_b_json, "B", "z", "pt2");
    return b;
}

}  // namespace pi0
