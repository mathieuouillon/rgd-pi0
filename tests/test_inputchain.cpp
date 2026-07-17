// Unit tests for the multi-input chain (src/stageB_bin/InputChain.hpp).
//
// --input repeats because the donor pool is only well filled when ONE Stage B
// run sees ALL of a target's slims (~30,866 files for LD2; one file fills 33 of
// the 224 shipped pool bins). These tests exist because chaining files makes four
// specific silent corruptions newly possible, and every one of them produces an
// output that looks entirely ordinary:
//
//   1. MIXED TARGETS. Nothing downstream of the read knows a photon's target -- a
//      DonorPhoton is four floats -- so an LD2 file chained with an Sn file
//      builds ONE pool in which Sn photons are the mixed background subtracted
//      from LD2 events. It corrupts the numerator AND the denominator of R_A. No
//      plot shows it. It must REFUSE, and the test is that it does.
//   2. A TRUNCATED INPUT. Stage A stamps `events.max_events_requested` saying its
//      output is a PREFIX and "no yield or normalisation from this file is
//      complete". Nothing has ever read it. Chained with full files, the missing
//      N_DIS is invisible while the spectra look complete.
//   3. AN ORDER-DEPENDENT POOL. Reservoir sampling is a function of the SEQUENCE
//      of offers (DonorPool.hpp), so `--input a --input b` and `--input b --input
//      a` would otherwise build different pools from the same data -- and a farm
//      glob promises no order.
//   4. A PER-FILE --max-events. Pass 1 truncates with RDataFrame::Range, which is
//      already chain-wide; pass 0 counts by hand. A budget scoped per file makes
//      the pool model a strictly larger sample than the spectra.
//
// WHAT IS TESTED HERE AND WHAT IS NOT. These are the pure decisions -- ordering,
// refusal, the seed, the budget arithmetic -- so they need no ROOT fixture and
// none is faked: a Provenance is a vector of pairs, so a chain is built in code.
// What CANNOT be reached from here is the wiring in stageB_bin/main.cxx: that the
// budget is threaded BY REFERENCE across the per-file loop, and that pass 1 hands
// RDataFrame the same canonical order. Those are checked at runtime instead, by
// the p0_events == n_events cross-check after both passes -- which is exactly the
// assertion a per-file budget fails -- and end-to-end against real slim files.
//
// Method note, as in test_donorpool: every digest here is a literal. The seed is
// a pure function of them, so an assertion about it is a claim about ONE known
// input and cannot flake.

#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "stageB_bin/InputChain.hpp"
#include "util/Provenance.hpp"

using namespace pi0::stageB;
using pi0::util::Provenance;
using Catch::Matchers::ContainsSubstring;

namespace {

/// A full, untruncated Stage A block for `target`, as stageA_skim writes one.
/// The fields are the ones validate_chain cross-checks; the values are Stage A's
/// real strings rather than plausible-looking ones, because two of the checks
/// turn on their exact text ("TRUNCATED", "all (...").
Provenance stagea(const std::string& target, const std::string& run = "22083",
                  const std::string& config_sha = "2aff789e3095c6e390502231e8ea794344fae06a96155fb4576308b905bdad35") {
    Provenance p;
    p.add("config.sha256", config_sha);
    p.add("target", target);
    p.add("polarity", "outbending");
    p.add("beam.energy_gev", "10.53");
    p.add("run", run);
    p.add("gbt.model", "RGA_inbending_pass1");
    p.add("gbt.fallback_used", "TRUE -- photons scored by a model trained on OTHER data");
    p.add("events.max_events_requested", "all (no --max-events; the whole input was read)");
    return p;
}

/// The same, but truncated -- Stage A's own wording on a --max-events run.
Provenance stagea_truncated(const std::string& target) {
    Provenance p = stagea(target);
    for (auto& [k, v] : p.entries) {
        if (k != "events.max_events_requested") continue;
        v = "20000 -- TRUNCATED RUN: --max-events was given, so events.input_tag0 below is a PREFIX of the "
            "input, not the whole of it, and no yield or normalisation from this file is complete";
    }
    return p;
}

/// A 64-hex-character digest that is `c` repeated. Real sha256 digests are
/// unreadable and the tests are about their ORDER, so they are made legible.
std::string digest(char c) { return std::string(64, c); }

ChainInput input(const std::string& path, char digest_char, const Provenance& p) {
    ChainInput ci;
    ci.path = path;
    ci.sha256 = digest(digest_char);
    ci.prov = p;
    return ci;
}

}  // namespace

// ===========================================================================
// the target check -- the one that corrupts R_A silently
// ===========================================================================

TEST_CASE("a chain of mixed targets is REFUSED", "[inputchain][validate]") {
    const std::vector<ChainInput> chain = {input("/data/ld2.root", 'a', stagea("LD2")),
                                           input("/data/sn.root", 'b', stagea("Sn"))};

    REQUIRE_THROWS_WITH(validate_chain(chain, false), ContainsSubstring("target"));
    // The message must name BOTH files and BOTH values: "the targets disagree" is
    // not actionable on a 30,866-file command line.
    REQUIRE_THROWS_WITH(validate_chain(chain, false), ContainsSubstring("/data/ld2.root"));
    REQUIRE_THROWS_WITH(validate_chain(chain, false), ContainsSubstring("/data/sn.root"));
    REQUIRE_THROWS_WITH(validate_chain(chain, false), ContainsSubstring("LD2"));
    REQUIRE_THROWS_WITH(validate_chain(chain, false), ContainsSubstring("Sn"));
    // --allow-truncated-inputs is about truncation and MUST NOT be a general
    // "accept this chain" escape hatch.
    REQUIRE_THROWS_AS(validate_chain(chain, true), std::runtime_error);
}

TEST_CASE("a chain of ONE target is accepted", "[inputchain][validate]") {
    const std::vector<ChainInput> chain = {input("/data/ld2_a.root", 'a', stagea("LD2", "22083")),
                                           input("/data/ld2_b.root", 'b', stagea("LD2", "22084")),
                                           input("/data/ld2_c.root", 'c', stagea("LD2", "22085"))};
    // Three different RUNS of one target is the NORMAL case and the entire reason
    // --input repeats. `run` must not be a chain invariant.
    REQUIRE_NOTHROW(validate_chain(chain, false));
}

TEST_CASE("every field that would become a false single value is refused", "[inputchain][validate]") {
    const auto chain_differing_in = [](const std::string& field, const std::string& other) {
        Provenance b = stagea("LD2");
        for (auto& [k, v] : b.entries) {
            if (k == field) v = other;
        }
        return std::vector<ChainInput>{input("/data/a.root", 'a', stagea("LD2")), input("/data/b.root", 'b', b)};
    };

    SECTION("config.sha256 -- different cuts selected the photons sharing one pool") {
        REQUIRE_THROWS_WITH(validate_chain(chain_differing_in("config.sha256", digest('9')), false),
                            ContainsSubstring("config.sha256"));
    }
    SECTION("gbt.model -- a different photon purity in one sample") {
        REQUIRE_THROWS_WITH(validate_chain(chain_differing_in("gbt.model", "RGC_Summer2022_pass1"), false),
                            ContainsSubstring("gbt.model"));
    }
    SECTION("gbt.fallback_used") {
        REQUIRE_THROWS_WITH(validate_chain(chain_differing_in("gbt.fallback_used", "false -- ..."), false),
                            ContainsSubstring("gbt.fallback_used"));
    }
    SECTION("beam.energy_gev -- compute_sidis reads it") {
        REQUIRE_THROWS_WITH(validate_chain(chain_differing_in("beam.energy_gev", "10.6"), false),
                            ContainsSubstring("beam.energy_gev"));
    }
    SECTION("polarity") {
        REQUIRE_THROWS_WITH(validate_chain(chain_differing_in("polarity", "inbending"), false),
                            ContainsSubstring("polarity"));
    }
}

TEST_CASE("a field one input declares and another omits cannot be shown to agree",
          "[inputchain][validate]") {
    // "" is a legal provenance value and Provenance::get() returns "" for an
    // absent key too, so a validator built on get() would read this as agreement
    // between a file that says LD2 and a file that says nothing.
    Provenance no_target;
    no_target.add("config.sha256", digest('2'));
    no_target.add("polarity", "outbending");

    const std::vector<ChainInput> chain = {input("/data/a.root", 'a', stagea("LD2")),
                                           input("/data/b.root", 'b', no_target)};
    REQUIRE_THROWS_WITH(validate_chain(chain, false), ContainsSubstring("<absent>"));
}

TEST_CASE("a chain in which NOBODY names a target is refused", "[inputchain][validate]") {
    // Every input agreeing that nothing is known is not a chain whose targets
    // were checked. One input is the old, warned-about "no provenance" case and
    // stays allowed; many is the R_A corruption with the evidence missing.
    Provenance blank;
    const std::vector<ChainInput> one = {input("/data/a.root", 'a', blank)};
    REQUIRE_NOTHROW(validate_chain(one, false));

    const std::vector<ChainInput> two = {input("/data/a.root", 'a', blank), input("/data/b.root", 'b', blank)};
    REQUIRE_THROWS_WITH(validate_chain(two, false), ContainsSubstring("target"));
}

// ===========================================================================
// truncation
// ===========================================================================

TEST_CASE("a TRUNCATED Stage A input is REFUSED", "[inputchain][truncation]") {
    SECTION("alone") {
        const std::vector<ChainInput> chain = {input("/data/smoke.root", 'a', stagea_truncated("LD2"))};
        REQUIRE_THROWS_WITH(validate_chain(chain, false), ContainsSubstring("TRUNCATED"));
        REQUIRE_THROWS_WITH(validate_chain(chain, false), ContainsSubstring("/data/smoke.root"));
    }
    SECTION("hidden among full files -- where N_DIS goes quietly wrong") {
        const std::vector<ChainInput> chain = {input("/data/full_a.root", 'a', stagea("LD2")),
                                               input("/data/smoke.root", 'b', stagea_truncated("LD2")),
                                               input("/data/full_c.root", 'c', stagea("LD2"))};
        REQUIRE_THROWS_WITH(validate_chain(chain, false), ContainsSubstring("/data/smoke.root"));
    }
    SECTION("--allow-truncated-inputs accepts it") {
        const std::vector<ChainInput> chain = {input("/data/smoke.root", 'a', stagea_truncated("LD2"))};
        REQUIRE_NOTHROW(validate_chain(chain, true));
    }
}

TEST_CASE("a full Stage A input is not mistaken for a truncated one", "[inputchain][truncation]") {
    const std::vector<ChainInput> chain = {input("/data/full.root", 'a', stagea("LD2"))};
    REQUIRE_NOTHROW(validate_chain(chain, false));
}

TEST_CASE("an input with no truncation flag at all is accepted", "[inputchain][truncation]") {
    // Absent is UNKNOWABLE, not false: a file written before the flag existed, or
    // by something other than this project's stageA_skim. Refusing it here would
    // reject files this program has always accepted. The existing "no provenance"
    // warning covers it.
    Provenance p = stagea("LD2");
    Provenance without;
    for (const auto& [k, v] : p.entries) {
        if (k != "events.max_events_requested") without.add(k, v);
    }
    const std::vector<ChainInput> chain = {input("/data/old.root", 'a', without)};
    REQUIRE_NOTHROW(validate_chain(chain, false));
}

// ===========================================================================
// the canonical order
// ===========================================================================

TEST_CASE("the chain is ordered by CONTENT, whatever order it was typed in",
          "[inputchain][order]") {
    std::vector<ChainInput> forward = {input("/z/first.root", 'a', stagea("LD2")),
                                       input("/a/second.root", 'b', stagea("LD2")),
                                       input("/m/third.root", 'c', stagea("LD2"))};
    std::vector<ChainInput> reversed = {forward[2], forward[1], forward[0]};

    order_canonically(forward);
    order_canonically(reversed);

    REQUIRE(forward.size() == reversed.size());
    for (std::size_t i = 0; i < forward.size(); ++i) {
        REQUIRE(forward[i].sha256 == reversed[i].sha256);
        REQUIRE(forward[i].path == reversed[i].path);
    }
    REQUIRE(forward[0].sha256 == digest('a'));
    // NOT sorted by path: "/a/second.root" would come first if it were, and the
    // pool would then depend on where somebody put the files -- the same defect
    // the seed refuses to have by hashing content rather than paths.
    REQUIRE(forward[0].path == "/z/first.root");
}

TEST_CASE("two inputs with identical content are REFUSED", "[inputchain][order]") {
    // The same events twice: N_DIS doubles and every yield doubles, so the ratio
    // of the two looks perfectly correct while the normalisation under it is not.
    std::vector<ChainInput> chain = {input("/data/a.root", 'a', stagea("LD2")),
                                     input("/data/a_copy.root", 'a', stagea("LD2"))};
    REQUIRE_THROWS_WITH(order_canonically(chain), ContainsSubstring("identical content"));
    REQUIRE_THROWS_WITH(order_canonically(chain), ContainsSubstring("/data/a_copy.root"));
}

// ===========================================================================
// the seed
// ===========================================================================

TEST_CASE("the seed does not depend on the command-line order", "[inputchain][seed]") {
    // THE PROPERTY THAT MAKES A FARM GLOB SAFE. The pool is a function of the
    // sequence of offers, so a seed that moved with the argv order would give two
    // backgrounds for one dataset.
    std::vector<ChainInput> forward = {input("/z/one.root", 'c', stagea("LD2")),
                                       input("/a/two.root", 'a', stagea("LD2")),
                                       input("/m/three.root", 'b', stagea("LD2"))};
    std::vector<ChainInput> reversed = {forward[2], forward[1], forward[0]};
    std::vector<ChainInput> shuffled = {forward[1], forward[0], forward[2]};

    order_canonically(forward);
    order_canonically(reversed);
    order_canonically(shuffled);

    const std::string d = chain_digest(forward);
    REQUIRE(chain_digest(reversed) == d);
    REQUIRE(chain_digest(shuffled) == d);
    REQUIRE(seed_from_digest(chain_digest(reversed)) == seed_from_digest(d));
    REQUIRE(seed_from_digest(chain_digest(shuffled)) == seed_from_digest(d));
}

TEST_CASE("a one-input chain's digest IS that input's digest", "[inputchain][seed]") {
    // NOT a special case for its own sake: re-hashing a lone digest would move
    // the seed, hence the pool, hence the mixed background of every result made
    // before --input could repeat. A single-input run must stay bit-for-bit what
    // it was.
    const std::vector<ChainInput> one = {input("/data/a.root", 'a', stagea("LD2"))};
    REQUIRE(chain_digest(one) == digest('a'));
}

TEST_CASE("the seed is the first 8 bytes of the digest, big-endian", "[inputchain][seed]") {
    // Pinned against a literal rather than against a re-derivation: this is the
    // number that reproduces every existing single-input result, and a test that
    // recomputes it the same way the code does would agree with any change.
    REQUIRE(seed_from_digest("644c6c4b39033dacae4a42e11306087f71b77ea1abeb2bf136dea50dcfb0b1da") ==
            0x644c6c4b39033dacULL);
    REQUIRE(seed_from_digest(std::string(64, '0')) == 0ULL);
    REQUIRE(seed_from_digest(std::string(64, 'f')) == 0xffffffffffffffffULL);
}

TEST_CASE("every distinct chain gets its own seed", "[inputchain][seed]") {
    // Two files, and either one alone, are three different datasets and must not
    // share a pool. In particular a chain's digest must not collapse to its first
    // input's -- which is what a fold that started from an empty accumulator, or
    // one that xor'd, could do.
    const std::vector<ChainInput> a = {input("/data/a.root", 'a', stagea("LD2"))};
    const std::vector<ChainInput> b = {input("/data/b.root", 'b', stagea("LD2"))};
    const std::vector<ChainInput> ab = {a[0], b[0]};

    REQUIRE(chain_digest(a) != chain_digest(b));
    REQUIRE(chain_digest(ab) != chain_digest(a));
    REQUIRE(chain_digest(ab) != chain_digest(b));
}

// ===========================================================================
// the budget
// ===========================================================================

TEST_CASE("--max-events is TOTAL across the chain, not per file", "[inputchain][budget]") {
    // THE DEFECT THIS PREVENTS: pass 1 truncates with RDataFrame::Range over the
    // whole chain -- N total, for free. Pass 0 counts by hand. A counter scoped
    // inside the per-file loop would give pass 0 N events PER FILE, so the pool
    // would be built from a strictly larger sample than the spectra it models,
    // with every printed number still plausible.
    //
    // This models the loop's use of budget_take: ONE budget, three files of four
    // events each. The reference-threading in main.cxx is what makes this shape
    // real, and that is checked at runtime by p0_events == n_events.
    long long budget = 5;
    std::vector<int> read_per_file;
    for (int file = 0; file < 3; ++file) {
        int n = 0;
        for (int entry = 0; entry < 4; ++entry) {
            if (!budget_take(budget)) break;
            ++n;
        }
        read_per_file.push_back(n);
    }

    REQUIRE(read_per_file == std::vector<int>{4, 1, 0});  // 5 TOTAL. Per-file would be {4, 4, 4}.
    REQUIRE(budget == 0);
}

TEST_CASE("a negative budget is unlimited and never decrements", "[inputchain][budget]") {
    // How "no --max-events" is spelled, matching make_grid's read_slim. It must
    // not wrap toward zero and start truncating on a long chain.
    long long budget = -1;
    for (int i = 0; i < 10000; ++i) REQUIRE(budget_take(budget));
    REQUIRE(budget == -1);
}

TEST_CASE("a spent budget stays spent", "[inputchain][budget]") {
    long long budget = 1;
    REQUIRE(budget_take(budget));
    REQUIRE(budget == 0);
    REQUIRE_FALSE(budget_take(budget));
    REQUIRE_FALSE(budget_take(budget));
    // It must not fall through zero into the negative, which is "unlimited".
    REQUIRE(budget == 0);
}

// ===========================================================================
// what the output provenance carries
// ===========================================================================

TEST_CASE("runs are carried forward as a LIST, in canonical order", "[inputchain][provenance]") {
    // Stage A pins ONE run per file by refusing a multi-run HIPO, because a header
    // naming one run "would be a lie". One output over N files means N runs, and
    // the honest record is all of them -- not the first one encountered.
    std::vector<ChainInput> chain = {input("/data/c.root", 'c', stagea("LD2", "22085")),
                                     input("/data/a.root", 'a', stagea("LD2", "22083")),
                                     input("/data/b.root", 'b', stagea("LD2", "22084"))};
    order_canonically(chain);
    REQUIRE(distinct_values(chain, "run") == std::vector<std::string>{"22083", "22084", "22085"});
}

TEST_CASE("a repeated run appears once", "[inputchain][provenance]") {
    std::vector<ChainInput> chain = {input("/data/a.root", 'a', stagea("LD2", "22083")),
                                     input("/data/b.root", 'b', stagea("LD2", "22083"))};
    order_canonically(chain);
    REQUIRE(distinct_values(chain, "run") == std::vector<std::string>{"22083"});
    REQUIRE(distinct_values(chain, "target") == std::vector<std::string>{"LD2"});
}

TEST_CASE("prov_lookup distinguishes an absent key from an empty value",
          "[inputchain][provenance]") {
    // The distinction Provenance::get() cannot make, and every check above turns
    // on it.
    Provenance p;
    p.add("empty", "");
    REQUIRE(prov_lookup(p, "empty").has_value());
    REQUIRE(prov_lookup(p, "empty")->empty());
    REQUIRE_FALSE(prov_lookup(p, "absent").has_value());
    REQUIRE(p.get("empty") == p.get("absent"));  // ... which is why get() will not do
}
