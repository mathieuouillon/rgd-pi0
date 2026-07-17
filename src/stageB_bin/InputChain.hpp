#pragma once

/// \file InputChain.hpp
/// \brief The N Stage A slim files of ONE Stage B run: their canonical order,
///        the cross-checks that must hold across them, and the seed they derive.
///        Header-only.
///
/// ---------------------------------------------------------------------------
/// WHY A STAGE B RUN TAKES MANY FILES AND WHY THAT IS DANGEROUS
/// ---------------------------------------------------------------------------
/// The donor pool is a reservoir per (Q^2, x_B, multiplicity) bin, and a bin the
/// input never populated yields NO mixed background at all -- not a thin one, a
/// missing one. One farm slim file fills 33 of the shipped 224 bins. The RG-D
/// production is ~30,866 slim files for LD2 alone, so the pool is only well
/// filled if ONE Stage B run sees ALL of a target's slims. Hence --input repeats.
///
/// That is also what makes this header necessary. Stage B reads a slim file's
/// provenance and propagates it; with one input, "the target" is whatever the
/// file said and there is nothing to get wrong. With N inputs there is exactly
/// one pool, and NOTHING in the reading code can tell an LD2 photon from an Sn
/// one -- a DonorPhoton is four floats. Chain an LD2 file with an Sn file and the
/// Sn photons become the mixed background subtracted from the LD2 spectra: the
/// numerator AND the denominator of R_A are corrupted, and the output is a set of
/// entirely ordinary-looking spectra. No plot would show it. No downstream check
/// could find it. So the checks live HERE, they run BEFORE the first read, and
/// they REFUSE rather than warn.
///
/// The rule this file follows throughout: a cross-input disagreement is refused,
/// never averaged, never last-wins, and never merely printed. A warning about a
/// pool that is already wrong is a warning nobody can act on after the fact.
///
/// ---------------------------------------------------------------------------
/// WHY THE ORDER IS THE CONTENT HASH AND NOT THE COMMAND LINE
/// ---------------------------------------------------------------------------
/// The pool is a reservoir sample, and reservoir sampling is a function of the
/// SEQUENCE of offers (DonorPool.hpp's BUILD phase; DonorPool.cpp calls
/// order-within-a-bin "the irreducible part"). So `--input a --input b` and
/// `--input b --input a` would build DIFFERENT pools from the SAME data. A farm
/// that globs a directory does not promise an order, and an analysis whose
/// background depends on the shell's glob order is the defect this whole rewrite
/// exists to end.
///
/// The inputs are therefore SORTED before anything reads them. The sort key is
/// each file's SHA-256 -- its CONTENT -- and NOT its path.
///
/// Sorting by path would be the obvious choice and it would be wrong, for the
/// exact reason stated at seed_from_file_content()'s old home in stageB_bin's
/// main.cxx: it would make the pool depend on where somebody put the files. Copy
/// a target's slims to a second mount, or rename them, and the sort order moves,
/// the offer sequence moves, and the same data gives a different background.
/// Ordering by content keeps the promise the seed already makes -- two identical
/// sets of files, anywhere, on any machine, under any names, give the same pool.
///
/// The sort key is the same digest the seed is derived from, so it costs nothing:
/// both need it and it is read once.
///
/// The resolved order is stamped into the output provenance regardless. A
/// deterministic order that is not recorded is only reproducible by someone who
/// re-derives it, and the point of the stamp is that the artefact says what it
/// did rather than that a reader trusts this comment.

#include <algorithm>
#include <cstdint>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "util/Provenance.hpp"
#include "util/Sha256.hpp"

namespace pi0::stageB {

/// One Stage A slim file of the chain, with everything read out of it that any
/// cross-check or the seed needs.
struct ChainInput {
    std::string path;
    std::string sha256;          ///< of the file's CONTENT. The sort key AND the seed's input.
    util::Provenance prov;       ///< Stage A's block, verbatim. Empty if the file carries none.
    long long n_events_used{};   ///< filled by the read pass; stamped per input, see below.
};

/// \return the value for `key`, or nullopt if the key is ABSENT.
///
/// Provenance::get() cannot make this distinction: it returns "" for an absent
/// key AND for a key whose value is "", and "" is a legal value. A validator
/// built on get() therefore cannot tell "this file does not say" from "this file
/// says nothing", and would silently pass a chain it never actually checked.
/// Every check below turns on that difference, so they search `entries` instead.
[[nodiscard]] inline std::optional<std::string> prov_lookup(const util::Provenance& p, const std::string& key) {
    for (const auto& [k, v] : p.entries) {
        if (k == key) return v;
    }
    return std::nullopt;
}

/// Consume one event from the chain-wide --max-events allowance.
///
/// \param budget  remaining allowance, decremented. NEGATIVE MEANS UNLIMITED.
/// \return false once the allowance is spent, i.e. "stop reading".
///
/// THE ALLOWANCE IS TOTAL ACROSS THE CHAIN, NOT PER FILE, and this function is
/// how that is kept true. It exists as a named thing rather than as two lines
/// open-coded in the read loop because the failure it prevents is invisible:
/// pass 1 truncates with RDataFrame::Range over the whole chain, which is "N
/// total" for free, while pass 0 counts by hand. Give pass 0 a counter scoped
/// inside the per-file loop and it truncates at N PER FILE -- so the pool is
/// built from a strictly larger sample than the spectra it has to model, on a
/// path where every number still looks plausible. The two implementations must
/// stay in lockstep; this is pass 0's half of that.
///
/// main.cxx re-checks the outcome (p0_events == n_events) after both passes.
/// That check is the cheapest test in the program and must not be weakened.
[[nodiscard]] inline bool budget_take(long long& budget) {
    if (budget == 0) return false;
    if (budget > 0) --budget;
    return true;
}

/// Put the chain into its canonical order: ascending by content hash.
///
/// \throws std::runtime_error if two inputs have the SAME content.
///
/// Duplicate content is refused rather than deduplicated or accepted. The same
/// slim file named twice -- a glob that matched a copy, a file list built by
/// concatenation -- is the same events twice: N_DIS doubles, every yield doubles,
/// and the ratio of the two looks entirely correct while the normalisation
/// underneath it is not. It cannot be a "warning" for the same reason the target
/// check cannot: by the time the output exists the damage is unreadable from it.
/// Refusing also makes the sort a TOTAL order with no tie-break to invent, which
/// is why this check lives here and not in the validator.
inline void order_canonically(std::vector<ChainInput>& inputs) {
    std::sort(inputs.begin(), inputs.end(),
              [](const ChainInput& a, const ChainInput& b) { return a.sha256 < b.sha256; });

    for (std::size_t i = 0; i + 1 < inputs.size(); ++i) {
        if (inputs[i].sha256 == inputs[i + 1].sha256) {
            throw std::runtime_error(
                "two inputs have identical content (sha256 " + inputs[i].sha256.substr(0, 16) + "...):\n  " +
                inputs[i].path + "\n  " + inputs[i + 1].path +
                "\nThe same events would be counted twice: N_DIS and every yield would both double, so the "
                "spectra would look right and the normalisation would not. Pass each slim file once.");
        }
    }
}

/// The fields that MUST agree across every input of one chain.
///
/// Each is here because one pool and one output provenance are built from all N
/// files, so a field that differs makes the output's single value a lie about
/// some of its inputs -- which is precisely what Stage A refuses to do when it
/// exits 5 on a HIPO spanning two runs ("the provenance header -- which records
/// a single run and a single model -- would be a lie"). Propagating one such
/// header across N chained files commits that same lie one stage later.
///
/// `run` is NOT in this list, deliberately: one slim file per run is Stage A's
/// contract, so chaining runs is the NORMAL case and the whole point of --input
/// repeating. Runs are carried forward as a LIST instead (see runs_in_chain).
[[nodiscard]] inline const std::vector<std::string>& chain_invariant_fields() {
    static const std::vector<std::string> f = {
        // THE ONE THAT CORRUPTS R_A SILENTLY. Nothing downstream of the read
        // knows a photon's target; the pool would mix Sn photons into LD2 events.
        "target",
        // Different cuts selected the photons that share one pool and one spectrum.
        "config.sha256",
        // Different GBT models = a different photon purity in one sample.
        "gbt.model",
        "gbt.fallback_used",
        // Both are stamped once into the output and are used by compute_sidis.
        "beam.energy_gev",
        "polarity",
    };
    return f;
}

/// Refuse anything about this chain that cannot be one Stage B run.
///
/// \param allow_truncated  honour --allow-truncated-inputs. See below.
/// \throws std::runtime_error naming the two disagreeing files and the field.
///
/// TRUNCATION. Stage A stamps `events.max_events_requested`, and on a --max-events
/// run it says the output is "a PREFIX of the input ... no yield or normalisation
/// from this file is complete". Nothing has ever read that flag. Chain one
/// truncated file with four full ones and N_DIS -- the normalisation denominator
/// -- is quietly wrong: the spectra are right, the events they are divided by are
/// not, and the result is a ratio that is wrong by a factor nobody can recover.
/// A truncated slim is not a production input, so it is refused. --allow-truncated
/// -inputs exists because a truncated slim IS the right input for a smoke test,
/// and it stamps its own loudness into the output rather than passing quietly.
///
/// A single input is exempt from the CROSS-file checks, and not as a shortcut:
/// there is no second file to disagree with, and this program's behaviour on one
/// --input must be exactly what it was before this file existed.
inline void validate_chain(const std::vector<ChainInput>& inputs, bool allow_truncated) {
    if (inputs.empty()) throw std::runtime_error("the chain is empty");

    // ---- truncation, per input, single or not ----------------------------
    for (const ChainInput& in : inputs) {
        const std::optional<std::string> mx = prov_lookup(in.prov, "events.max_events_requested");
        // Absent = a file with no provenance, or one written before the flag
        // existed. Unknowable, not false: it is left to the existing "no
        // provenance" warning rather than turned into a refusal here, which
        // would reject files this program has always accepted.
        if (!mx.has_value()) continue;
        if (mx->find("TRUNCATED") == std::string::npos) continue;
        if (allow_truncated) continue;
        throw std::runtime_error(
            "this input is a TRUNCATED Stage A run and is not a valid production input:\n  " + in.path +
            "\n  Stage A says events.max_events_requested = " + *mx +
            "\nIts event count is a PREFIX of its HIPO, so N_DIS -- the normalisation denominator -- would be "
            "short by an unknown amount while its spectra look complete. Chained with full files the shortfall "
            "is invisible. Re-skim it without --max-events, or pass --allow-truncated-inputs if this is a smoke "
            "test whose yields nobody will quote (the output will say so).");
    }

    if (inputs.size() == 1) return;

    // ---- the cross-file invariants ----------------------------------------
    // Everything is compared against inputs[0] rather than pairwise: the fields
    // must all be EQUAL, so equality with one representative is equality with
    // each other, and the message names a pair a human can actually diff.
    const ChainInput& ref = inputs.front();
    for (const std::string& field : chain_invariant_fields()) {
        const std::optional<std::string> a = prov_lookup(ref.prov, field);
        for (std::size_t i = 1; i < inputs.size(); ++i) {
            const std::optional<std::string> b = prov_lookup(inputs[i].prov, field);
            if (a.has_value() && b.has_value() && *a == *b) continue;
            if (!a.has_value() && !b.has_value()) continue;  // neither file says; see below

            std::ostringstream os;
            os << "the inputs disagree on \"" << field << "\", so they cannot be one Stage B run:\n  "
               << ref.path << "\n    " << field << " = " << (a ? *a : std::string("<absent>")) << "\n  "
               << inputs[i].path << "\n    " << field << " = " << (b ? *b : std::string("<absent>"));
            if (field == "target") {
                os << "\nONE donor pool is built from every input, and nothing downstream of the read knows a "
                      "photon's target -- a donor is four floats. Chaining these would make one target's photons "
                      "the mixed background subtracted from the other's events, corrupting both the numerator and "
                      "the denominator of R_A, and the spectra would look entirely ordinary.";
            } else {
                os << "\nThis run would stamp ONE value for that field into an output built from both files, "
                      "which would make the output's provenance false for at least one of them.";
            }
            os << "\nRun Stage B once per " << (field == "target" ? "target" : "value") << ".";
            throw std::runtime_error(os.str());
        }
    }

    // A chain in which NOBODY names a target is not a chain whose targets were
    // checked -- the loop above passes it by agreeing that nothing is known. For
    // one input that is the old, warned-about "no provenance" case and is
    // allowed; for many it is the R_A corruption above with the evidence missing.
    if (!prov_lookup(ref.prov, "target").has_value()) {
        throw std::runtime_error(
            "no input names a target (no Stage A `target` in its provenance), so this chain's targets cannot be "
            "shown to match. One pool is built from all " +
            std::to_string(inputs.size()) +
            " files and nothing downstream can separate their photons. Re-skim with a Stage A that stamps "
            "provenance, or run one input at a time.");
    }
}

/// The distinct values of `key` across the chain, in canonical (chain) order.
///
/// For `run` this is the union that the output provenance carries INSTEAD of a
/// single run. Stage A pins one run per file by refusing a multi-run HIPO,
/// precisely because a header naming one run would be a lie about the rest; a
/// chain of 30,866 files spans 30,866 runs, and the honest record of that is a
/// list, not the first one encountered.
[[nodiscard]] inline std::vector<std::string> distinct_values(const std::vector<ChainInput>& inputs,
                                                              const std::string& key) {
    std::vector<std::string> out;
    for (const ChainInput& in : inputs) {
        const std::optional<std::string> v = prov_lookup(in.prov, key);
        if (!v.has_value() || v->empty()) continue;
        if (std::find(out.begin(), out.end(), *v) == out.end()) out.push_back(*v);
    }
    return out;
}

/// The chain's content digest: what the seed is derived from.
///
/// \param inputs  ALREADY in canonical order (order_canonically).
///
/// ONE INPUT'S CHAIN DIGEST IS THAT INPUT'S DIGEST, not a hash of it. That is
/// not a special case for its own sake: it is what keeps a single-input run
/// bit-for-bit identical to every single-input run made before --input could
/// repeat. Re-hashing a lone digest would silently move the seed, hence the
/// pool, hence the mixed background of every existing result, for no gain.
///
/// For many inputs the digest is the SHA-256 of the concatenated per-input
/// digests in canonical order. Since that order is itself derived from the
/// content, the seed does not depend on the order the files were typed on the
/// command line -- `--input a --input b` and `--input b --input a` give the same
/// seed and the same pool, which is the property that makes a farm glob safe.
///
/// The digests are concatenated rather than xor'd or summed: those commute, so
/// they would collide on any permutation of the same files, which is exactly the
/// case this has to distinguish from -- a chain is a SEQUENCE of offers.
[[nodiscard]] inline std::string chain_digest(const std::vector<ChainInput>& inputs) {
    if (inputs.empty()) throw std::runtime_error("the chain is empty; there is no digest");
    if (inputs.size() == 1) return inputs.front().sha256;

    util::Sha256 h;
    for (const ChainInput& in : inputs) {
        h.update(reinterpret_cast<const std::uint8_t*>(in.sha256.data()), in.sha256.size());
    }
    return h.final_hex();
}

/// The donor pool's seed: the first 8 bytes of `digest`, big-endian.
///
/// Split out of the old seed_from_file_content() so that "which bytes become the
/// seed" is one function with one answer whether the chain holds one file or
/// thirty thousand. cuts.json's mixing.seed_mode is "file_hash" and Cuts::load
/// refuses any other mode; this, plus chain_digest above, is what "file_hash"
/// means. Only 64 of the 256 bits are used, which is fine: this seeds a
/// std::mt19937_64, and the question is "do two different chains collide by
/// accident", not "can someone forge a seed".
[[nodiscard]] inline std::uint64_t seed_from_digest(const std::string& digest) {
    if (digest.size() < 16) throw std::runtime_error("a sha256 digest is 64 hex characters; got " + digest);
    std::uint64_t seed = 0;
    for (int i = 0; i < 16; ++i) {  // 16 hex chars = 8 bytes
        const char c = digest[static_cast<std::size_t>(i)];
        const std::uint64_t nibble = (c >= '0' && c <= '9') ? static_cast<std::uint64_t>(c - '0')
                                                            : static_cast<std::uint64_t>(c - 'a' + 10);
        seed = (seed << 4) | nibble;
    }
    return seed;
}

}  // namespace pi0::stageB
