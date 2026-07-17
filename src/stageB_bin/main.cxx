// ---------------------------------------------------------------------------
// stageB_bin -- Stage A slim TTree -> binned m_gg spectra and count-weighted sums.
//
// usage:
//   stageB_bin --input <slim.root> [--input <more.root> ...]
//              --output <binned.root> --config <cuts.json>
//              --grid-a <grid_A_q2_xb.json> --grid-b <grid_B_z_pt2.json>
//              [--max-events N] [--allow-truncated-inputs]
//
// ===========================================================================
// ONE RUN, MANY INPUTS, ONE POOL -- AND ONE TARGET
// ===========================================================================
// --input REPEATS, and the N files are ONE run: one donor pool drawn from all of
// them, one set of spectra, one output. That is not a convenience. The pool is a
// reservoir per (Q^2, x_B, multiplicity) bin, and a bin no input populated has NO
// mixed background rather than a thin one -- one farm slim file fills 33 of the
// shipped 224 bins, while the RG-D production is ~30,866 slim files for LD2
// alone. The pool is only well filled if one run sees all of a target's slims.
//
// EVERYTHING THAT MAKES THAT SAFE IS IN src/stageB_bin/InputChain.hpp, and it
// runs BEFORE the first event is read: the chain is put into a canonical order
// (by CONTENT hash, not by path -- see there), it is refused if the inputs
// disagree about their target, their cuts, their GBT model, their beam energy or
// their polarity, and it is refused if any input is a truncated Stage A run.
// Those are refusals and not warnings for one reason: nothing downstream of the
// read knows a photon's target -- a donor is four floats -- so an LD2 file
// chained with an Sn file builds one pool in which Sn photons are the mixed
// background subtracted from LD2 events. It corrupts the numerator AND the
// denominator of R_A and it produces entirely ordinary-looking spectra. There is
// no plot on which it shows up and no later stage that could find it.
//
// ===========================================================================
// WHAT THIS PROGRAM IS FOR. READ THIS BEFORE CHANGING ANY ACCUMULATOR.
// ===========================================================================
// The superseded analysis reported every bin at the GEOMETRIC MIDPOINT of its
// box (note sec:binning-caveat). For an interior box that is harmless. For the
// OUTERMOST box in each dimension -- which runs to the kinematic limit and is far
// wider than its neighbours -- it is nowhere near the data. Measured from the
// real result files:
//
//   * 52.6% of bins sat in the top box of at least one dimension;
//   * the "z ~ 0.70" point was really an integral over z in [0.37, 1.0],
//     dominated by z ~ 0.45;
//   * 15.5% of bins implied y > 0.85 -- violating the DIS cut their own events
//     had passed;
//   * 6.9% implied nu > E_beam. KINEMATICALLY IMPOSSIBLE. Up to nu = 21.7 GeV
//     against a 10.53 GeV beam.
//
// It made the measured attenuation look ~3x weaker than it is, and it could not
// be repaired after the fact, because the per-bin means were never accumulated.
// The histograms alone cannot reconstruct them. That is the whole reason this
// program exists, and the equal-statistics grid does NOT fix it: an
// equal-statistics rectangular grid has exactly the same unbounded outermost bin
// as the kd-tree did. Measured on the 799-event smoke file by make_grid: the top
// z bin is 9.1x the median interior width and its midpoint (0.658) misses the
// true count-weighted mean (0.486) by 1.35x; top p_T^2 misses by 1.76x.
//
// THEREFORE, and this is the single non-negotiable thing in this file:
//
//     per (4D bin, m_gg bin) this program accumulates
//         counts, sum_q2, sum_xb, sum_z, sum_pt2
//
// PER (4D bin, m_gg bin) -- not per 4D bin. The yield is extracted from a
// BACKGROUND-SUBTRACTED spectrum, so the mean that belongs on the abscissa is the
// mean over exactly the events under the fitted peak. Accumulating per m_gg bin
// lets the downstream sideband-subtract the SUMS the same way it subtracts the
// yield. A mean over the whole (deliberately wide) mass window would be
// contaminated by the combinatorial background, which populates the bin
// differently from the signal -- so it would be a mean over the wrong events, and
// it would look perfectly reasonable.
//
// The same argument applies twice more, and both are implemented here:
//   * pT broadening: counts, sum_pt2, sum_pt4 per (3D bin, m_gg bin). The old
//     analysis accumulated these with NO background subtraction over a +-200 MeV
//     window, which is why its pT broadening is diluted by an unknown factor.
//   * the BSA: counts per (4D bin, m_gg bin, phi_h bin, helicity). Binning the
//     BSA in m_gg is what lets the dilution S/(S+B) be MEASURED rather than
//     ignored.
//
// core/Binning.hpp deliberately offers no bin_centre() and no box-midpoint
// accessor. There is nothing to call because nothing should call it. This file
// is where the thing that should be called instead gets filled.
//
// ===========================================================================
// TWO PASSES
// ===========================================================================
// PASS 0 -- build the frozen donor pool. SEQUENTIAL, by DonorPool's contract:
//   reservoir sampling is a function of the SEQUENCE of offers, so a parallel
//   pre-pass would make the pool depend on thread scheduling again, which is the
//   exact defect DonorPool exists to remove.
// PASS 1 -- RDataFrame over the same files, in the same canonical order. The pool
//   is frozen and const, so mixing is a pure per-entry lookup with no state, no
//   locks and no ordering.
//
// Both passes read the same PREFIX of the CHAIN under --max-events -- N events in
// total across every input, not N per file. They must: a pool built from more
// events than were binned would be a background estimated from a sample the
// spectra never saw. The two implementations of that truncation are independent
// (pass 0 counts by hand through InputChain's budget_take; pass 1 hands the whole
// chain to RDataFrame::Range, which is already chain-wide) and the p0_events ==
// n_events check below is what keeps them in lockstep.
//
// WHY PASS 1 IS NOT MULTITHREADED, THOUGH THE POOL WAS DESIGNED TO ALLOW IT
// ------------------------------------------------------------------------
// The pool is frozen, const and lock-free precisely so that this pass CAN be
// parallel, and that door stays open. It is not opened here, because the
// accumulators are sums of DOUBLES: with per-slot partials merged at the end, the
// summation order depends on how the threads happened to divide the file, so
// sum_q2 would differ in its last bits between a 4-thread and an 8-thread run.
// Bit-reproducibility is exactly what this rewrite is about -- an analysis whose
// binning changed between two passes over the same data is what it replaces --
// and buying speed with a non-deterministic abscissa would be trading away the
// one thing this program is for. Enabling MT needs a deterministic reduction
// (fixed per-slot partition summed in slot order, or fixed-point accumulation),
// which is a real change and is not in scope. The counts, being integers, would
// be fine; the sums are the problem.
//
// ===========================================================================
// THE MIXING ASYMMETRY, AND THE INVARIANT THAT MAKES IT CORRECT
// ===========================================================================
// A mixed pair is one photon from THIS event and one from the frozen pool. The
// e-gamma cut (pairing.e_gamma_min_angle_deg) is applied ONLY to the current
// event's photon. That is not a corner cut -- it is correct by construction: the
// donor carries no meaningful angle to THIS event's electron, and it was already
// tested against its OWN event's electron, which is the physically meaningful
// reference. DonorPool.hpp says the same thing and says not to "fix" it.
//
// What made the old code wrong was a different thing entirely: it inserted ALL of
// the current event's photons into the pool, INCLUDING the ones that had FAILED
// their own event's e-gamma cut. So a photon rejected from the same-event
// spectrum still contributed to the mixed one, distorting the background shape
// near the electron direction -- silently, because the result is an
// ordinary-looking spectrum.
//
//     *** THE INVARIANT: this program offers to the pool ONLY photons that
//     *** passed their own event's e-gamma cut. It is enforced at exactly one
//     *** place -- build_photons() below -- which is the only thing in this file
//     *** that turns a slim-file entry into a photon list, and both passes call
//     *** it. There is no path on which a failing photon reaches offer().
//
// DonorPool cannot check this for you: from inside offer(), a photon that failed
// the cut is indistinguishable from one that passed. Its multiplicity tripwire
// catches a caller who filtered the photons but not the count; it catches nothing
// if you get both wrong the same way.
//
// EVERYTHING ELSE ABOUT A MIXED PAIR IS IDENTICAL TO A SAME-EVENT PAIR: the same
// admissibility rule (pi0::admissible_pair -- the SAME function find_gg_pairs
// calls, so there is one copy of the cut and the two spectra cannot drift apart),
// the same z window, the same grids, the same m_gg axis. The one structural
// difference is the pairing SHAPE: find_gg_pairs is greedy and exclusive over one
// list, while mixing is exhaustive and non-exclusive across two -- see the
// mix_event() comment.
//
// ===========================================================================
// HELICITY
// ===========================================================================
// helicity == 0 means UNDEFINED (REC::Event.helicity, "0=UDF"), not a third
// state. It is dropped from the BSA histogram and from NOTHING ELSE. Those events
// still enter the m_gg spectra, the kinematic sums, N_DIS and the donor pool,
// because every one of those is helicity-blind. Dropping them from the pool would
// shrink the donor sample for no physics reason and -- worse -- would give the
// mixed background a different parent sample from the same-event spectrum it has
// to subtract from. A mixed pair has two events' helicities and therefore none;
// the mixed spectrum is not helicity-resolved and cannot be.
// ---------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include "config/Cuts.hpp"
#include "core/Binning.hpp"
#include "core/Constants.hpp"
#include "core/Kinematics.hpp"
#include "core/Pairing.hpp"
#include "core/Types.hpp"
#include "stageB_bin/DonorPool.hpp"
#include "stageB_bin/InputChain.hpp"
#include "util/Provenance.hpp"
#include "util/Sha256.hpp"

// Set from meson.project_version(). See the provenance block in main().
#ifndef PI0_CODE_VERSION
#error "PI0_CODE_VERSION is not defined. It must be injected by meson (see src/stageB_bin/meson.build)."
#endif

namespace {

using RVecD = ROOT::VecOps::RVec<double>;
using pi0::stageB::ChainInput;
using pi0::util::Provenance;
using pi0::util::utc_now;

// ===========================================================================
// CLI
// ===========================================================================

struct Args {
    /// Every --input, in the order given. REORDERED before use: the pool depends
    /// on the sequence of offers, so the command line's order must not reach it.
    /// See InputChain.hpp.
    std::vector<std::string> inputs;
    std::string output;
    std::string config;
    std::string grid_a;
    std::string grid_b;
    /// Stop after this many entries of the CHAIN -- N in total across every
    /// input, not N per input. nullopt = all of them. `unsigned int` because that
    /// is what RDataFrame::Range takes; the budget is a long long everywhere else.
    std::optional<unsigned int> max_events;
    /// Accept a Stage A file that says it is a truncated run. See InputChain.hpp.
    bool allow_truncated_inputs{false};
};

void print_usage(const char* argv0) {
    std::cerr << "usage: " << argv0 << " --input <slim.root> [--input <more.root> ...]\n"
              << "                       --output <binned.root>\n"
              << "                       --config <cuts.json>\n"
              << "                       --grid-a <grid_A_q2_xb.json> --grid-b <grid_B_z_pt2.json>\n"
              << "                       [--max-events <N>] [--allow-truncated-inputs]\n"
              << "\n"
              << "  --input       a Stage A slim ROOT file (TTree \"events\"). Repeat for more; the N\n"
              << "                files are ONE run over ONE donor pool, which is the point -- a pool\n"
              << "                bin no input filled has NO mixed background, and one farm file fills\n"
              << "                33 of the 224 shipped bins. Give a target ALL of its slims.\n"
              << "                The order does not matter: the chain is sorted by content hash before\n"
              << "                anything reads it, so a glob is safe and the pool is the same either\n"
              << "                way. The resolved order is stamped into the output.\n"
              << "                THE INPUTS MUST BE ONE TARGET. A chain disagreeing about its target,\n"
              << "                cuts, GBT model, beam energy or polarity is REFUSED, not warned about:\n"
              << "                nothing downstream knows a photon's target, so a mixed chain would\n"
              << "                subtract one target's background from another's events and say nothing.\n"
              << "  --output      binned ROOT file to create (overwritten if it exists).\n"
              << "  --config      the cut configuration. Every threshold and every axis comes from here.\n"
              << "  --grid-a      the frozen (Q^2, x_B) grid. Its hash is stamped into the output.\n"
              << "  --grid-b      the frozen (z, p_T^2) grid. Its hash is stamped into the output.\n"
              << "  --max-events  stop after the first N entries of the CHAIN (default: all). N TOTAL\n"
              << "                ACROSS ALL INPUTS, not N per input: the inputs are read in canonical\n"
              << "                order until the allowance is spent, so later ones may go unread. BOTH\n"
              << "                passes are truncated identically -- a pool built from more events than\n"
              << "                were binned would estimate the background from a sample the spectra\n"
              << "                never saw. The output is a PREFIX of the chain, not a sample of it;\n"
              << "                the truncation is stamped into the provenance so a partial run cannot\n"
              << "                be mistaken for a full one.\n"
              << "  --allow-truncated-inputs\n"
              << "                accept an input whose own Stage A provenance says it is a TRUNCATED\n"
              << "                run. Refused by default: such a file's event count is a prefix of its\n"
              << "                HIPO, so N_DIS -- the normalisation denominator -- is short by an\n"
              << "                unknown amount while its spectra look complete, and chained with full\n"
              << "                files the shortfall is invisible. For smoke tests whose yields nobody\n"
              << "                will quote. Stamped loudly into the output.\n";
}

/// \throws std::runtime_error on any malformed or missing argument.
[[nodiscard]] Args parse_args(int argc, char* argv[]) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        const std::string flag = argv[i];
        auto value = [&]() -> std::string {
            if (i + 1 >= argc) throw std::runtime_error(flag + " requires a value");
            return argv[++i];
        };
        if (flag == "--input") a.inputs.push_back(value());
        else if (flag == "--allow-truncated-inputs") a.allow_truncated_inputs = true;
        else if (flag == "--output") a.output = value();
        else if (flag == "--config") a.config = value();
        else if (flag == "--grid-a") a.grid_a = value();
        else if (flag == "--grid-b") a.grid_b = value();
        else if (flag == "--max-events") {
            const std::string v = value();
            unsigned long long n = 0;
            try {
                // The '-' check is not decoration: stoull turns "-1" into
                // ULLONG_MAX, which would clamp to "all events" -- a typo that
                // silently means the opposite of what it says. Same guard, and
                // the same reasoning, as stageA_skim.
                if (v.find('-') != std::string::npos) throw std::invalid_argument("negative");
                std::size_t pos = 0;
                n = std::stoull(v, &pos);
                if (pos != v.size()) throw std::invalid_argument("trailing characters");
            } catch (const std::exception&) {
                throw std::runtime_error("--max-events wants a positive integer, got \"" + v + "\"");
            }
            if (n == 0) {
                throw std::runtime_error(
                    "--max-events 0 would read nothing and write an empty file. Omit the flag to process "
                    "the whole input.");
            }
            if (n > std::numeric_limits<unsigned int>::max()) {
                // The cap is RDataFrame::Range's, not this program's, and it is
                // reachable now that one run chains a whole target: 30,866 slim
                // files pass 2^32 events long before they run out. It is a
                // refusal rather than a clamp for the reason the '-' check above
                // exists -- a silently ignored limit reads as a full run.
                throw std::runtime_error("--max-events " + v + " exceeds the largest range this program can take (" +
                                         std::to_string(std::numeric_limits<unsigned int>::max()) +
                                         "). Omit the flag to process every event of every input.");
            }
            a.max_events = static_cast<unsigned int>(n);
        } else {
            throw std::runtime_error("unknown argument \"" + flag + "\"");
        }
    }
    if (a.inputs.empty()) throw std::runtime_error("at least one --input is required");
    if (a.output.empty()) throw std::runtime_error("--output is required");
    if (a.config.empty()) throw std::runtime_error("--config is required");
    if (a.grid_a.empty()) throw std::runtime_error("--grid-a is required");
    if (a.grid_b.empty()) throw std::runtime_error("--grid-b is required");
    return a;
}

[[nodiscard]] std::string fmt(double v) {
    std::ostringstream os;
    os << v;
    return os.str();
}

// ===========================================================================
// axes
// ===========================================================================

/// Build a uniform Grid1D of `n` bins over [lo, hi] from a config block.
///
/// A Grid1D rather than the obvious `int((v - lo) / w * n)` arithmetic, and this
/// is deliberate. Grid1D is the project's ONE axis type, and its semantics are
/// written down and pinned by tests: half-open [lo, hi), a value exactly on the
/// TOP edge lands in the LAST bin rather than nowhere (the old find_1d_bin()
/// returned -1 there and so dropped every p_T^2 equal to the top edge),
/// out-of-range and NaN give -1, no flow bins. Open-coding the arithmetic here
/// would give the m_gg and phi_h axes a SECOND set of edge conventions, disagreeing
/// with the analysis grids about which bin an edge value belongs to -- and that
/// disagreement would be invisible in every histogram it corrupted.
///
/// The cost is a binary search over ~200 edges instead of a multiply, on a path
/// taken once per pair. It is not measurable next to the four-vector algebra that
/// produced the value being binned.
///
/// The edges are computed as lo + i * (hi - lo) / n rather than by repeated
/// addition, so rounding cannot accumulate and the last edge is exactly `hi`.
[[nodiscard]] pi0::Grid1D uniform_axis(std::string name, double lo, double hi, int n) {
    if (n < 1) throw std::runtime_error("uniform_axis(" + name + "): needs at least one bin, got " + std::to_string(n));
    if (!(lo < hi)) throw std::runtime_error("uniform_axis(" + name + "): lo must be below hi");

    pi0::Grid1D g;
    g.name = std::move(name);
    g.edges.resize(static_cast<std::size_t>(n) + 1);
    for (int i = 0; i <= n; ++i) {
        g.edges[static_cast<std::size_t>(i)] = lo + (hi - lo) * static_cast<double>(i) / static_cast<double>(n);
    }
    g.edges.back() = hi;  // exact, not lo + n*(hi-lo)/n
    return g;
}

// ===========================================================================
// the seed
// ===========================================================================
//
// THE CONTENT, NOT THE PATH -- and now the CHAIN's content, not one file's.
// pi0::stageB::chain_digest() + seed_from_digest() are what cuts.json's
// mixing.seed_mode = "file_hash" means; both live in InputChain.hpp, next to the
// canonical ordering they share a digest with, because the ordering and the seed
// have to make the same promise or neither does: the same data, anywhere, under
// any names, in any command-line order, gives the same pool.
//
// Hashing the path would be cheaper and WRONG in a specific, nasty way -- copying
// a file to another mount would change the seed and therefore the mixed
// background, so the same data would give two answers from two directories. The
// same argument is why the chain is sorted by content and not by path.
//
// Cost: one streaming read of each input. Real -- tens of seconds for a multi-GB
// slim -- but this program then reads every input TWICE more, so it is a fraction
// of the runtime and it buys the property the whole design rests on.

/// One input's photons, e-gamma cut APPLIED, as DonorPhotons.
///
/// *** PASS 0's HALF OF THE e-gamma INVARIANT, IN ONE PLACE. *** This is the only
/// thing on the pass 0 path that turns a slim entry into photons, so it is called
/// once per file rather than copy-pasted per input -- N copies of this filter
/// would be N chances for the pool to be drawn from a different photon population
/// than the spectra it models. It already had to agree with build_photons(); a
/// per-input loop must not multiply that.
[[nodiscard]] std::vector<pi0::DonorPhoton> build_donor_photons(const TTreeReaderArray<double>& gpx,
                                                                const TTreeReaderArray<double>& gpy,
                                                                const TTreeReaderArray<double>& gpz,
                                                                const TTreeReaderArray<double>& gang,
                                                                const pi0::Cuts& cuts) {
    std::vector<pi0::DonorPhoton> out;
    const std::size_t n = std::min({gpx.GetSize(), gpy.GetSize(), gpz.GetSize(), gang.GetSize()});
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (!(gang[i] > cuts.pairing.e_gamma_min_angle_deg)) continue;
        const double px = gpx[i], py = gpy[i], pz = gpz[i];
        const double p = std::sqrt(px * px + py * py + pz * pz);
        pi0::DonorPhoton d{};
        d.px = static_cast<float>(px);
        d.py = static_cast<float>(py);
        d.pz = static_cast<float>(pz);
        // Photons are massless: E = |p|, exactly as core/Types.hpp and Stage A
        // define them.
        d.e = static_cast<float>(p);
        d.inv_p = static_cast<float>(p > 0.0 ? 1.0 / p : 0.0);
        out.push_back(d);
    }
    return out;
}

/// What pass 0 counted, over the whole chain.
struct PoolPassStats {
    long long events{};
    long long offered_events{};
    long long no_class{};
};

/// Offer ONE input's events to the pool, sequentially.
///
/// \param budget  the CHAIN-WIDE --max-events allowance, decremented across
///                inputs. Negative = unlimited. Threaded by reference, exactly as
///                make_grid's read_slim does, because the alternative -- a counter
///                scoped to this function -- would truncate at N PER FILE while
///                pass 1's RDataFrame::Range truncates at N over the chain, and
///                the pool would then model a strictly different sample from the
///                spectra with every printed number still plausible.
/// \param n_events_used  this input's own contribution, for the provenance. Under
///                --max-events a file may be reached with the allowance already
///                spent and contribute 0; stamping it per input is what makes
///                "which files actually fed this pool" auditable rather than
///                inferable from the flag.
void read_pool_pass(const std::string& path, const pi0::Cuts& cuts, pi0::DonorPool& pool, PoolPassStats& stats,
                    long long& budget, long long& n_events_used) {
    std::unique_ptr<TFile> in(TFile::Open(path.c_str(), "READ"));
    if (in == nullptr || in->IsZombie()) throw std::runtime_error("cannot open the input file: " + path);

    TTreeReader reader("events", in.get());
    if (reader.GetTree() == nullptr) {
        throw std::runtime_error("the input has no TTree named \"events\": " + path);
    }
    TTreeReaderValue<double> r_q2(reader, "q2");
    TTreeReaderValue<double> r_xb(reader, "xb");
    TTreeReaderArray<double> r_gpx(reader, "gpx");
    TTreeReaderArray<double> r_gpy(reader, "gpy");
    TTreeReaderArray<double> r_gpz(reader, "gpz");
    TTreeReaderArray<double> r_gang(reader, "g_e_gamma_deg");

    bool stopped_on_budget = false;
    while (reader.Next()) {
        if (!pi0::stageB::budget_take(budget)) {
            stopped_on_budget = true;
            break;
        }
        ++stats.events;
        ++n_events_used;

        // Same filter as pass 1, same cut, one function. If this ever diverges
        // from build_photons(), the pool is drawn from a different photon
        // population than the spectra it has to model.
        const std::vector<pi0::DonorPhoton> passing = build_donor_photons(r_gpx, r_gpy, r_gpz, r_gang, cuts);

        const int b = pool.pool_bin(*r_q2, *r_xb, passing.size());
        if (b < 0) {
            // Off the pool grid in Q^2 or x_B, or in no multiplicity class --
            // which INCLUDES zero passing photons. Not an error: an event with no
            // photon has nothing to donate.
            ++stats.no_class;
            continue;
        }
        pool.offer(b, passing);
        ++stats.offered_events;
    }

    // Ported from make_grid, and chaining is what makes it worth having: a
    // schema-divergent file is now one file in N rather than the only file, so a
    // missing branch must surface as an error (kEntryDictionaryError /
    // kEntryChainSetupError) rather than as a column of silent zeroes in the pool.
    //
    // ONLY WHEN THE LOOP RAN TO THE END OF THE TREE. The status describes the
    // read that ended the loop, so it is evidence only when the READER ended it.
    // Break on the budget and the last Next() DID succeed -- the status is
    // kEntryValid and we simply chose not to use the entry -- so testing it
    // unconditionally throws on a healthy file the moment --max-events runs out
    // mid-file. (make_grid's copy had exactly that defect, on the stated
    // assumption that a clean break "leaves kEntryNotFound behind on some ROOT
    // versions"; it does not on this one. Fixed there too, the same way.)
    //
    // kEntryBeyondEnd is the ONLY clean end: per TTreeReader.h it is "last entry
    // loop has reached its end", while kEntryNotFound is "the tree entry number
    // does not exist" -- a real failure sitting alongside kEntryChainSetupError
    // and kEntryDictionaryError. It was only ever excused here as the guess about
    // the budget break; `stopped_on_budget` states that directly, so excusing it
    // as well would mask the error it names.
    if (!stopped_on_budget && reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd) {
        throw std::runtime_error("failed reading \"events\" from " + path + " (TTreeReader status " +
                                 std::to_string(static_cast<int>(reader.GetEntryStatus())) +
                                 "). Check the branch list matches the Stage A schema.");
    }
}

// ===========================================================================
// photons
// ===========================================================================

/// Turn one slim-file entry's photon columns into the event's photon list,
/// APPLYING THE e-gamma CUT.
///
/// *** THE ONE PLACE THE e-gamma CUT IS APPLIED. *** Both passes call this and
/// nothing else builds a photon list, so "the pool only ever received photons
/// that passed their own event's e-gamma cut" is a property of the call graph
/// rather than of a comment. See this file's header for why that invariant is the
/// difference between this and the old code.
///
/// Stage A deliberately does NOT apply this cut -- it ships g_e_gamma_deg per
/// photon and leaves the cut to the consumer -- so this is where it happens, and
/// `n_photons` everywhere downstream (the pool's multiplicity class included)
/// means "photons that passed it".
///
/// The four columns are RVecs of the same length by construction; std::min over
/// their sizes rather than trusting that, because a truncated file should drop
/// photons rather than read past an array.
[[nodiscard]] std::vector<pi0::Photon> build_photons(const RVecD& gpx, const RVecD& gpy, const RVecD& gpz,
                                                     const RVecD& g_e_gamma_deg, const pi0::Cuts& cuts) {
    std::vector<pi0::Photon> out;
    const std::size_t n = std::min({gpx.size(), gpy.size(), gpz.size(), g_e_gamma_deg.size()});
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (!(g_e_gamma_deg[i] > cuts.pairing.e_gamma_min_angle_deg)) continue;
        out.push_back(pi0::Photon{gpx[i], gpy[i], gpz[i], g_e_gamma_deg[i]});
    }
    return out;
}

// ===========================================================================
// accumulators
// ===========================================================================

/// One (4D bin, m_gg bin) cell.
///
/// n_same IS the count that weights the sums. They are the same number and they
/// live in the same struct so that they cannot be made to disagree: the abscissa
/// is sum_x / n_same over whatever m_gg range the downstream integrates, and the
/// yield in that same range is a sum of the same n_same. A separate "counts"
/// column beside the spectrum would be two names for one quantity and an
/// invitation for one of them to be filled on a slightly different condition --
/// which is how the old analysis ended up with a Q2-vs-xB histogram that looked
/// like a normalisation and was not one.
struct SpectrumCell {
    std::int64_t n_same{};   ///< same-event pairs. The abscissa's weight AND the yield's raw count.
    std::int64_t n_mixed{};  ///< mixed-event pairs. NOT weighted into any sum below.
    // SAME-EVENT ONLY. A mixed pair is not a pi0 -- it is one photon from this
    // event and one from an unrelated one -- so its kinematics describe nothing
    // and must never reach the abscissa. This is why the sums sit next to n_same
    // and not next to n_mixed.
    double sum_q2{};
    double sum_xb{};
    double sum_z{};
    double sum_pt2{};
};

/// One (3D bin, m_gg bin) cell, for pT broadening. SAME-EVENT ONLY, for the same
/// reason.
///
/// sum_pt4 is here so that the VARIANCE of p_T^2 -- and hence the uncertainty on
/// Delta<p_T^2> -- can be formed from sideband-subtracted sums, exactly like the
/// mean. Accumulating only sum_pt2 would give a subtractable mean with an
/// unsubtractable error bar.
struct Ptb3dCell {
    std::int64_t counts{};
    double sum_pt2{};
    double sum_pt4{};
};

}  // namespace

int main(int argc, char* argv[]) {
    // ---- arguments -------------------------------------------------------
    Args args;
    try {
        args = parse_args(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n\n";
        print_usage(argv[0]);
        return 1;
    }

    try {
        // ---- configuration ------------------------------------------------
        const pi0::Cuts cuts = pi0::Cuts::load(args.config);
        const pi0::Binning binning = pi0::Binning::load(args.grid_a, args.grid_b);
        const std::string config_sha = pi0::util::sha256_file(args.config);

        const pi0::Grid1D mgg_axis =
            uniform_axis("mgg", cuts.mgg_histogram.min_gev, cuts.mgg_histogram.max_gev, cuts.mgg_histogram.bins);
        const pi0::Grid1D phi_axis =
            uniform_axis("phi_h", cuts.bsa.phi_min_deg, cuts.bsa.phi_max_deg, cuts.bsa.n_phi_bins);

        const int n_mgg = mgg_axis.nbins();
        const int n_phi = phi_axis.nbins();
        const int n4d = binning.n4d();
        const int n3d = binning.n3d();
        const int n_a = binning.A.ncells();

        // ---- the chain -----------------------------------------------------
        // Read every input's identity and provenance, order the chain, and
        // REFUSE anything that cannot be one run -- ALL of it before a single
        // event is read. Two full passes over 30,866 files is not a thing to
        // discover a mixed target at the end of, and a chain that should not have
        // run must not leave a plausible output behind it.
        std::vector<ChainInput> chain;
        chain.reserve(args.inputs.size());
        for (const std::string& path : args.inputs) {
            ChainInput ci;
            ci.path = path;
            ci.sha256 = pi0::util::sha256_file(path);
            {
                std::unique_ptr<TFile> in(TFile::Open(path.c_str(), "READ"));
                if (in == nullptr || in->IsZombie()) throw std::runtime_error("cannot open the input file: " + path);
                ci.prov = Provenance::read(*in);
            }
            chain.push_back(std::move(ci));
        }
        pi0::stageB::order_canonically(chain);
        pi0::stageB::validate_chain(chain, args.allow_truncated_inputs);

        // The paths in canonical order. Pass 1 hands this to RDataFrame, which
        // reads a chain in the order of the vector it is given -- so both passes
        // see one sequence, and --max-events cuts both at the same prefix.
        std::vector<std::string> ordered_paths;
        ordered_paths.reserve(chain.size());
        for (const ChainInput& ci : chain) ordered_paths.push_back(ci.path);

        // ---- the seed ------------------------------------------------------
        const std::string chain_sha = pi0::stageB::chain_digest(chain);
        const std::uint64_t seed = pi0::stageB::seed_from_digest(chain_sha);

        std::cout << "stageB_bin\n"
                  << "  inputs     : " << chain.size() << " (in canonical order, by content hash)\n";
        for (const ChainInput& ci : chain) {
            std::cout << "               " << ci.sha256.substr(0, 16) << "...  " << ci.path << '\n';
        }
        std::cout << "  chain sha  : " << chain_sha.substr(0, 16)
                  << (chain.size() == 1 ? "...  (one input: the chain's digest IS its digest)\n"
                                        : "...  (sha256 of the per-input digests, in canonical order)\n")
                  << "  output     : " << args.output << '\n'
                  << "  config     : " << args.config << "  (sha256 " << config_sha.substr(0, 16) << "...)\n"
                  << "  grid A     : " << args.grid_a << '\n'
                  << "  grid B     : " << args.grid_b << '\n'
                  << "  grid hash  : " << binning.provenance_hash() << '\n'
                  << "  bins       : " << n4d << " 4D, " << n3d << " 3D, " << n_a << " Grid A cells, " << n_mgg
                  << " m_gg, " << n_phi << " phi_h\n"
                  << "  m_gg axis  : [" << cuts.mgg_histogram.min_gev << ", " << cuts.mgg_histogram.max_gev << "] GeV\n"
                  << "  mixing     : donors_per_bin = " << cuts.mixing.donors_per_bin
                  << ", seed_mode = " << cuts.mixing.seed_mode << '\n'
                  << "  seed       : 0x" << std::hex << std::setw(16) << std::setfill('0') << seed << std::dec
                  << std::setfill(' ') << "  (first 8 bytes of the CHAIN's sha256 -- its content, not its paths)\n"
                  << "  max events : "
                  << (args.max_events.has_value()
                          ? std::to_string(*args.max_events) +
                                "   *** TRUNCATED RUN (--max-events) *** -- TOTAL across all inputs"
                          : std::string("all"))
                  << "\n\n";

        // ---- Stage A's provenance ------------------------------------------
        // Already read above, before the chain was validated, so a slim file
        // without a provenance block is reported now rather than after two full
        // passes over the whole chain.
        //
        // The WARNINGS below are asserted of the chain as a whole and read from
        // chain.front(). That is not a shortcut and not a last-wins: validate_chain
        // has already REFUSED any chain whose inputs disagree about config.sha256
        // or gbt.fallback_used, so the front's value is every input's value. The
        // things that legitimately differ across a chain -- run, above all -- are
        // carried forward as LISTS below rather than being read off one file.
        const Provenance& stagea_prov = chain.front().prov;
        if (stagea_prov.entries.empty()) {
            std::cout << "WARNING: the input has no /provenance directory. It was not written by this\n"
                      << "         project's stageA_skim, or it was written before provenance existed.\n"
                      << "         Nothing about its cuts, its GBT model or its target can be propagated,\n"
                      << "         and this output will therefore not be traceable past this stage.\n\n";
        } else {
            // Not decoration. A slim file skimmed against a different cuts.json
            // was selected by different cuts from the ones being applied now, and
            // nothing downstream would ever notice. It warns rather than refuses
            // because re-binning an older skim with a newer config is a real and
            // legitimate thing to do -- deliberately, and having been told. (That
            // the CHAIN agrees with ITSELF on this hash is a different question,
            // and that one is a refusal.)
            const std::string a_sha = stagea_prov.get("config.sha256");
            if (!a_sha.empty() && a_sha != config_sha) {
                std::cout << "WARNING: the inputs were skimmed against a DIFFERENT cuts.json.\n"
                          << "           Stage A config.sha256 : " << a_sha << '\n'
                          << "           this run's            : " << config_sha << '\n'
                          << "         The photon and electron selections in those files are not the ones in\n"
                          << "         " << args.config << ". Both hashes are recorded in the output.\n\n";
            }
            if (stagea_prov.get("gbt.fallback_used").rfind("TRUE", 0) == 0) {
                std::cout << "WARNING: the inputs' photons were scored by an RG-A fallback GBT model trained\n"
                          << "         on OTHER data (Stage A provenance: gbt.fallback_used). Any plot made\n"
                          << "         from this output must say so. Propagated to the output.\n\n";
            }
        }
        if (args.allow_truncated_inputs) {
            std::cout << "WARNING: --allow-truncated-inputs was given. At least one input may be a PREFIX of\n"
                      << "         its HIPO, so N_DIS -- the normalisation denominator -- can be short by an\n"
                      << "         unknown amount while the spectra look complete. NO YIELD, RATIO OR\n"
                      << "         NORMALISATION FROM THIS OUTPUT MAY BE QUOTED. Stamped into the output.\n\n";
        }

        // ---- PASS 0: the frozen donor pool ---------------------------------
        // SEQUENTIAL, via TTreeReader rather than RDataFrame. DonorPool's build
        // phase is not thread-safe and not merely as an omission: reservoir
        // sampling is order-dependent by nature, so a parallel pre-pass would make
        // the pool depend on scheduling -- the exact defect it exists to remove.
        //
        // ONE pool over the WHOLE chain: it is constructed here, before the loop,
        // and frozen after it. That is already the natural shape -- the pool never
        // needed to know how many files fed it -- and it is the entire point of
        // --input repeating, since a bin no input reached has no background at all.
        std::cout << "--- pass 0: building the donor pool (sequential, over the whole chain) ---\n";
        pi0::DonorPool pool(cuts.mixing.pool_grid, seed, cuts.mixing.donors_per_bin);

        PoolPassStats p0;
        // The allowance is CHAIN-WIDE and is threaded through every input by
        // reference. Negative = unlimited, matching make_grid's read_slim.
        long long budget = args.max_events.has_value() ? static_cast<long long>(*args.max_events) : -1;
        for (ChainInput& ci : chain) {
            read_pool_pass(ci.path, cuts, pool, p0, budget, ci.n_events_used);
            std::cout << "  read " << std::setw(9) << ci.n_events_used << " events from " << ci.path << '\n';
            if (budget == 0) {
                std::cout << "  (--max-events reached; the remaining inputs are not read. They are stamped\n"
                          << "   into the provenance with n_events_used = 0 -- they are in the chain, hence\n"
                          << "   in the seed, and contributed nothing.)\n";
                break;
            }
        }
        pool.freeze();
        const long long p0_events = p0.events, p0_offered_events = p0.offered_events;
        std::cout << "  events read            : " << p0_events << "  (total, over " << chain.size() << " input(s))\n"
                  << "  events offered         : " << p0_offered_events << '\n'
                  << "  events in no pool bin  : " << p0.no_class
                  << "  (off the pool grid, or no e-gamma-passing photon)\n"
                  << "  pool bins filled       : " << pool.n_filled() << " of " << pool.n_bins() << '\n';
        pool.report_underfilled(std::clog);
        std::cout << '\n';

        // ---- accumulators ---------------------------------------------------
        std::vector<SpectrumCell> spectra(static_cast<std::size_t>(n4d) * static_cast<std::size_t>(n_mgg));
        std::vector<Ptb3dCell> ptb(static_cast<std::size_t>(n3d) * static_cast<std::size_t>(n_mgg));
        std::vector<std::int64_t> n_dis(static_cast<std::size_t>(n_a), 0);

        // Dense in memory, sparse on disk. 1400 * 200 * 12 * 2 = 6.72e6 counters
        // = 54 MB, which is nothing to hold and far cheaper to fill than a hash
        // map on a per-pair path. Only the non-zero rows are written out; see the
        // schema block below.
        //
        //   bsa_index = ((bin4d * n_mgg + imgg) * n_phi + iphi) * 2 + ihel
        //   ihel = 0 for helicity > 0, 1 for helicity < 0
        //
        // That formula is written down HERE, and again in the output's provenance,
        // because the superseded analysis's leaf index formula was written down
        // NOWHERE and had to be reverse-engineered out of its own output files.
        std::vector<std::int64_t> bsa(static_cast<std::size_t>(n4d) * static_cast<std::size_t>(n_mgg) *
                                          static_cast<std::size_t>(n_phi) * 2,
                                      0);

        // ---- PASS 1: RDataFrame --------------------------------------------
        std::cout << "--- pass 1: binning (single-threaded; see the header on why) ---\n";

        long long n_events = 0, n_off_grid_a = 0;
        long long n_gamma_total = 0, n_gamma_pass = 0;
        long long n_same_pairs = 0, n_same_z = 0, n_same_binned = 0;
        long long n_mixed_admissible = 0, n_mixed_z = 0, n_mixed_binned = 0;
        long long n_bsa_filled = 0, n_hel_undef = 0;
        long long n_no_pool_bin = 0;

        // The CANONICAL order, not the command line's: RDataFrame reads a chain in
        // the order of the vector it is given, and pass 0 offered the pool in that
        // same order. Two passes over two orders would truncate at two different
        // prefixes under --max-events and the p0_events check below would catch it
        // -- but only when --max-events is given, which is not when it matters.
        ROOT::RDataFrame df("events", ordered_paths);
        ROOT::RDF::RNode node = df;
        // Range over a multi-file RDataFrame is already N TOTAL across the chain,
        // not N per file -- the same allowance pass 0 spends through budget_take.
        // Range is single-thread only. This program never enables implicit MT, so
        // that costs nothing -- and see the header for why MT is off.
        if (args.max_events.has_value()) node = node.Range(0u, *args.max_events);

        node.Foreach(
            [&](double q2, double xb, double nu, double w, double y, double ex, double ey, double ez, double ee,
                int helicity, const RVecD& gpx, const RVecD& gpy, const RVecD& gpz, const RVecD& gang) {
                ++n_events;

                // ---- N_DIS: ONCE PER EVENT, NOT PER pi0 --------------------
                // *** THE INCLUSIVE DIS DENOMINATOR. *** It is filled here, at the
                // top of the event, before any photon is looked at, and it is the
                // only thing in this program filled per event. The old analysis
                // had a Q2-vs-xB histogram filled PER pi0 CANDIDATE that was
                // structurally indistinguishable from this one and must NOT be
                // used as a normalisation. This program writes no per-candidate
                // (Q^2, x_B) histogram at all -- that is the strongest available
                // way to keep them from being confused. If anyone ever needs the
                // per-candidate count, it is sum(spectra.n_same) over the A cell,
                // which is self-evidently a pi0 count and not a DIS count.
                const int cell_a = binning.A.find(q2, xb);
                if (cell_a >= 0) {
                    ++n_dis[static_cast<std::size_t>(cell_a)];
                } else {
                    ++n_off_grid_a;
                }

                const pi0::DisKin dis{q2, nu, xb, w, y};

                n_gamma_total += static_cast<long long>(gpx.size());
                const std::vector<pi0::Photon> passing = build_photons(gpx, gpy, gpz, gang, cuts);
                n_gamma_pass += static_cast<long long>(passing.size());

                // ---- same-event pairs ---------------------------------------
                const std::vector<pi0::GGPair> pairs = find_gg_pairs(passing, cuts.pairing);
                n_same_pairs += static_cast<long long>(pairs.size());

                for (const pi0::GGPair& p : pairs) {
                    const pi0::SidisKin s =
                        pi0::kin::compute_sidis(p.px, p.py, p.pz, p.e, dis, ex, ey, ez, ee, cuts.beam.energy_gev);

                    // The z window, applied EXPLICITLY rather than left to Grid B's
                    // range to imply. make_grid applies the same cut before it
                    // fits the edges, so the events binned here are the population
                    // the edges were fitted to. Leaving it implicit would make a
                    // cut in cuts.json into a side effect of a machine-generated
                    // edge array -- and the day someone widened the grid, the
                    // selection would change with it.
                    if (!(s.z > cuts.sidis.z_min && s.z < cuts.sidis.z_max)) continue;
                    ++n_same_z;

                    const int imgg = mgg_axis.find(p.mgg);
                    if (imgg < 0) continue;
                    const int bin4d = binning.find_4d(q2, xb, s.z, s.pt2);
                    if (bin4d < 0) continue;
                    ++n_same_binned;

                    // *** THE ABSCISSA FIX. Per (4D bin, m_gg bin). ***
                    SpectrumCell& c =
                        spectra[static_cast<std::size_t>(bin4d) * static_cast<std::size_t>(n_mgg) +
                                static_cast<std::size_t>(imgg)];
                    ++c.n_same;
                    c.sum_q2 += q2;
                    c.sum_xb += xb;
                    c.sum_z += s.z;
                    c.sum_pt2 += s.pt2;

                    // pT broadening, per (3D bin, m_gg bin). find_3d rather than
                    // bin4d / n_pt2: they agree whenever both are in range, but
                    // deriving one from the other is a coincidence to depend on
                    // (core/Binning.hpp says so explicitly).
                    const int bin3d = binning.find_3d(q2, xb, s.z);
                    if (bin3d >= 0) {
                        Ptb3dCell& t = ptb[static_cast<std::size_t>(bin3d) * static_cast<std::size_t>(n_mgg) +
                                           static_cast<std::size_t>(imgg)];
                        ++t.counts;
                        t.sum_pt2 += s.pt2;
                        t.sum_pt4 += s.pt2 * s.pt2;
                    }

                    // BSA, per (4D bin, m_gg bin, phi_h bin, helicity).
                    if (helicity == 0) {
                        // UNDEFINED, not a third state. Dropped from THIS histogram
                        // and from nothing else: this event has already been
                        // counted in N_DIS, in the spectrum and in the sums above.
                        ++n_hel_undef;
                    } else {
                        const int iphi = phi_axis.find(s.phi_h_deg);
                        if (iphi >= 0) {
                            const std::size_t ihel = (helicity > 0) ? 0u : 1u;
                            const std::size_t idx =
                                ((static_cast<std::size_t>(bin4d) * static_cast<std::size_t>(n_mgg) +
                                  static_cast<std::size_t>(imgg)) *
                                     static_cast<std::size_t>(n_phi) +
                                 static_cast<std::size_t>(iphi)) *
                                    2u +
                                ihel;
                            ++bsa[idx];
                            ++n_bsa_filled;
                        }
                    }
                }

                // ---- mixed-event pairs ---------------------------------------
                // EXHAUSTIVE and NON-EXCLUSIVE: every current-event photon against
                // every donor in this event's pool bin. That is NOT what
                // find_gg_pairs does -- it pairs one list against itself, greedily,
                // consuming each photon at most once -- so find_gg_pairs cannot
                // serve this loop and is not called here. What IS shared is the
                // thing that must be: pi0::admissible_pair, the exact predicate
                // find_gg_pairs applies, so the same-event and mixed spectra are
                // cut by one copy of one rule.
                //
                // Exhaustive is also what the old analysis did (cuts.json's
                // /mixing/_comment), so this is not a change of method: the mixed
                // statistics greatly exceed the same-event statistics and the
                // background shape is not statistics-limited.
                //
                // The pool bin is keyed on the count of e-gamma-PASSING photons --
                // the same definition pass 0 used and the same one that governs how
                // many pairs this event makes in the same-event spectrum. One
                // definition everywhere, or the multiplicity matching is decorative.
                const int pb = pool.pool_bin(q2, xb, passing.size());
                if (pb < 0) {
                    ++n_no_pool_bin;
                    return;
                }
                const std::vector<pi0::DonorPhoton>& donors = pool.donors(pb);

                for (const pi0::Photon& g : passing) {
                    for (const pi0::DonorPhoton& d : donors) {
                        // The donor is NOT re-tested against this event's electron.
                        // Correct by construction -- it has no meaningful angle to
                        // it, and it was already tested against its own event's.
                        // e_gamma_deg is set to the donor's own value being
                        // unavailable here; nothing below reads it, and
                        // admissible_pair does not look at it.
                        const pi0::Photon dp{static_cast<double>(d.px), static_cast<double>(d.py),
                                             static_cast<double>(d.pz), 0.0};

                        const std::optional<pi0::GGPair> mp = pi0::admissible_pair(g, dp, cuts.pairing);
                        if (!mp.has_value()) continue;
                        ++n_mixed_admissible;

                        // Same kinematics, computed against THIS event's virtual
                        // photon -- the mixed pair is a fake pi0 in this event.
                        const pi0::SidisKin s = pi0::kin::compute_sidis(mp->px, mp->py, mp->pz, mp->e, dis, ex, ey, ez,
                                                                        ee, cuts.beam.energy_gev);
                        if (!(s.z > cuts.sidis.z_min && s.z < cuts.sidis.z_max)) continue;
                        ++n_mixed_z;

                        const int imgg = mgg_axis.find(mp->mgg);
                        if (imgg < 0) continue;
                        const int bin4d = binning.find_4d(q2, xb, s.z, s.pt2);
                        if (bin4d < 0) continue;
                        ++n_mixed_binned;

                        // n_mixed ONLY. No sum_* is touched here, and that is the
                        // point: a mixed pair's kinematics describe no pi0 and must
                        // never reach the abscissa.
                        ++spectra[static_cast<std::size_t>(bin4d) * static_cast<std::size_t>(n_mgg) +
                                  static_cast<std::size_t>(imgg)]
                              .n_mixed;
                    }
                }
            },
            {"q2", "xb", "nu", "w", "y", "ex", "ey", "ez", "ee", "helicity", "gpx", "gpy", "gpz", "g_e_gamma_deg"});

        // ---- consistency checks ---------------------------------------------
        // Cheap, so checked rather than assumed. A wrong number here is worse than
        // no number, because it is the one everyone quotes.
        {
            std::int64_t ndis_total = 0;
            for (const std::int64_t v : n_dis) ndis_total += v;
            if (ndis_total + n_off_grid_a != n_events) {
                throw std::runtime_error("N_DIS bookkeeping is inconsistent: " + std::to_string(ndis_total) +
                                         " binned + " + std::to_string(n_off_grid_a) + " off Grid A != " +
                                         std::to_string(n_events) + " events read.");
            }
            // THE CHEAPEST TEST OF THE WHOLE MULTI-INPUT PATH, and it must not be
            // weakened. The two passes truncate through INDEPENDENT machinery --
            // pass 0 counts by hand across the per-file loop (budget_take), pass 1
            // hands the chain to RDataFrame::Range -- so this equality is what
            // proves they agree about what "N events" means. A budget scoped per
            // file rather than per chain fails here and nowhere else.
            if (p0_events != n_events) {
                throw std::runtime_error(
                    "the two passes read different numbers of events (" + std::to_string(p0_events) + " vs " +
                    std::to_string(n_events) + ") over " + std::to_string(chain.size()) +
                    " input(s). The donor pool would then be drawn from a different sample than the spectra it "
                    "models. Under --max-events this is what a PER-FILE truncation looks like where a CHAIN-WIDE "
                    "one was meant.");
            }
        }

        // ---- output ----------------------------------------------------------
        std::unique_ptr<TFile> out(TFile::Open(args.output.c_str(), "RECREATE"));
        if (out == nullptr || out->IsZombie()) throw std::runtime_error("cannot create output file: " + args.output);

        // =====================================================================
        // THE OUTPUT SCHEMA -- THE CONTRACT WITH THE PYTHON
        // =====================================================================
        // FOUR TTrees OF FLAT SCALAR BRANCHES, not TH2/TH3.
        //
        // WHY TTrees. The requirement is that uproot read this WITHOUT ROOT, and
        // uproot reads a TTree of scalar branches into numpy arrays with no ROOT
        // installed and no dictionary. That is true of TH2 as well, but a histogram
        // is the wrong shape for what this is: the bins are FLAT indices into a
        // factorized grid, not points on two axes, so a TH2's axes would be
        // meaningless integers and its edges would encode nothing. A TTree names
        // its columns; the reader gets `bin4d` and `imgg` rather than "x" and "y",
        // and the decode is the formula in core/Binning.hpp.
        //
        // The BSA settles it on size alone: dense, it is 1400 x 200 x 12 x 2 =
        // 6.72e6 cells, ~99.9% of them zero. As a histogram that is 54 MB of zeros.
        //
        // ALWAYS DECODE VIA THE INDEX COLUMNS. Every tree carries its own bin
        // indices, so the correct read is a scatter:
        //
        //     import uproot, numpy as np
        //     f = uproot.open("binned.root")
        //     t = f["spectra"].arrays(library="np")
        //     n_same = np.zeros((n4d, n_mgg))
        //     n_same[t["bin4d"], t["imgg"]] = t["n_same"]
        //
        // `spectra`, `ptb3d` and `n_dis` are additionally written DENSE and in
        // index order (row = bin4d * n_mgg + imgg, and so on), so a positional
        // reshape happens to work on them too. `bsa` is SPARSE -- rows with a zero
        // count are omitted -- so on that tree a positional read is silently wrong.
        // Use the scatter everywhere and the distinction never bites.
        //
        // COUNTS ARE INTEGERS (/L), not doubles. The subtraction's variance is
        // sigma_S^2 = N_same + alpha^2 * N_mixed, which is a statement about
        // counting statistics; a count that has been through a float is a count
        // that can no longer be trusted to be one. n_mixed can exceed 2^31 on a
        // production file, hence 64-bit.
        // =====================================================================

        // --- spectra: per (4D bin, m_gg bin) ---------------------------------
        {
            int b_bin4d = 0, b_imgg = 0;
            std::int64_t b_n_same = 0, b_n_mixed = 0;
            double b_sum_q2 = 0, b_sum_xb = 0, b_sum_z = 0, b_sum_pt2 = 0;

            auto* t = new TTree("spectra",
                                "same/mixed m_gg spectra and the SAME-EVENT count-weighted kinematic sums, "
                                "per (4D bin, m_gg bin)");
            t->Branch("bin4d", &b_bin4d, "bin4d/I");
            t->Branch("imgg", &b_imgg, "imgg/I");
            t->Branch("n_same", &b_n_same, "n_same/L");
            t->Branch("n_mixed", &b_n_mixed, "n_mixed/L");
            t->Branch("sum_q2", &b_sum_q2, "sum_q2/D");
            t->Branch("sum_xb", &b_sum_xb, "sum_xb/D");
            t->Branch("sum_z", &b_sum_z, "sum_z/D");
            t->Branch("sum_pt2", &b_sum_pt2, "sum_pt2/D");

            for (int b = 0; b < n4d; ++b) {
                for (int m = 0; m < n_mgg; ++m) {
                    const SpectrumCell& c =
                        spectra[static_cast<std::size_t>(b) * static_cast<std::size_t>(n_mgg) +
                                static_cast<std::size_t>(m)];
                    b_bin4d = b;
                    b_imgg = m;
                    b_n_same = c.n_same;
                    b_n_mixed = c.n_mixed;
                    b_sum_q2 = c.sum_q2;
                    b_sum_xb = c.sum_xb;
                    b_sum_z = c.sum_z;
                    b_sum_pt2 = c.sum_pt2;
                    t->Fill();
                }
            }
            t->Write();
        }

        // --- ptb3d: per (3D bin, m_gg bin) -----------------------------------
        {
            int b_bin3d = 0, b_imgg = 0;
            std::int64_t b_counts = 0;
            double b_sum_pt2 = 0, b_sum_pt4 = 0;

            auto* t = new TTree("ptb3d",
                                "pT-broadening sums per (3D bin = Grid A cell x z bin, m_gg bin). SAME-EVENT ONLY; "
                                "sideband-subtract these exactly like the yield");
            t->Branch("bin3d", &b_bin3d, "bin3d/I");
            t->Branch("imgg", &b_imgg, "imgg/I");
            t->Branch("counts", &b_counts, "counts/L");
            t->Branch("sum_pt2", &b_sum_pt2, "sum_pt2/D");
            t->Branch("sum_pt4", &b_sum_pt4, "sum_pt4/D");

            for (int b = 0; b < n3d; ++b) {
                for (int m = 0; m < n_mgg; ++m) {
                    const Ptb3dCell& c =
                        ptb[static_cast<std::size_t>(b) * static_cast<std::size_t>(n_mgg) + static_cast<std::size_t>(m)];
                    b_bin3d = b;
                    b_imgg = m;
                    b_counts = c.counts;
                    b_sum_pt2 = c.sum_pt2;
                    b_sum_pt4 = c.sum_pt4;
                    t->Fill();
                }
            }
            t->Write();
        }

        // --- n_dis: per Grid A cell ------------------------------------------
        {
            int b_cell_a = 0;
            std::int64_t b_n_dis = 0;
            auto* t = new TTree("n_dis",
                                "INCLUSIVE DIS COUNT per Grid A cell, filled ONCE PER EVENT. This is the "
                                "normalisation denominator. It is NOT a pi0 count and there is deliberately no "
                                "per-pi0-candidate (Q2, xB) histogram in this file");
            t->Branch("cell_a", &b_cell_a, "cell_a/I");
            t->Branch("n_dis", &b_n_dis, "n_dis/L");
            for (int a = 0; a < n_a; ++a) {
                b_cell_a = a;
                b_n_dis = n_dis[static_cast<std::size_t>(a)];
                t->Fill();
            }
            t->Write();
        }

        // --- bsa: per (4D bin, m_gg bin, phi_h bin, helicity) ----------------
        long long bsa_rows = 0;
        {
            int b_bin4d = 0, b_imgg = 0, b_iphi = 0, b_helicity = 0;
            std::int64_t b_counts = 0;
            auto* t = new TTree("bsa",
                                "BSA counts per (4D bin, m_gg bin, phi_h bin, helicity). SPARSE: zero-count cells "
                                "are omitted, so decode via the index columns, never positionally. helicity is +1 "
                                "or -1 only; helicity == 0 is UNDEFINED and enters no row here (it does enter every "
                                "other tree)");
            t->Branch("bin4d", &b_bin4d, "bin4d/I");
            t->Branch("imgg", &b_imgg, "imgg/I");
            t->Branch("iphi", &b_iphi, "iphi/I");
            t->Branch("helicity", &b_helicity, "helicity/I");
            t->Branch("counts", &b_counts, "counts/L");

            for (int b = 0; b < n4d; ++b) {
                for (int m = 0; m < n_mgg; ++m) {
                    for (int ph = 0; ph < n_phi; ++ph) {
                        for (int h = 0; h < 2; ++h) {
                            const std::size_t idx =
                                ((static_cast<std::size_t>(b) * static_cast<std::size_t>(n_mgg) +
                                  static_cast<std::size_t>(m)) *
                                     static_cast<std::size_t>(n_phi) +
                                 static_cast<std::size_t>(ph)) *
                                    2u +
                                static_cast<std::size_t>(h);
                            if (bsa[idx] == 0) continue;
                            b_bin4d = b;
                            b_imgg = m;
                            b_iphi = ph;
                            b_helicity = (h == 0) ? +1 : -1;
                            b_counts = bsa[idx];
                            t->Fill();
                            ++bsa_rows;
                        }
                    }
                }
            }
            t->Write();
        }

        // ---- provenance ------------------------------------------------------
        // N + 1 directories, and none of them are merged.
        //
        //   /provenance              -- what THIS stage did and stands behind.
        //   /provenance_stageA_000   -- what INPUT 0 said about itself, VERBATIM.
        //   /provenance_stageA_001   -- input 1, and so on, in CANONICAL order.
        //
        // The stage's own block is not merged with an input's, because both have
        // an `input.path` and a `config.sha256` and they mean different files.
        // Merging would need renaming, and a renamed provenance value is a value
        // somebody has already interpreted for you. Kept apart, an inherited value
        // can never be mistaken for one this stage asserts. A Stage C repeats the
        // pattern.
        //
        // ONE /provenance_stageA_NNN PER INPUT, and they are not merged either.
        // The number is the CANONICAL index -- the order the pool was actually
        // offered in, which is not the command line's. Merging N Stage A blocks
        // into one is impossible rather than merely untidy: every block has its
        // own `input.path`, `run` and `config.sha256`, Provenance::add appends
        // with no dedup (duplicate keys become ROOT cycles ;1 ;2, and get() would
        // then silently answer with whichever ROOT handed back first), and a
        // single mkdir("provenance_stageA") called N times collides. Stage A pins
        // ONE run per file by refusing a multi-run HIPO, because a header naming
        // one run "would be a lie"; flattening N of those headers into one block
        // would tell exactly that lie one stage later.
        Provenance prov;
        prov.add("stageB_bin.code_version", PI0_CODE_VERSION);
        prov.add("stageB_bin.created_utc", utc_now());
        prov.add("inputs.n", std::to_string(chain.size()));
        // THE ORDERED LIST, and the order is load-bearing rather than cosmetic:
        // the pool is a reservoir, reservoir sampling is a function of the
        // SEQUENCE of offers, so this list -- not the command line -- is what
        // produced the mixed background in this file.
        {
            std::ostringstream os;
            for (std::size_t i = 0; i < chain.size(); ++i) {
                os << (i ? "\n" : "") << i << ": " << chain[i].sha256.substr(0, 16) << "...  " << chain[i].path
                   << "  (n_events_used=" << chain[i].n_events_used << ")";
            }
            prov.add("inputs.ordered", os.str());
        }
        for (std::size_t i = 0; i < chain.size(); ++i) {
            const std::string k = "inputs." + std::to_string(i);
            prov.add(k + ".path", chain[i].path);
            prov.add(k + ".sha256", chain[i].sha256);
            // Per input, not just the chain total. Under --max-events an input can
            // be reached with the allowance spent and contribute NOTHING while
            // still folding into the seed -- so "which files actually fed this
            // pool" must be readable off the artefact rather than re-derived from
            // the flag and the file sizes.
            prov.add(k + ".n_events_used", std::to_string(chain[i].n_events_used));
        }
        prov.add("inputs.order",
                 "ASCENDING BY CONTENT SHA-256, not by path and not as typed. The pool is a reservoir sample and "
                 "reservoir sampling is a function of the SEQUENCE of offers, so `--input a --input b` and "
                 "`--input b --input a` would otherwise build different pools from the same data -- and a farm glob "
                 "promises no order. The key is the CONTENT because sorting by path would make the pool depend on "
                 "where somebody put the files: copy a target's slims to another mount and the order, the offer "
                 "sequence and the background would all move. The list above is the order this run actually used.");
        prov.add("inputs.target",
                 pi0::stageB::prov_lookup(stagea_prov, "target").value_or("unknown (the inputs carry no provenance)") +
                     " -- VALIDATED IDENTICAL ACROSS EVERY INPUT. A chain whose inputs disagreed about their target "
                     "is REFUSED, not warned about: one donor pool is built from all of them and nothing downstream "
                     "of the read knows a photon's target (a donor is four floats), so a mixed chain would make one "
                     "target's photons the mixed background subtracted from another's events -- corrupting both the "
                     "numerator and the denominator of R_A while producing ordinary-looking spectra.");
        {
            // A LIST, because this is the field that legitimately differs across a
            // chain and is the reason --input repeats at all. Stage A refuses a
            // multi-run HIPO (exit 5) precisely so that one file means one run;
            // one output over 30,866 files means 30,866 runs, and the honest
            // record of that is all of them, not the first one encountered.
            const std::vector<std::string> runs = pi0::stageB::distinct_values(chain, "run");
            std::ostringstream os;
            for (std::size_t i = 0; i < runs.size(); ++i) os << (i ? "," : "") << runs[i];
            if (runs.empty()) os << "unknown (the inputs carry no provenance)";
            prov.add("inputs.runs", os.str());
            prov.add("inputs.n_runs", std::to_string(runs.size()));
        }
        prov.add("config.path", args.config);
        prov.add("config.sha256", config_sha);
        prov.add("config.sha256_matches_stageA",
                 stagea_prov.get("config.sha256").empty()
                     ? "unknown (the input carries no provenance)"
                     : (stagea_prov.get("config.sha256") == config_sha
                            ? "true"
                            : "FALSE -- the input was skimmed against a DIFFERENT cuts.json; its electron and "
                              "photon selections are not the ones in this run's config. Both hashes are in this "
                              "file (see provenance_stageA/config.sha256)."));

        prov.add("grid_a.path", args.grid_a);
        prov.add("grid_b.path", args.grid_b);
        // The single most important pair of lines for traceability. The superseded
        // production's edges were never archived and its binning is unrecoverable;
        // this is the number that makes a result matchable back to the geometry
        // that made it.
        prov.add("binning.provenance_hash", binning.provenance_hash());
        prov.add("binning.shape", std::to_string(binning.A.x.nbins()) + " x " + std::to_string(binning.A.y.nbins()) +
                                      " (Grid A: Q2 x xB) * " + std::to_string(binning.B.x.nbins()) + " x " +
                                      std::to_string(binning.B.y.nbins()) + " (Grid B: z x pT2) = " +
                                      std::to_string(n4d) + " 4D bins; " + std::to_string(n3d) + " 3D bins");
        prov.add("binning.index_formula_4d",
                 "bin4d = (i_q2 * n_xb + i_xb) * (n_z * n_pt2) + i_z * n_pt2 + i_pt2. Decode: b = bin4d % (n_z*n_pt2); "
                 "a = bin4d / (n_z*n_pt2); i_pt2 = b % n_pt2; i_z = b / n_pt2; i_xb = a % n_xb; i_q2 = a / n_xb. "
                 "See src/core/Binning.hpp -- written down here because the superseded analysis's leaf formula was "
                 "written down NOWHERE and had to be reverse-engineered out of its output files.");
        prov.add("binning.index_formula_3d", "bin3d = (i_q2 * n_xb + i_xb) * n_z + i_z. NOT bin4d / n_pt2.");
        prov.add("binning.index_formula_bsa",
                 "dense index = ((bin4d * n_mgg + imgg) * n_phi + iphi) * 2 + ihel, ihel = 0 for helicity +1 and 1 "
                 "for helicity -1. The `bsa` tree is SPARSE over this index: a cell absent from the tree is zero. "
                 "Decode via the bin4d/imgg/iphi/helicity columns, never positionally.");
        prov.add("binning.grid_a_source",
                 "see provenance/source inside " + args.grid_a +
                     " -- if it says 'placeholder', these are cuts.json's hand-authored product edges and NOT "
                     "equal-statistics edges fitted to data. No result may be quoted from a placeholder grid.");

        prov.add("mgg_axis", std::to_string(cuts.mgg_histogram.bins) + " uniform bins over [" +
                                 fmt(cuts.mgg_histogram.min_gev) + ", " + fmt(cuts.mgg_histogram.max_gev) +
                                 "] GeV, from cuts.json /mgg_histogram");
        prov.add("phi_h_axis", std::to_string(cuts.bsa.n_phi_bins) + " uniform bins over [" +
                                   fmt(cuts.bsa.phi_min_deg) + ", " + fmt(cuts.bsa.phi_max_deg) +
                                   "] DEGREES, Trento convention, from cuts.json /bsa");
        prov.add("axis.semantics",
                 "every axis in this file is a pi0::Grid1D: half-open [lo, hi), a value exactly on the TOP edge "
                 "lands in the LAST bin, out-of-range and NaN give -1, no flow bins.");

        prov.add("mixing.seed_mode", cuts.mixing.seed_mode);
        {
            std::ostringstream os;
            os << "0x" << std::hex << std::setw(16) << std::setfill('0') << seed;
            prov.add("mixing.seed", os.str());
        }
        prov.add("mixing.chain_sha256", chain_sha);
        prov.add("mixing.seed_derivation",
                 std::string("the first 8 bytes, big-endian, of mixing.chain_sha256 above. ") +
                     (chain.size() == 1
                          ? "This run had ONE input, and a one-input chain's digest IS that input's sha256 -- not a "
                            "hash of it. That is deliberate and load-bearing: it makes a single-input run "
                            "bit-for-bit identical to one made before --input could repeat, so re-hashing does not "
                            "silently move the seed, the pool and the mixed background of every existing result."
                          : "For a chain the digest is the SHA-256 of the per-input sha256 digests concatenated in "
                            "CANONICAL order (inputs.order above). Since that order is itself derived from the "
                            "content, the seed does not depend on the order the files were typed: `--input a "
                            "--input b` and `--input b --input a` give the same seed and the same pool. The digests "
                            "are concatenated rather than xor'd or summed -- those commute, so they would collide "
                            "on any permutation of the same files, which is the case this must distinguish.") +
                     " THE CONTENT, NOT THE PATHS: hashing a path would make the mixed background depend on where "
                     "somebody put the files, so the same data would give two answers from two directories. The "
                     "same data, anywhere, under any names, gives the same pool.");
        prov.add("mixing.donors_per_bin", std::to_string(cuts.mixing.donors_per_bin));
        prov.add("mixing.pool_bins", std::to_string(pool.n_bins()) + " (Q2 x xB x n_photons class, from cuts.json "
                                                                     "/mixing/pool_grid -- NOT Grid A)");
        prov.add("mixing.pool_bins_filled", std::to_string(pool.n_filled()) + " of " + std::to_string(pool.n_bins()) +
                                                " (the shortfall is bins this file never populated; an empty bin "
                                                "means NO mixed background there, which is a different problem from "
                                                "a thin one -- see the stderr report from report_underfilled)");
        prov.add("mixing.pool_events_offered", std::to_string(p0_offered_events) + " of " + std::to_string(p0_events) +
                                                   " events read");
        prov.add("mixing.e_gamma_invariant",
                 "ONLY photons passing their own event's e-gamma cut (pairing.e_gamma_min_angle_deg = " +
                     fmt(cuts.pairing.e_gamma_min_angle_deg) +
                     " deg) were ever offered to the pool. The superseded analysis pooled photons that had FAILED "
                     "that cut, so a photon rejected from the same-event spectrum still shaped the mixed one.");
        prov.add("mixing.asymmetry",
                 "for a mixed pair the e-gamma cut is applied to the CURRENT event's photon ONLY. Correct by "
                 "construction: the donor carries no meaningful angle to this event's electron, and was already cut "
                 "against its own event's. Everything else -- admissibility (pi0::admissible_pair, the same function "
                 "find_gg_pairs calls), the z window, the grids, the m_gg axis -- is identical between the two "
                 "spectra.");
        prov.add("mixing.pairing_shape",
                 "same-event: greedy EXCLUSIVE (find_gg_pairs), each photon consumed at most once. Mixed: EXHAUSTIVE "
                 "and non-exclusive, every current photon against every donor -- as the superseded analysis also did "
                 "(cuts.json /mixing/_comment). The two share the admissibility predicate, not the pairing loop.");
        prov.add("mixing.residual_self_mixing",
                 "an event sampled INTO the pool is later mixed against the pool and can pair with its own photons, "
                 "putting a real pi0 correlation into the background for that event. Bounded by O(donors_per_bin / "
                 "N_events_in_bin), stated in DonorPool.hpp, NOT quantified.");

        prov.add("threading", "single-threaded, deliberately. The frozen pool makes pass 1 safely parallel, but the "
                              "kinematic sums are sums of doubles: per-thread partials merged at the end would make "
                              "sum_q2 depend on how the file was divided, i.e. on the thread count. A reproducible "
                              "abscissa is what this stage is for.");

        prov.add("abscissa.contract",
                 "THE POINT OF THIS STAGE. spectra/sum_{q2,xb,z,pt2} are SAME-EVENT ONLY and are accumulated per "
                 "(4D bin, m_gg bin), weighted by spectra/n_same in the same row. Report a bin at "
                 "sum_X / n_same summed over the SIDEBAND-SUBTRACTED m_gg range -- NOT at the box midpoint, and NOT "
                 "at a mean over the full mass window (which the combinatorial background contaminates). The old "
                 "analysis reported box midpoints: 52.6% of bins sat in a top box, 15.5% implied y > 0.85 -- "
                 "violating the DIS cut their own events passed -- and 6.9% implied nu > E_beam, up to 21.7 GeV "
                 "against a 10.53 GeV beam. It made the attenuation look ~3x weaker than it is.");
        prov.add("abscissa.mixed_excluded",
                 "no mixed pair contributes to any sum_* column. A mixed pair is not a pi0 and its kinematics "
                 "describe nothing.");
        prov.add("n_dis.definition",
                 "INCLUSIVE DIS events per Grid A cell, filled ONCE PER EVENT before any photon is looked at. THE "
                 "NORMALISATION DENOMINATOR. This file contains no per-pi0-candidate (Q2, xB) histogram at all -- "
                 "the old analysis had one that looked exactly like this and must never be used as a normalisation. "
                 "The per-candidate count, if ever wanted, is sum(spectra/n_same) over the A cell.");

        prov.add("events.max_events_requested",
                 args.max_events.has_value()
                     ? std::to_string(*args.max_events) +
                           " -- TRUNCATED RUN: both passes read only this PREFIX of the CHAIN -- N events in TOTAL "
                           "across all inputs, in inputs.order, not N per input -- so no yield or normalisation in "
                           "this file is complete. Inputs late in the order may have contributed nothing; see each "
                           "inputs.N.n_events_used"
                     : std::string("all (no --max-events; every event of every input was read)"));
        // A SEPARATE KEY FROM events.max_events_requested, and deliberately not
        // folded into it. That one is about what THIS stage read; this is about
        // what its INPUTS were, which is a defect this stage cannot undo and did
        // not cause. Stamped whether or not the flag was given, so that "no input
        // was truncated" is an assertion in the file rather than an inference from
        // an absent key -- Provenance::get() cannot tell an absent key from an
        // empty value, so a key that is present only sometimes is a key nobody
        // can check.
        prov.add("inputs.allow_truncated",
                 args.allow_truncated_inputs
                     ? std::string(
                           "TRUE -- *** --allow-truncated-inputs WAS GIVEN. AT LEAST ONE INPUT MAY BE A TRUNCATED "
                           "STAGE A RUN: its event count is a PREFIX of its HIPO, so N_DIS -- the normalisation "
                           "denominator -- can be short by an unknown amount while the spectra look complete. NO "
                           "YIELD, RATIO OR NORMALISATION FROM THIS FILE MAY BE QUOTED. *** This flag exists for "
                           "smoke tests; a production run must not carry it. Which inputs were truncated is in each "
                           "provenance_stageA_NNN/events.max_events_requested.")
                     : std::string("false (the default: every input asserted a full Stage A run, or carried no "
                                   "provenance at all to assert with)"));
        prov.add("events.read", std::to_string(n_events) + " (total, over " + std::to_string(chain.size()) +
                                    " input(s); the per-input split is in inputs.N.n_events_used)");
        prov.add("events.off_grid_a",
                 std::to_string(n_off_grid_a) + " of " + std::to_string(n_events) +
                     " events fell outside Grid A's (Q2, xB) range and are in NO N_DIS cell. They are real DIS "
                     "events that passed every cut; the grid simply does not cover them. N_DIS summed over Grid A "
                     "equals events.read MINUS this number.");
        prov.add("pi0.same_event_binned", std::to_string(n_same_binned));
        prov.add("pi0.helicity_undefined",
                 std::to_string(n_hel_undef) +
                     " binned pi0 had helicity == 0 (UNDEFINED). They are absent from the `bsa` tree and present in "
                     "every other tree.");
        prov.add("pairs.mixed_binned", std::to_string(n_mixed_binned));

        prov.add("schema",
                 "4 TTrees of flat scalar branches, readable by uproot without ROOT. "
                 "spectra(bin4d, imgg, n_same, n_mixed, sum_q2, sum_xb, sum_z, sum_pt2); "
                 "ptb3d(bin3d, imgg, counts, sum_pt2, sum_pt4); "
                 "n_dis(cell_a, n_dis); "
                 "bsa(bin4d, imgg, iphi, helicity, counts). "
                 "spectra/ptb3d/n_dis are DENSE and in index order; `bsa` is SPARSE (zero-count cells omitted). "
                 "ALWAYS decode via the index columns -- on `bsa` a positional read is silently wrong.");

        out->cd();
        prov.write(*out, "provenance");
        // ONE DIRECTORY PER INPUT, numbered by CANONICAL index. Zero-padded to
        // three digits so that a 30,866-file chain lists in ROOT's alphabetical
        // key order in the order it ran, and so that `provenance_stageA_000` is
        // greppable as "the first file offered to the pool". The index is the
        // subscript of inputs.N.path in /provenance -- one number joins them.
        for (std::size_t i = 0; i < chain.size(); ++i) {
            std::ostringstream dir;
            dir << "provenance_stageA_" << std::setw(3) << std::setfill('0') << i;
            out->cd();
            chain[i].prov.write(*out, dir.str().c_str());
        }
        out->cd();
        prov.write_text(*out, "provenance_text");
        // The text form of every input's block, concatenated with a header naming
        // the file each came from -- the one place a human reads the whole chain's
        // inheritance at once, and the reason the per-input blocks can stay
        // separate without becoming unreadable.
        {
            std::ostringstream os;
            for (std::size_t i = 0; i < chain.size(); ++i) {
                os << "=== input " << i << " : " << chain[i].path << " (sha256 " << chain[i].sha256 << ")\n"
                   << chain[i].prov.as_text() << (i + 1 < chain.size() ? "\n" : "");
            }
            out->cd();
            TObjString text(os.str().c_str());
            text.Write("provenance_stageA_text");
        }
        out->Close();

        // ---- report -----------------------------------------------------------
        std::cout << "  events read                    : " << n_events << '\n'
                  << "  events in a Grid A cell (N_DIS): " << (n_events - n_off_grid_a) << '\n'
                  << "  events off Grid A              : " << n_off_grid_a << '\n'
                  << "  photons (Stage A selected)     : " << n_gamma_total << '\n'
                  << "  photons passing e-gamma        : " << n_gamma_pass << "  (> "
                  << cuts.pairing.e_gamma_min_angle_deg << " deg)\n"
                  << "\n  same-event pairs found         : " << n_same_pairs << "  (greedy exclusive, on the "
                  << "e-gamma-passing list)\n"
                  << "  ... in the z window            : " << n_same_z << "  (" << cuts.sidis.z_min << " < z < "
                  << cuts.sidis.z_max << ")\n"
                  << "  ... binned (on the 4D grid)    : " << n_same_binned << '\n'
                  << "  ... with helicity == 0         : " << n_hel_undef << "  (absent from `bsa` ONLY)\n"
                  << "  BSA entries filled             : " << n_bsa_filled << '\n'
                  << "\n  mixed pairs admissible         : " << n_mixed_admissible << '\n'
                  << "  ... in the z window            : " << n_mixed_z << '\n'
                  << "  ... binned                     : " << n_mixed_binned << '\n'
                  << "  events with no pool bin        : " << n_no_pool_bin << '\n'
                  << "\n  bsa rows written (sparse)      : " << bsa_rows << " of "
                  << (static_cast<long long>(n4d) * n_mgg * n_phi * 2) << " dense cells\n";

        std::cout << "\nprovenance written to " << args.output << ":/provenance\n" << prov.as_text();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 2;
    }
}
