// Unit tests for the progress reporter's pure logic: the two output forms, and
// the throttle that keeps a batch log from filling with redraw frames.
//
// The terminal behaviour itself (carriage returns, isatty) is exercised by the
// end-to-end runs; here we test what can be tested without a TTY, using an
// injected clock and a forced interactive/batch mode.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <sstream>
#include <string>

#include "util/Progress.hpp"

using Catch::Matchers::ContainsSubstring;
using pi0::format_progress;
using pi0::Progress;

TEST_CASE("batch form is a single grep-friendly line, no carriage return", "[progress]") {
    const std::string s = format_progress("skim", 2000, 5000, 320.0, /*interactive=*/false);
    CHECK(s.find('\r') == std::string::npos);
    CHECK_THAT(s, ContainsSubstring("40%"));
    CHECK_THAT(s, ContainsSubstring("2000/5000"));
    CHECK_THAT(s, ContainsSubstring("eta"));
}

TEST_CASE("interactive form redraws in place with a bar", "[progress]") {
    const std::string s = format_progress("skim", 2500, 5000, 100.0, /*interactive=*/true);
    CHECK(s.front() == '\r');           // returns to line start
    CHECK_THAT(s, ContainsSubstring("#"));   // has a bar
    CHECK_THAT(s, ContainsSubstring("50%"));
}

TEST_CASE("percentage is clamped and a zero total does not divide by zero", "[progress]") {
    CHECK_THAT(format_progress("x", 10, 0, 1.0, false), ContainsSubstring("0%"));   // total 0
    CHECK_THAT(format_progress("x", 7000, 5000, 1.0, false), ContainsSubstring("100%"));  // over
    CHECK_THAT(format_progress("x", -5, 5000, 1.0, true), ContainsSubstring("0%"));  // under
}

TEST_CASE("batch mode emits on a coarse step, not on every update", "[progress]") {
    // A frozen clock, so nothing is emitted on the 30s timer -- only the 5%
    // percentage steps fire. 100 updates over a total of 100 should therefore
    // produce ~20 lines (one per 5% bucket), not 100.
    std::ostringstream os;
    double t = 1000.0;  // never advances: isolate the percentage trigger
    Progress p("load", 100, os, /*interactive=*/0, [&] { return t; });
    for (int i = 1; i <= 100; ++i) p.set(i);
    p.finish();
    std::size_t lines = 0;
    for (char ch : os.str()) if (ch == '\n') ++lines;
    INFO("emitted:\n" << os.str());
    CHECK(lines >= 18);   // ~20 buckets
    CHECK(lines <= 25);   // definitely not ~100
}

TEST_CASE("interactive mode throttles to about ten frames a second", "[progress]") {
    // Advance the clock 10 ms per update over 1 s of updates (100 of them). At a
    // 100 ms throttle that is ~10 redraws, plus the forced first and the finish.
    std::ostringstream os;
    double t = 0.0;
    Progress p("run", 1000, os, /*interactive=*/1, [&] { return t; });
    for (int i = 1; i <= 100; ++i) {
        t += 0.01;
        p.set(i * 10);
    }
    p.finish();
    // Count carriage returns as redraw frames.
    std::size_t frames = 0;
    for (char ch : os.str()) if (ch == '\r') ++frames;
    INFO("frames: " << frames);
    CHECK(frames >= 8);
    CHECK(frames <= 14);
}

TEST_CASE("finish always draws the final state even if throttled", "[progress]") {
    std::ostringstream os;
    double t = 0.0;
    Progress p("run", 10, os, /*interactive=*/0, [&] { return t; });
    p.set(10);        // 100%, but the clock has not moved and no bucket may have flushed
    p.finish();
    CHECK_THAT(os.str(), ContainsSubstring("100%"));
    CHECK_THAT(os.str(), ContainsSubstring("10/10"));
}

TEST_CASE("add() accumulates the same as an equivalent set()", "[progress]") {
    std::ostringstream os;
    double t = 0.0;
    Progress p("run", 100, os, /*interactive=*/0, [&] { return t; });
    for (int i = 0; i < 40; ++i) p.add(1);
    p.finish();
    CHECK_THAT(os.str(), ContainsSubstring("40/100"));
}
