#pragma once

/// \file Progress.hpp
/// \brief A progress reporter that behaves ONE way on a terminal and a DIFFERENT
///        way in a batch log, because the same output cannot serve both.
///        Header-only.
///
/// WHY THE TWO MODES ARE NOT OPTIONAL. On a terminal you want a live bar that
/// redraws in place with a carriage return. In a SWIF2/Slurm job that same
/// stream is captured to log.txt, where carriage returns and thousands of redraw
/// frames make the log unreadable and un-greppable -- and the RG-D train files
/// are ~200 GB, so a job runs for a long time producing a lot of frames. So:
///   * interactive (stderr is a TTY): a single line, redrawn in place, throttled.
///   * non-interactive: an occasional plain line ("... 40% (2.0M/5.0M) 320s"),
///     newline-terminated, emitted on a coarse step so the log stays a log.
///
/// THREAD-SAFE BY CONSTRUCTION. The count is atomic, and the redraw is guarded by
/// a try_lock: a worker that cannot get the lock just bumps the counter and moves
/// on, so progress reporting never serialises the workers and never blocks. Only
/// one thread ever writes to the stream at a time. This matters because the whole
/// point of adding this alongside multi-threading is that N workers call add()
/// concurrently.
///
/// The clock and the is-a-tty decision are injectable so the throttle and the
/// formatting can be unit-tested without a real terminal or wall-clock. Nothing
/// here touches ROOT.

#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <mutex>
#include <ostream>
#include <string>
#include <unistd.h>

namespace pi0 {

/// Render a count/total/elapsed into the one line a progress reporter prints.
/// Free function so it can be tested directly. `width` is the bar width in
/// characters for the interactive form; `interactive` picks the form.
[[nodiscard]] inline std::string format_progress(const std::string& label, std::int64_t done,
                                                 std::int64_t total, double elapsed_s,
                                                 bool interactive, int width = 28) {
    char buf[512];
    const double frac = (total > 0) ? static_cast<double>(done) / static_cast<double>(total) : 0.0;
    const double clamped = frac < 0.0 ? 0.0 : (frac > 1.0 ? 1.0 : frac);
    const double rate = (elapsed_s > 0.0) ? done / elapsed_s : 0.0;
    const double eta = (rate > 0.0 && total > done) ? (total - done) / rate : 0.0;

    if (!interactive) {
        // A grep-friendly line. No carriage return, no bar.
        std::snprintf(buf, sizeof buf, "%s %3.0f%% (%lld/%lld) %.0fs elapsed, eta %.0fs",
                      label.c_str(), clamped * 100.0, static_cast<long long>(done),
                      static_cast<long long>(total), elapsed_s, eta);
        return std::string(buf);
    }

    // An in-place bar. The leading '\r' returns to the start of the line; the
    // caller is responsible for the final newline (finish()).
    std::string bar(static_cast<std::size_t>(width), ' ');
    const int filled = static_cast<int>(clamped * width + 0.5);
    for (int i = 0; i < filled && i < width; ++i) bar[static_cast<std::size_t>(i)] = '#';
    std::snprintf(buf, sizeof buf, "\r%s [%s] %3.0f%% %lld/%lld  eta %.0fs   ",
                  label.c_str(), bar.c_str(), clamped * 100.0, static_cast<long long>(done),
                  static_cast<long long>(total), eta);
    return std::string(buf);
}

/// A live progress reporter. Construct with a label and a total; call add()/set()
/// as work completes from any thread; call finish() once at the end.
class Progress {
   public:
    using Clock = std::function<double()>;  ///< seconds since some fixed origin

    /// \param label        shown at the head of every line
    /// \param total        the count 100% corresponds to (0 => percentage omitted)
    /// \param os           where to write (default std::cerr; progress is not data)
    /// \param interactive  tri-state: -1 = auto-detect (isatty on the stream's fd
    ///                     when it is stdout/stderr), 0 = force batch, 1 = force bar
    /// \param clock        injected for testing; default is steady_clock seconds
    explicit Progress(std::string label, std::int64_t total, std::ostream& os,
                      int interactive = -1, Clock clock = default_clock())
        : label_(std::move(label)), total_(total), os_(&os), clock_(std::move(clock)) {
        interactive_ = (interactive < 0) ? stderr_is_tty() : (interactive != 0);
        start_ = clock_();
        last_draw_ = start_ - 1.0;  // force the first draw
    }

    Progress(const Progress&) = delete;
    Progress& operator=(const Progress&) = delete;

    void add(std::int64_t n = 1) { set(done_.fetch_add(n, std::memory_order_relaxed) + n); }

    void set(std::int64_t done) {
        done_.store(done, std::memory_order_relaxed);
        maybe_draw(false);
    }

    /// Draw the final state and, interactively, terminate the line.
    void finish() {
        std::lock_guard<std::mutex> lk(draw_mu_);
        draw_locked(true);
        if (interactive_) (*os_) << "\n";
        os_->flush();
    }

    [[nodiscard]] bool interactive() const { return interactive_; }

   private:
    static Clock default_clock() {
        return [] {
            return std::chrono::duration<double>(
                       std::chrono::steady_clock::now().time_since_epoch())
                .count();
        };
    }

    static bool stderr_is_tty() { return ::isatty(fileno(stderr)) != 0; }

    void maybe_draw(bool force) {
        // A try_lock, not a lock: a worker that loses the race just skips this
        // frame. Progress must never be a synchronisation point for the workers.
        std::unique_lock<std::mutex> lk(draw_mu_, std::try_to_lock);
        if (!lk.owns_lock()) return;

        const double now = clock_();
        if (!force) {
            // Interactive: redraw at most ~10/s. Batch: emit on a coarse step so
            // a long job writes a handful of lines, not thousands.
            if (interactive_) {
                if (now - last_draw_ < 0.1) return;
            } else {
                const std::int64_t done = done_.load(std::memory_order_relaxed);
                if (!batch_step_reached(done, now)) return;
            }
        }
        draw_locked(force);
    }

    /// Batch mode emits a line when EITHER the percentage crosses a 5% step OR
    /// 30 s have passed -- whichever first -- so a fast run is quiet and a slow
    /// one still shows life.
    bool batch_step_reached(std::int64_t done, double now) {
        bool step = false;
        if (total_ > 0) {
            const int pct = static_cast<int>(100.0 * static_cast<double>(done) / static_cast<double>(total_));
            const int bucket = pct / 5;
            if (bucket > last_bucket_) {
                last_bucket_ = bucket;
                step = true;
            }
        }
        if (now - last_draw_ >= 30.0) step = true;
        return step;
    }

    void draw_locked(bool /*final*/) {
        const double now = clock_();
        const std::int64_t done = done_.load(std::memory_order_relaxed);
        (*os_) << format_progress(label_, done, total_, now - start_, interactive_);
        if (!interactive_) (*os_) << "\n";
        os_->flush();
        last_draw_ = now;
    }

    std::string label_;
    std::int64_t total_;
    std::ostream* os_;
    Clock clock_;
    bool interactive_ = false;
    double start_ = 0.0;
    double last_draw_ = 0.0;
    int last_bucket_ = -1;
    std::atomic<std::int64_t> done_{0};
    std::mutex draw_mu_;
};

}  // namespace pi0
