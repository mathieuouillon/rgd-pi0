"""Tests for pi0.progress: the two output forms and the throttle, without a TTY."""

from __future__ import annotations

import io

from pi0.progress import Progress, _fmt, track


def test_batch_form_is_one_grep_friendly_line():
    s = _fmt("skim", 2000, 5000, 320.0, interactive=False)
    assert "\r" not in s
    assert "40%" in s and "2000/5000" in s and "eta" in s


def test_interactive_form_redraws_in_place_with_a_bar():
    s = _fmt("skim", 2500, 5000, 100.0, interactive=True)
    assert s.startswith("\r")
    assert "#" in s and "50%" in s


def test_zero_total_does_not_divide_by_zero():
    assert "0%" in _fmt("x", 10, 0, 1.0, interactive=False)


def test_percentage_is_clamped():
    assert "100%" in _fmt("x", 7000, 5000, 1.0, interactive=False)


def test_batch_mode_emits_on_a_coarse_step_not_every_update():
    buf = io.StringIO()
    t = [1000.0]  # frozen clock: isolate the 5% percentage trigger from the 30s timer
    p = Progress("load", 100, stream=buf, interactive=False, clock=lambda: t[0])
    for i in range(1, 101):
        p.set(i)
    p.finish()
    lines = buf.getvalue().count("\n")
    assert 18 <= lines <= 25, f"expected ~20 buckets, got {lines}"


def test_interactive_throttles_to_about_ten_frames_a_second():
    buf = io.StringIO()
    t = [0.0]
    p = Progress("run", 1000, stream=buf, interactive=True, clock=lambda: t[0])
    for i in range(1, 101):
        t[0] += 0.01
        p.set(i * 10)
    p.finish()
    frames = buf.getvalue().count("\r")
    assert 8 <= frames <= 14, f"expected ~10 frames, got {frames}"


def test_finish_always_draws_the_final_state():
    buf = io.StringIO()
    t = [0.0]
    p = Progress("run", 10, stream=buf, interactive=False, clock=lambda: t[0])
    p.set(10)
    p.finish()
    assert "100%" in buf.getvalue() and "10/10" in buf.getvalue()


def test_track_advances_one_per_item():
    buf = io.StringIO()
    t = [0.0]
    seen = list(track(range(5), "x", stream=buf, interactive=False, clock=lambda: t[0]))
    assert seen == [0, 1, 2, 3, 4]
    assert "5/5" in buf.getvalue()


def test_add_is_thread_safe():
    import threading

    buf = io.StringIO()
    t = [0.0]
    p = Progress("run", 4000, stream=buf, interactive=False, clock=lambda: t[0])

    def worker():
        for _ in range(1000):
            p.add()

    threads = [threading.Thread(target=worker) for _ in range(4)]
    for th in threads:
        th.start()
    for th in threads:
        th.join()
    p.finish()
    assert "4000/4000" in buf.getvalue()
