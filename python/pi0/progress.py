"""A progress reporter that behaves one way on a terminal and another in a log.

The Python twin of ``src/util/Progress.hpp``, and for the same reason: a live
in-place bar is right on a terminal and wrong in a captured log, where carriage
returns and thousands of redraw frames make the file unreadable. So:

* interactive (stderr is a TTY): one line redrawn in place, throttled to ~10/s;
* non-interactive: an occasional newline-terminated line, on a coarse step.

Progress goes to **stderr**, never stdout -- stdout may carry results a caller
pipes or parses. Depends on nothing outside the standard library.

Usage::

    with progress("skim", total=len(files)) as p:
        for f in files:
            work(f)
            p.add()

or, for a known iterable::

    for f in track(files, "skim"):
        work(f)
"""

from __future__ import annotations

import sys
import time
from contextlib import contextmanager
from typing import Iterable, Iterator, TypeVar

__all__ = ["Progress", "progress", "track"]

T = TypeVar("T")


def _fmt(label: str, done: int, total: int, elapsed: float, interactive: bool, width: int = 28) -> str:
    frac = (done / total) if total > 0 else 0.0
    frac = 0.0 if frac < 0 else (1.0 if frac > 1 else frac)
    rate = done / elapsed if elapsed > 0 else 0.0
    eta = (total - done) / rate if rate > 0 and total > done else 0.0
    if not interactive:
        return f"{label} {frac*100:3.0f}% ({done}/{total}) {elapsed:.0f}s elapsed, eta {eta:.0f}s"
    filled = int(frac * width + 0.5)
    bar = "#" * filled + " " * (width - filled)
    return f"\r{label} [{bar}] {frac*100:3.0f}% {done}/{total}  eta {eta:.0f}s   "


class Progress:
    """Live progress. Call ``add()``/``set()`` as work completes; ``finish()`` once.

    Thread-safe: ``add`` is a single ``+=`` guarded by a lock, and the redraw is
    throttled, so many worker threads calling it stay cheap. Interactivity is
    auto-detected from ``stderr.isatty()`` and can be forced for testing.
    """

    def __init__(self, label: str, total: int, *, stream=None, interactive: bool | None = None,
                 clock=time.monotonic) -> None:
        self.label = label
        self.total = total
        self._stream = stream if stream is not None else sys.stderr
        self._clock = clock
        self.interactive = (self._stream.isatty() if interactive is None
                            and hasattr(self._stream, "isatty") else bool(interactive))
        self._done = 0
        self._start = clock()
        self._last_draw = self._start - 1.0
        self._last_bucket = -1
        import threading
        self._lock = threading.Lock()

    def add(self, n: int = 1) -> None:
        with self._lock:
            self._done += n
            self._maybe_draw()

    def set(self, done: int) -> None:
        with self._lock:
            self._done = done
            self._maybe_draw()

    def _maybe_draw(self) -> None:
        now = self._clock()
        if self.interactive:
            if now - self._last_draw < 0.1:
                return
        else:
            # A line per 5% bucket, or every 30 s -- whichever first.
            step = False
            if self.total > 0:
                bucket = int(100 * self._done / self.total) // 5
                if bucket > self._last_bucket:
                    self._last_bucket = bucket
                    step = True
            if now - self._last_draw >= 30.0:
                step = True
            if not step:
                return
        self._draw(now)

    def _draw(self, now: float) -> None:
        line = _fmt(self.label, self._done, self.total, now - self._start, self.interactive)
        self._stream.write(line)
        if not self.interactive:
            self._stream.write("\n")
        self._stream.flush()
        self._last_draw = now

    def finish(self) -> None:
        with self._lock:
            self._draw(self._clock())
            if self.interactive:
                self._stream.write("\n")
            self._stream.flush()


@contextmanager
def progress(label: str, total: int, **kw) -> Iterator[Progress]:
    p = Progress(label, total, **kw)
    try:
        yield p
    finally:
        p.finish()


def track(it: Iterable[T], label: str, total: int | None = None, **kw) -> Iterator[T]:
    """Wrap an iterable, advancing one step per item."""
    if total is None:
        try:
            total = len(it)  # type: ignore[arg-type]
        except TypeError:
            total = 0
    with progress(label, total, **kw) as p:
        for item in it:
            yield item
            p.add()
