// Sidebar callouts for highlighting findings, warnings, and key results.
// All callouts use a thin coloured left bar and a faint tint behind the body.
// The `note` variant uses the document's primary teal so callouts feel native
// to the rest of the theme.

#import "@preview/showybox:2.0.4": showybox

#let _bar(title, body, accent, bg, glyph) = showybox(
  frame: (
    border-color: accent,
    title-color: bg,
    body-color: bg,
    thickness: (left: 3pt, rest: 0pt),
    radius: 2pt,
    inset: (x: 11pt, y: 8pt),
  ),
  title-style: (
    color: accent.darken(15%),
    weight: "bold",
    boxed-style: none,
  ),
  title: [#text(font: ("Inter", "Helvetica Neue", "Helvetica"))[#glyph] #h(0.45em) #title],
  body,
)

// Teal — same accent as the rest of the theme.
#let note-box(title: "Note", body) = _bar(
  title, body,
  rgb("#0d7e80"),
  rgb("#e6f3f3"),
  "ⓘ",
)

// Amber — for caveats / things-to-revisit.
#let important-box(title: "Important", body) = _bar(
  title, body,
  rgb("#b88600"),
  rgb("#fff7d6"),
  "!",
)

// Deep green — for confirmed results / measurement values.
#let result-box(title: "Result", body) = _bar(
  title, body,
  rgb("#1a7f37"),
  rgb("#dafbe1"),
  "✓",
)

// Crimson — for known problems / blockers.
#let warning-box(title: "Warning", body) = _bar(
  title, body,
  rgb("#b42318"),
  rgb("#fee4e2"),
  "⚠",
)
