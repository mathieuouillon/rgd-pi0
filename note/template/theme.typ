// Visual theme for an HEP analysis note. Clean, professional, lightly accented
// with a deep teal. Booktabs-style tables, hierarchical headings, running
// header showing the current section, muted line numbers, and a polished
// title and outline.

// ---------- Color palette --------------------------------------------------
#let primary       = rgb("#0d7e80")   // deep teal — links, refs, rules
#let primary-dark  = rgb("#0a5d5e")   // section labels, level-3 heading text
#let primary-pale  = rgb("#e6f3f3")   // table-header tint, callout body
#let ink           = rgb("#1a2530")   // body text, level-1 / level-2 heading
#let muted         = rgb("#5a6970")   // captions, secondary metadata
#let rule-c        = rgb("#d6dcdf")   // hairlines, dividers
#let chip-bg       = rgb("#f1f3f5")   // inline raw background

#let accent  = primary                // back-compat alias used by lib.typ
#let teal    = primary

// ---------- Fonts ----------------------------------------------------------
// STIX Two Text / Math is the AMS's house font: a modern, highly legible
// serif designed specifically for scientific publishing. Used here for body,
// title, and all heading levels so the document reads with one consistent
// voice. JetBrains Mono is the gold-standard monospace for code-like inline
// tokens.
#let body-font    = ("STIX Two Text", "Libertinus Serif", "New Computer Modern")
#let heading-font = ("STIX Two Text", "Libertinus Serif", "New Computer Modern")
#let mono-font    = ("JetBrains Mono", "Menlo", "DejaVu Sans Mono")
#let math-font    = ("STIX Two Math", "New Computer Modern Math")

// Document-name string used by the running header. Override via show rule.
#let _doc-name = state("doc-name", "RG-D Common Analysis Note")
#let set-doc-name(name) = _doc-name.update(name)

// Whether to show the running header / page number. Off for title and TOC,
// flipped on for the body.
#let _header-active = state("header-active", false)
#let activate-header()   = _header-active.update(true)
#let deactivate-header() = _header-active.update(false)

// ---------- Theme show-rule -----------------------------------------------
#let theme(body) = {
  // Page geometry & running header / footer
  set page(
    paper: "us-letter",
    margin: (x: 1in, top: 1.15in, bottom: 1in),
    numbering: none,
    header-ascent: 25%,
    footer-descent: 30%,
    header: context {
      // Suppress header on title page and TOC.
      if not _header-active.get() { return }
      // Pick the most recent level-1 heading on or before this page.
      // Prefer one on the current page (catches headings that start the page).
      let here-page = here().page()
      let all = query(heading.where(level: 1))
      let on-page    = all.filter(h => h.location().page() == here-page)
      let before-pg  = all.filter(h => h.location().page() <  here-page)
      let current    = if on-page.len() > 0 { on-page.first() }
                       else if before-pg.len() > 0 { before-pg.last() }
                       else { none }
      let title-text = if current != none {
        let num = if current.numbering != none {
          numbering(current.numbering, ..counter(heading).at(current.location()))
        } else { none }
        if num != none [#num — #current.body] else [#current.body]
      } else { "" }
      block(spacing: 0pt)[
        #grid(
          columns: (1fr, auto),
          align: (left, right),
          text(size: 8.5pt, fill: muted, font: heading-font, smallcaps[#title-text]),
          text(size: 8.5pt, fill: muted, font: heading-font)[#_doc-name.get()],
        )
        #v(2pt)
        #line(length: 100%, stroke: 0.4pt + primary)
      ]
    },
    footer: context {
      if not _header-active.get() { return }
      align(center, text(size: 9pt, fill: muted, font: heading-font)[
        #counter(page).display("1")
      ])
    },
  )

  // Body text
  set text(font: body-font, size: 10.5pt, lang: "en", fill: ink)
  set par(justify: true, leading: 0.65em, first-line-indent: 0pt, spacing: 0.95em)

  // ---------- Headings ----------------------------------------------------
  set heading(numbering: "1.1.1")

  show heading.where(level: 1): it => {
    pagebreak(weak: true)
    v(0.4em)
    block(below: 0.4em)[
      #set text(font: heading-font, size: 18pt, weight: "bold", fill: ink)
      #if it.numbering != none {
        text(fill: primary)[#counter(heading).display(it.numbering)]
        h(0.55em)
      }
      #it.body
    ]
    line(length: 100%, stroke: 1.4pt + primary)
    v(0.6em)
  }

  show heading.where(level: 2): it => {
    block(above: 1.8em, below: 0pt, sticky: true)[
      #set text(font: heading-font, size: 13pt, weight: "semibold", fill: ink)
      #if it.numbering != none {
        text(fill: primary)[#counter(heading).display(it.numbering)]
        h(0.45em)
      }
      #it.body
    ]
    v(0.9em)
  }

  show heading.where(level: 3): it => {
    block(above: 1.4em, below: 0pt, sticky: true)[
      #set text(font: heading-font, size: 10.8pt, weight: "semibold", fill: primary-dark)
      #if it.numbering != none {
        counter(heading).display(it.numbering)
        h(0.45em)
      }
      #it.body
    ]
    v(0.7em)
  }

  // ---------- Links & refs ------------------------------------------------
  show link: set text(fill: primary)
  show ref:  set text(fill: primary)

  // ---------- Captions ----------------------------------------------------
  show figure.caption: it => block(width: 94%)[
    #set text(size: 9.5pt, fill: muted)
    #set par(justify: false, leading: 0.55em)
    #strong(text(fill: primary)[
      #it.supplement~#context it.counter.display(it.numbering).
    ]) #h(0.2em) #it.body
  ]

  // ---------- Inline raw / `code` ----------------------------------------
  show raw.where(block: false): it => box(
    fill: chip-bg,
    inset: (x: 3.5pt, y: 1pt),
    outset: (y: 2pt),
    radius: 2.5pt,
    stroke: 0.4pt + rule-c,
    text(font: mono-font, size: 0.88em, fill: ink)[#it],
  )

  // ---------- Tables: booktabs-style ------------------------------------
  set table(
    stroke: (x, y) => (
      top: if y == 0 { 0.9pt + ink } else if y == 1 { 0.4pt + ink } else { none },
      bottom: 0pt,
      left: 0pt, right: 0pt,
    ),
    align: (x, y) => if y == 0 { center + horizon } else { left + horizon },
    inset: (x: 7pt, y: 5.5pt),
    fill: (x, y) => if y == 0 { primary-pale } else { none },
  )
  // Bottom rule on the last row of every table — applied via a show rule.
  show table: set block(spacing: 1.1em)
  show figure.where(kind: table): set figure(supplement: [Table])
  show table.cell.where(y: 0): set text(weight: "bold", fill: ink, size: 10pt)

  // ---------- Equations -------------------------------------------------
  set math.equation(numbering: n => text(fill: muted, [(#n)]), supplement: [Eq.])
  show math.equation: set text(font: math-font)
  show math.equation: set block(spacing: 0.85em)

  // ---------- Outline (TOC) ---------------------------------------------
  show outline.entry.where(level: 1): set text(weight: "semibold", fill: ink, size: 10.5pt)
  show outline.entry.where(level: 1): set block(above: 0.8em)
  show outline.entry: it => link(it.element.location())[#it.indented(it.prefix(), it.inner())]

  body
}
