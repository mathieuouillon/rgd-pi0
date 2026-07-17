// Public surface of the analysis-note template.
// Body files import from here: `#import "../template/lib.typ": *`

#import "theme.typ": theme, accent, muted, primary, primary-dark, primary-pale, ink, rule-c, set-doc-name, activate-header, deactivate-header, heading-font
#import "figures.typ": wide-figure, subfig2, subfig3, subfig2x2, stack-figure
#import "callouts.typ": note-box, important-box, result-box, warning-box
#import "glossary.typ": acronyms, glossary-section
#import "refs.typ": refs
#import "tables.typ": sf-param-table
#import "@preview/glossarium:0.5.10": make-glossary, register-glossary, print-glossary, gls, glspl
#import "@preview/physica:0.9.8": *

// Document show-rule. Used in main.typ as
//   `#show: note.with(title: ..., author: ..., email: ...)`
#let note(
  title: "Analysis Note",
  subtitle: none,
  author: "",
  affiliation: "",
  email: "",
  date: datetime.today(),
  version: none,                 // e.g. "v1.0", "Draft"
  doc-name: none,                 // running-header text; defaults to title
  line-numbers: true,
  glossary-at-end: true,
  body,
) = {
  show: theme

  if doc-name != none {
    set-doc-name(doc-name)
  } else {
    set-doc-name(title)
  }

  // Subtle line numbering, suppressed inside figures, headings, and tables.
  if line-numbers {
    set par.line(numbering: n => text(size: 6.5pt, fill: rule-c)[#n])
    show figure: set par.line(numbering: none)
    show heading: set par.line(numbering: none)
    show table:   set par.line(numbering: none)
  }

  // ---------- Title page ------------------------------------------------
  v(2.5em)
  align(center)[
    #line(length: 55%, stroke: 0.8pt + primary)
    #v(0.7em)
    #text(size: 22pt, weight: "bold", fill: ink, font: heading-font)[#title]
    #if subtitle != none {
      v(0.55em)
      text(size: 13pt, style: "italic", fill: muted)[#subtitle]
    }
    #v(0.5em)
    #line(length: 55%, stroke: 0.8pt + primary)
    #v(2.2em)
    #text(size: 11.5pt, weight: "semibold")[#author]
    #if affiliation != "" {
      linebreak()
      text(size: 10pt, fill: muted)[#affiliation]
    }
    #if email != "" {
      linebreak()
      text(size: 10pt)[#link("mailto:" + email)[#email]]
    }
    #v(0.7em)
    #text(size: 9.5pt, fill: muted)[#date.display("[month repr:long] [day], [year]")]
    #if version != none {
      linebreak()
      text(size: 9pt, fill: muted, style: "italic")[#version]
    }
  ]

  pagebreak()

  // ---------- Outline ---------------------------------------------------
  // Render the title manually as a styled block (NOT a heading) so it
  // doesn't get picked up by the running-header heading query.
  block(spacing: 0pt)[
    #set text(font: heading-font, size: 14pt, weight: "bold", fill: ink)
    Contents
    #v(0.3em)
    #line(length: 100%, stroke: 1pt + primary)
  ]
  v(0.6em)
  outline(title: none, indent: auto, depth: 3)

  pagebreak(weak: true)

  // Re-enable running header & page numbers for the body
  activate-header()
  counter(page).update(1)

  body

  if glossary-at-end {
    pagebreak(weak: true)
    set heading(numbering: none)
    heading(level: 1)[Acronyms]
    print-glossary(acronyms, show-all: true, disable-back-references: true)
  }
}
