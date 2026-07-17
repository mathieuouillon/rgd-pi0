// Helpers for the appendix sampling-fraction tables. Each table has 4 targets
// × 6 sectors × 8 numerical columns; built-in `table.cell(rowspan: 6)` covers
// the LaTeX `\multirow{6}{*}{LD2}` pattern.
//
// Pass each row as a length-9 array: (sector, a_mu, b_mu, c_mu, d_mu,
// a_sigma, b_sigma, c_sigma, d_sigma). Numeric values can be strings to keep
// scientific-notation formatting identical to the LaTeX source.

#let sf-target-block(target, rows) = {
  let n = rows.len()
  (
    table.cell(rowspan: n, align: horizon)[#target],
    ..rows.flatten(),
  )
}

// Builds a complete sampling-fraction parameter table for one polarity.
// `targets` is an array of (target-name, rows-array) tuples.
// 10 columns of scientific-notation values are dense; we drop the body text
// to 6.5pt (roughly LaTeX `\tiny`) and use tight cell padding so the table
// fits within the page text width.
#let sf-param-table(targets, caption, label) = [
  #figure(
    block[
      #set text(size: 6.5pt)
      #set par(justify: false)
      #table(
        columns: 10,
        align: center + horizon,
        inset: (x: 3pt, y: 3pt),
        // booktabs: thick top, thin under header, none in body
        stroke: (x, y) => (
          top: if y == 0 { 0.9pt } else if y == 1 { 0.4pt } else { none },
          left: 0pt, right: 0pt, bottom: 0pt,
        ),
        fill: (x, y) => if y == 0 { rgb("#e6f3f3") } else { none },
        table.header(
          [*Target*], [*Sector*],
          [$a_mu$], [$b_mu$], [$c_mu$], [$d_mu$],
          [$a_sigma$], [$b_sigma$], [$c_sigma$], [$d_sigma$],
        ),
        ..targets.map(((name, rows)) => sf-target-block(name, rows)).flatten(),
        table.hline(stroke: 0.9pt),
      )
    ],
    caption: caption,
    kind: table,
    supplement: [Table],
  ) #label
]
