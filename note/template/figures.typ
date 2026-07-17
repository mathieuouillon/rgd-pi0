// Figure helpers built on top of @preview/subpar for native sub-labels.
//
// Each `images` entry is a tuple:
//   (path, caption)            -> no sub-label
//   (path, caption, sublabel)  -> sublabel is a Typst label literal, e.g. <fig:foo-a>
//
// Reference the whole figure with `@fig:foo`, and individual panels with
// the explicit sub-label you provided.

#import "@preview/subpar:0.2.2"

// Internal: turn an image entry into the variadic `figure, <label>, ...`
// sequence subpar.grid expects.
#let _entries(images) = {
  let out = ()
  for entry in images {
    let path = entry.at(0)
    let cap  = entry.at(1)
    out.push(figure(image(path, width: 100%), caption: cap))
    if entry.len() > 2 {
      out.push(entry.at(2))
    }
  }
  out
}

// Single full-width image with explicit width fraction.
// `page` selects a page from a multi-page PDF (default 1).
#let wide-figure(path, caption, label, width: 70%, page: 1) = [
  #figure(image(path, width: width, page: page), caption: caption) #label
]

// 2-up subfigures.
#let subfig2(images, caption, label, gutter: 1em) = subpar.grid(
  ..(_entries(images)),
  columns: (1fr, 1fr),
  gutter: gutter,
  caption: caption,
  label: label,
)

// 3-up subfigures.
#let subfig3(images, caption, label, gutter: 0.6em) = subpar.grid(
  ..(_entries(images)),
  columns: (1fr, 1fr, 1fr),
  gutter: gutter,
  caption: caption,
  label: label,
)

// 2x2 subfigures.
#let subfig2x2(images, caption, label, gutter: 0.6em) = subpar.grid(
  ..(_entries(images)),
  columns: (1fr, 1fr),
  gutter: gutter,
  caption: caption,
  label: label,
)

// Stacked subfigures with constrained-width images.
#let stack-figure(images, caption, label, img-width: 200pt, gutter: 0.7em) = {
  let out = ()
  for entry in images {
    let path = entry.at(0)
    let cap  = entry.at(1)
    out.push(figure(image(path, width: img-width), caption: cap))
    if entry.len() > 2 {
      out.push(entry.at(2))
    }
  }
  subpar.grid(
    ..out,
    columns: (1fr,),
    gutter: gutter,
    caption: caption,
    label: label,
  )
}
