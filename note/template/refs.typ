// Multi-reference helper for the cleveref `\cref{a,b,c}` pattern.
// Typst's native `ref` only takes one label, so we render a comma-and-and list.
// Output mirrors the cleveref convention: "Figs. 5, 6 and 7" (Cref capitalises
// the first word; ref is lower-case-friendly via Typst's figure supplement).

#let refs(..items) = {
  let labels = items.pos()
  let n = labels.len()
  if n == 0 { return }
  if n == 1 { return ref(labels.at(0)) }
  for (i, l) in labels.enumerate() {
    ref(l)
    if i == n - 2 [ and ] else if i < n - 2 [, ]
  }
}
