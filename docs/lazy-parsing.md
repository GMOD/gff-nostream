# Lazy parsing

`parseLinesLazy`, `parseRecordsLazy`, and `parseFeatureLazy` parse the eight
fixed columns exactly as the eager functions do, but leave column 9 as raw text
on `feature.attributeString`. They read only `ID` and `Parent`, because the
parser needs those two to build the parent/child tree.

## Trade-off

On annotation-grade input such as GENCODE, NCBI or Ensembl, with 15-20
attributes a line, turning column 9 into object keys is about two thirds of the
cost of parsing a line. Deferring it helps when most of those attributes are
never read.

We measured both sides on a generated corpus of distinct genes, with one process
per arm:

| workload                            | ~1.2 attrs/line | ~10 attrs/line |
| ----------------------------------- | --------------- | -------------- |
| parse only                          | 1.2-1.5x        | 1.9-2.2x       |
| parse + the reads a render performs | 1.0-1.2x        | 1.0-1.2x       |

The second row matches real use: a render calls `getAttribute` several times per
feature, each call rescans the attribute string, and the lazy side loses most of
its speedup. A caller that reads every attribute should use the eager functions.

For a caller that holds a whole file resident, the lazy side retains less heap,
and that difference lasts longer than the parse speedup. On ~10 attrs/line, the
eager side retains 16.1MB against 7.0MB lazy (2.3x); on ~1.2 attrs/line, 8.0MB
against 7.1MB.

## Gotchas

`attributeString` is a normal enumerable property. A spread or `JSON.stringify`
of a lazy feature contains that string and _not_ the attributes, so call
`getAttributes` first to parse them into an object.

`getAttribute` takes the _parsed_ key, lowercased and `2`-suffixed if it
collides with a fixed field. `getAttribute(f, 'name')` finds `Name=`, and
`getAttribute(f, 'start')` always returns `undefined`.
