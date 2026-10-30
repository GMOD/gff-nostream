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

How much faster the lazy side is, measured with 5.6.4 on node 24.21, one process
per arm, median of 12 runs:

- maize NAM's gene models: 1.1M lines, ~4 attrs/line
- a GENCODE-shaped file: 135k lines, ~10 attrs/line

| workload                            | maize | GENCODE-shaped |
| ----------------------------------- | ----- | -------------- |
| parse only                          | 1.9x  | 3.2x           |
| parse + `name`, `id`, `gbkey` reads | 1.4x  | 1.8x           |

The second row is closer to real use. A render reads a few attributes per
feature, and each `getAttribute` call rescans the attribute string. A caller
that reads every attribute should use the eager functions. The scripts behind
these figures are not in the repo.

For a caller that holds a whole file resident, the lazy side retains less heap,
and that difference lasts longer than the parse speedup. On an earlier generated
corpus with ~10 attrs/line, the eager side retained 16.1MB against 7.0MB lazy
(2.3x); at ~1.2 attrs/line, 8.0MB against 7.1MB.

## Gotchas

`attributeString` is a normal enumerable property. A spread or `JSON.stringify`
of a lazy feature contains that string and _not_ the attributes, so call
`getAttributes` first to parse them into an object.

`getAttribute` takes the _parsed_ key, lowercased and `2`-suffixed if it
collides with a fixed field. `getAttribute(f, 'name')` finds `Name=`, and
`getAttribute(f, 'start')` always returns `undefined`.
