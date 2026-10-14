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

The second row matches real use. Every `getAttribute` scans the attribute
string, and a render calls it several times per feature, so by the time a render
has read its attributes the lazy side has lost most of its parsing speedup. A
caller that reads every attribute of every feature should use the eager
functions, which parse each attribute string once rather than once per lookup.

For a caller that holds a whole file resident, the lazy side retains less heap,
and that difference lasts longer than the parse speedup. On ~10 attrs/line, the
eager side retains 16.1MB against 7.0MB lazy (2.3x); on ~1.2 attrs/line, 8.0MB
against 7.1MB.

## Gotchas

`attributeString` is a normal enumerable property. A spread or `JSON.stringify`
of a lazy feature contains that string and _not_ the attributes, so call
`getAttributes` first to parse them into an object.

`getAttribute` takes the _parsed_ key, not the raw GFF3 tag. The parsed key is
lowercased, and suffixed with `2` if it collides with a fixed field.
`getAttribute(f, 'name')` finds `Name=`; a reserved key such as
`getAttribute(f, 'start')` always returns `undefined`, since the parser never
stores an attribute under that key.
