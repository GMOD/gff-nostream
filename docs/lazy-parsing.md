# Lazy parsing

`parseLinesLazy`, `parseRecordsLazy`, and `parseFeatureLazy` parse the eight
fixed columns exactly as the eager functions do, but leave column 9 as raw text
on `feature.attributeString`. They read only `ID` and `Parent`, because nothing
can build the parent/child tree without those two.

## When it pays

On annotation-grade input (GENCODE, NCBI, Ensembl — 15-20 attributes a line)
turning column 9 into object keys is about two thirds of the cost of parsing a
line. Deferring it pays when most of those attributes are never read.

Measured against a generated corpus of distinct genes, one process per arm:

| workload                            | ~1.2 attrs/line | ~10 attrs/line |
| ----------------------------------- | --------------- | -------------- |
| parse only                          | 1.2-1.5x        | 1.9-2.2x       |
| parse + the reads a render performs | 1.0-1.2x        | 1.0-1.2x       |

The second row is the one to keep in mind. Every `getAttribute` is a scan, and a
render does several per feature, so at the point of use the lazy side gives most
of the parse win back. Deferring work only pays if the consumer does not then
pay for it on the way out. A caller that reads every attribute of every feature
should use the eager functions, which parse each attribute string once rather
than once per lookup.

Retained heap is the more durable effect for a caller that holds a whole file
resident: on ~10 attrs/line, 16.1MB eager against 7.0MB lazy (2.3x); on ~1.2
attrs/line, 8.0MB against 7.1MB.

## Gotchas

`attributeString` is a normal enumerable property. Spreading or
`JSON.stringify`-ing a lazy feature yields that string and _not_ the attributes
— call `getAttributes` to materialize them first.

`getAttribute` takes the _parsed_ key, not the raw GFF3 tag: lowercased, and
suffixed with `2` if it collides with a fixed field. `getAttribute(f, 'name')`
finds `Name=`; a reserved key such as `getAttribute(f, 'start')` always returns
`undefined`, since no attribute ever lands there.

## Benchmarking

`benchmarks/lazy.bench.ts` runs under `pnpm benchonly`, but read its numbers
with care — the file's header comment explains why vitest's harness flatters
whichever arm runs second, and why sizing a GFF3 fixture up by concatenating a
real file produces a tree shape that flatters the lazy arm specifically. The
figures above come from separate processes over a corpus of distinct genes, not
from that file.
