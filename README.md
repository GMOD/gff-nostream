# gff-nostream

[![NPM version](https://img.shields.io/npm/v/gff-nostream.svg?style=flat-square)](https://npmjs.org/package/gff-nostream)
[![Build Status](https://img.shields.io/github/actions/workflow/status/GMOD/gff-nostream/publish.yml?branch=main)](https://github.com/GMOD/gff-nostream/actions/workflows/publish.yml)

Parse GFF3 data. A simplified version of
[@gmod/gff](https://github.com/GMOD/gff-js) with no Node.js stream dependency.

## Install

```sh
pnpm add gff-nostream
```

## Usage

```js
import { readFileSync } from 'node:fs'
import { parseStringSync } from 'gff-nostream'

const features = parseStringSync(readFileSync('my_annotations.gff3', 'utf8'))
```

There is no filesystem dependency — in the browser, pass any GFF3 string, such
as one from `fetch(…).then(r => r.text())`.

## Object format

Features are returned as flat objects with coordinates converted to 0-based
half-open, `strand` as a number (`1`/`-1`/`0`), attributes spread as lowercase
top-level keys, single-valued attributes unwrapped from their array, and child
features nested under `subfeatures`.

A gene with an mRNA child:

```json
{
  "refName": "ctg123",
  "source": null,
  "type": "gene",
  "start": 999,
  "end": 9000,
  "strand": 1,
  "subfeatures": [
    {
      "refName": "ctg123",
      "source": null,
      "type": "mRNA",
      "start": 1049,
      "end": 9000,
      "strand": 1,
      "subfeatures": [],
      "id": "mRNA00001",
      "parent": "gene00001",
      "name": "EDEN.1"
    }
  ],
  "id": "gene00001",
  "name": "EDEN"
}
```

The fixed fields are `refName`, `source`, `type`, `start`, `end`, `score`,
`strand`, `phase`, and `subfeatures`. An attribute whose lowercased name would
land on one of them is suffixed with `2` — `Start=` becomes `start2`. `seq_id`
and `refname` are reserved the same way, so a `Seq_id=` attribute cannot sit
beside `refName`.

## Parsing behavior

These apply to every parse function below.

Comments, directives, and `##FASTA` sections are ignored.

Multi-location features — the same ID on several lines, such as a CDS spanning
several segments — are not merged. Each line is its own flat feature, attached
to its parent (or kept as a top-level item) independently.

A feature whose `Parent` is never defined in the input is returned as a
top-level feature, after the features that appeared in the input, rather than
being dropped. This is common when parsing a slice of a file, e.g. a tabix
region query that cuts off the parent line.

## Lazy parsing

Every parse function has a `…Lazy` counterpart that leaves column 9 as raw text
on `feature.attributeString` instead of spreading it into keys. Only `ID` and
`Parent` are read, since the tree cannot be built without them.

```js
import { getAttribute, parseLinesLazy } from 'gff-nostream'

const features = parseLinesLazy(lines)
const names = features.map(f => getAttribute(f, 'name'))
```

Worth it when most attributes are never read, and not when they are — each
`getAttribute` rescans the string. See
[docs/lazy-parsing.md](docs/lazy-parsing.md) for the measured trade-off and the
gotchas.

## API

### Parsing

#### `parseStringSync(str: string): GffFeature[]`

Parse a GFF3 string.

#### `parseLines(lines: readonly string[]): GffFeature[]`

Parse an array of raw GFF3 feature lines, for a caller that has split and
filtered the file itself — a tabix region query, or a whole-file scan grouping
lines by reference sequence. The lines must already be free of blanks, comments,
and any `##FASTA` section.

#### `parseRecords<R extends LineRecord>(records: readonly R[]): ParsedRecord<R>[]`

Parse an array of records wrapping raw GFF3 lines. Each top-level feature is
returned paired with the record it came from, so a caller can attach its own
stable id (a byte offset, a hash, …) without the parser stamping anything onto
the feature. Records may carry extra fields (`R` is inferred), which pass
through untouched on `record`.

```ts
const records = lines.map((line, i) => ({ line, offset: offsets[i] }))
const features = parseRecords(records).map(({ feature, record }) => ({
  ...feature,
  id: record.offset,
}))
```

#### `parseLinesLazy` / `parseRecordsLazy`

`parseLines` and `parseRecords` returning `LazyGffFeature`s. See
[Lazy parsing](#lazy-parsing).

#### `parseFeatureLazy(line: string): LazyGffFeature`

Parse a single line, with no tree building.

### Attributes

#### `getAttribute(feature: LazyGffFeature, key: string): unknown`

The value of one attribute, by its parsed key (lowercased, `2`-suffixed if
reserved). A single-valued attribute comes back as a string, a multi-valued one
as a string array, and an absent one as `undefined` — the same shapes the eager
parser produces.

#### `getAttributes(feature: LazyGffFeature): Record<string, unknown>`

Every attribute, as the eager parser would have spread them onto the feature.

#### `getLinkAttributes(feature: LazyGffFeature): { id: unknown; parent: unknown }`

`ID` and `Parent` in one pass. This is what the lazy parsers use to build the
tree; it is exported for callers doing their own linking.

### Other

#### `extractType(line: string): string`

Extract the feature type (GFF3 column 3) from a raw line without fully splitting
it. Returns `''` for a line with fewer than two tabs.

### Types

```ts
interface LineRecord {
  line: string
}

interface ParsedRecord<R extends LineRecord = LineRecord> {
  feature: GffFeature
  record: R // the input record this top-level feature was parsed from
}
```

`ParsedLazyRecord` is `ParsedRecord` with a `LazyGffFeature`. `GffFeature` and
`LazyGffFeature` are also exported.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development and release steps.
