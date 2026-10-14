# @gmod/gff-nostream

[![NPM version](https://img.shields.io/npm/v/@gmod/gff-nostream.svg?style=flat-square)](https://npmjs.org/package/@gmod/gff-nostream)
[![Build Status](https://img.shields.io/github/actions/workflow/status/GMOD/gff-nostream/publish.yml?branch=main)](https://github.com/GMOD/gff-nostream/actions/workflows/publish.yml)

Parses GFF3 data, as a simplified version of
[@gmod/gff](https://github.com/GMOD/gff-js) with no Node.js stream dependency.

## Install

```sh
pnpm add @gmod/gff-nostream
```

## Usage

```js
import { readFileSync } from 'node:fs'
import { parseStringSync } from '@gmod/gff-nostream'

const features = parseStringSync(readFileSync('my_annotations.gff3', 'utf8'))
```

The package has no filesystem dependency.

- In the browser, pass any GFF3 string, such as one from
  `fetch(…).then(r => r.text())`.
- For a remote `.gff3.gz`, query it with
  [@gmod/tabix](https://github.com/GMOD/tabix-js) and read through a filehandle
  from
  [`@gmod/range-cache-filehandle`](https://github.com/GMOD/range-cache-filehandle)
  so its byte ranges are cached. Parse the lines tabix returns with the
  functions here.

## Object format

The parser returns flat objects:

- Coordinates are 0-based half-open.
- `strand` is a number (`1`/`-1`/`0`).
- Attributes become lowercase top-level keys.
- Single-valued attributes are unwrapped from their array.
- Child features nest under `subfeatures`.

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
`strand`, `phase`, and `subfeatures`.

- The parser appends `2` to an attribute whose lowercased name matches a fixed
  field, so `Start=` becomes `start2`.
- It reserves `seq_id` and `refname` the same way, so `Seq_id=` becomes
  `seq_id2` and never sits beside `refName`.

## Parsing behavior

The rules in this section apply to every parse function below.

- The parser ignores comments, directives, and `##FASTA` sections.
- The parser trims spaces around attribute tags before matching them, so
  `Name=A; ID=x` has an `ID`, as it does under @gmod/gff.
- The parser keeps a feature whose `Parent` never appears in the input, and
  returns it as a top-level feature after the ones that did. A slice of a file
  routinely contains such features, for example when a tabix region query cuts
  off the parent line.

### Multi-location features

A multi-location feature has the same ID on several lines, such as a CDS
spanning several segments.

- With a parent, the parser attaches each line to the parent as a separate flat
  feature.
- With no parent, the first line becomes the feature, spanning every segment.
  The parser adds each segment, the first included, under it as a subfeature of
  the same type.
- An NCBI `cDNA_match`, written one line per aligned block, therefore comes back
  as one feature that holds its blocks, while a single-line feature gains no
  child.
- The parser treats a repeated ID whose first line already has children, or
  whose lines differ in type, as a duplicate, and returns each line as a
  separate top-level item.

The parser merges the lines of a repeated ID that has a parent in the same way,
once another line names that ID as its parent.

- NCBI writes a ribosomal-frameshift polyprotein (SARS-CoV-2 ORF1ab, HIV-1
  gag-pol) as one CDS line per reading frame, all under one ID, and gives its
  mature peptides that ID as their parent.
- The parser returns one CDS spanning both frames. Its subfeatures are its
  segments in file order, then the peptides, and the gene holds that one CDS.
- A shared-ID CDS that nothing names as a parent keeps one line per segment, as
  does a further line that names different parents from the first. Every
  multi-exon CDS in an NCBI or Ensembl transcript is the first case.

## Lazy parsing

Every parse function has a `…Lazy` counterpart that leaves column 9 as raw text
on `feature.attributeString` instead of spreading it into keys. It reads only
`ID` and `Parent`, since it cannot build the tree without them.

```js
import { getAttribute, parseLinesLazy } from '@gmod/gff-nostream'

const features = parseLinesLazy(lines)
const names = features.map(f => getAttribute(f, 'name'))
```

Lazy parsing is worth it when most attributes are never read, and not when they
are, because each `getAttribute` rescans the string. See
[docs/lazy-parsing.md](docs/lazy-parsing.md) for the measured trade-off and the
gotchas.

## API

### Parsing

#### `parseStringSync(str: string): GffFeature[]`

Parse a GFF3 string.

#### `parseLines(lines: readonly string[]): GffFeature[]`

Parse an array of raw GFF3 feature lines, for a caller that has split and
filtered the file itself, such as a tabix region query or a whole-file scan that
groups lines by reference sequence. The lines must already be free of blanks,
comments, and any `##FASTA` section.

#### `parseRecords<R extends LineRecord>(records: readonly R[]): ParsedRecord<R>[]`

Parse an array of records wrapping raw GFF3 lines.

- The parser returns each top-level feature paired with its source record and
  adds nothing to the feature itself, so a caller can attach its own stable id
  (a byte offset, a hash, …).
- Records may have extra fields (`R` is inferred), and the parser returns them
  unchanged on `record`.

```ts
const records = lines.map((line, i) => ({ line, offset: offsets[i] }))
const features = parseRecords(records).map(({ feature, record }) => ({
  ...feature,
  id: record.offset,
}))
```

#### `parseLinesLazy` / `parseRecordsLazy`

`parseLinesLazy` and `parseRecordsLazy` work like `parseLines` and
`parseRecords`, and return `LazyGffFeature`s. See [Lazy parsing](#lazy-parsing).

#### `parseFeatureLazy(line: string): LazyGffFeature`

Parse a single line, with no tree building.

### Attributes

#### `getAttribute(feature: LazyGffFeature, key: string): unknown`

The value of one attribute, by its parsed key (lowercased, `2`-suffixed if
reserved). `getAttribute` returns a single-valued attribute as a string, a
multi-valued one as a string array, and an absent one as `undefined`, the same
types the eager parser produces.

#### `getAttributes(feature: LazyGffFeature): Record<string, unknown>`

Every attribute, as the eager parser would have spread them onto the feature.

#### `getLinkAttributes(feature: LazyGffFeature): { id: unknown; parent: unknown }`

`ID` and `Parent` in one pass. The lazy parsers build the tree with it, and the
package exports it for callers doing their own linking.

### Other

#### `extractType(line: string): string`

Extract the feature type (GFF3 column 3) from a raw line without fully splitting
it. Returns `''` for a line with fewer than two tabs.

#### `hasIdAttribute(line: string): boolean`

Whether a raw line has an `ID` the parser would link children to. A line with no
`ID` can have no children. The check uses the parser's own tag scan, so it
treats spacing and case the same way as the linker. A region reader that widens
its fetch to complete a subfeature list can use it to bound the widening.

### Types

```ts
interface LineRecord {
  line: string
}

interface ParsedRecord<R extends LineRecord = LineRecord> {
  feature: GffFeature
  record: R // the input record this top-level feature came from
}
```

`ParsedLazyRecord` is `ParsedRecord` with a `LazyGffFeature`. The package also
exports `GffFeature` and `LazyGffFeature`.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development and release steps.
