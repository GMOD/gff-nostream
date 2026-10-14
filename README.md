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

In the browser, pass any GFF3 string, such as one from
`fetch(…).then(r => r.text())`. For a remote `.gff3.gz`, query it with
[@gmod/tabix](https://github.com/GMOD/tabix-js) and pass the returned lines to
`parseLines`.

## Object format

- Coordinates are 0-based half-open.
- `strand` is a number (`1`/`-1`/`0`).
- Attributes become lowercase top-level keys, unwrapped from their array when
  single-valued.
- Child features nest under `subfeatures`.

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

The parser appends `2` to an attribute whose lowercased name matches a fixed
field (`refName`, `source`, `type`, `start`, `end`, `score`, `strand`, `phase`,
`subfeatures`, `seq_id`, `refname`), so `Start=` becomes `start2`.

## Parsing behavior

- The parser ignores comments, directives, and `##FASTA` sections.
- The parser trims spaces around attribute tags, so `Name=A; ID=x` has an `ID`.
- A feature whose `Parent` is missing from the input, as when a tabix query cuts
  off the parent line, comes back as a top-level feature.

### Multi-location features

A multi-location feature has the same ID on several lines, such as a CDS
spanning several segments.

- With no parent, the first line becomes one feature spanning every segment,
  with each segment as a subfeature of the same type. An NCBI `cDNA_match` comes
  back this way.
- With a parent, each line stays a separate feature under it, as for a
  multi-exon CDS. The parser still merges an ID that another line names as its
  parent, such as an NCBI frameshift CDS (SARS-CoV-2 ORF1ab) with mature
  peptides under it.
- Lines that differ in type, or whose first line already has children, stay
  separate duplicates.

## Lazy parsing

Every parse function has a `…Lazy` counterpart that leaves column 9 as raw text
on `feature.attributeString` instead of spreading it into keys. It reads only
`ID` and `Parent`, since it cannot build the tree without them.

```js
import { getAttribute, parseLinesLazy } from '@gmod/gff-nostream'

const features = parseLinesLazy(lines)
const names = features.map(f => getAttribute(f, 'name'))
```

Lazy parsing pays off only when most attributes go unread, since each
`getAttribute` rescans the string. See
[docs/lazy-parsing.md](docs/lazy-parsing.md).

## API

### Parsing

#### `parseStringSync(str: string): GffFeature[]`

Parse a GFF3 string.

#### `parseLines(lines: readonly string[]): GffFeature[]`

Parse GFF3 feature lines, such as a tabix query returns. The lines must contain
no blanks, comments, or `##FASTA` section.

#### `parseRecords<R extends LineRecord>(records: readonly R[]): ParsedRecord<R>[]`

Parse an array of records wrapping raw GFF3 lines.

The parser pairs each top-level feature with its source record, extra fields
included, so a caller can attach its own stable id.

```ts
const records = lines.map((line, i) => ({ line, offset: offsets[i] }))
const features = parseRecords(records).map(({ feature, record }) => ({
  ...feature,
  id: record.offset,
}))
```

#### `parseLinesLazy` / `parseRecordsLazy`

Like `parseLines` and `parseRecords`, returning `LazyGffFeature`s. See
[Lazy parsing](#lazy-parsing).

#### `parseFeatureLazy(line: string): LazyGffFeature`

Parse a single line, with no tree building.

### Attributes

#### `getAttribute(feature: LazyGffFeature, key: string): unknown`

One attribute by its parsed key (lowercased, `2`-suffixed if reserved), in the
shape the eager parser produces.

#### `getAttributes(feature: LazyGffFeature): Record<string, unknown>`

Every attribute, as the eager parser would spread them.

#### `getLinkAttributes(feature: LazyGffFeature): { id: unknown; parent: unknown }`

`ID` and `Parent` in one pass, for callers doing their own linking.

### Other

#### `extractType(line: string): string`

The feature type (column 3) of a raw line, or `''` for a line with fewer than
two tabs.

#### `hasIdAttribute(line: string): boolean`

Whether a raw line has an `ID` the parser would link children to, by the
parser's own tag scan.

### Types

```ts
interface LineRecord {
  line: string
}

interface ParsedRecord<R extends LineRecord = LineRecord> {
  feature: GffFeature
  record: R
}
```

`ParsedLazyRecord` is `ParsedRecord` with a `LazyGffFeature`.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development and release steps.
