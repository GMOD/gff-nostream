# gff-nostream

[![NPM version](https://img.shields.io/npm/v/gff-nostream.svg?style=flat-square)](https://npmjs.org/package/gff-nostream)
![Build Status](https://img.shields.io/github/actions/workflow/status/GMOD/gff-nostream/publish.yml?branch=main)

Parse GFF3 data. A simplified version of
[@gmod/gff](https://github.com/GMOD/gff-js) with no Node.js stream dependency.

## Install

    $ pnpm add gff-nostream

## Usage

```js
import { parseStringSync } from 'gff-nostream'
import fs from 'fs'

const features = parseStringSync(fs.readFileSync('my_annotations.gff3', 'utf8'))
```

In the browser or other non-Node environments, pass any GFF3 string directly —
for example from `fetch`:

```js
import { parseStringSync } from 'gff-nostream'

const text = await fetch('my_annotations.gff3').then(r => r.text())
const features = parseStringSync(text)
```

## Object format

Features are returned as flat objects with coordinates converted to 0-based
half-open, `strand` as a number (`1`/`-1`/`0`), attributes spread as lowercase
top-level keys, single-valued attributes unwrapped from their array, and child
features nested under `subfeatures`. An attribute whose lowercased name collides
with a built-in field (`start`, `end`, `seq_id`, `refname`, `score`, `type`,
`source`, `phase`, `strand`, `subfeatures`) is suffixed with `2` (`start2`,
`type2`).

A gene with an mRNA child:

```json
{
  "refName": "ctg123",
  "source": null,
  "type": "gene",
  "start": 999,
  "end": 9000,
  "strand": 1,
  "id": "gene00001",
  "name": "EDEN",
  "subfeatures": [
    {
      "refName": "ctg123",
      "source": null,
      "type": "mRNA",
      "start": 1049,
      "end": 9000,
      "strand": 1,
      "id": "mRNA00001",
      "parent": "gene00001",
      "subfeatures": []
    }
  ]
}
```

Multi-location features (the same ID on multiple lines, such as a CDS spanning
several segments) are not merged — each line is its own flat feature, attached
to its parent (or kept as a top-level item) independently.

A feature whose `Parent` is never defined in the input — common when parsing a
slice of a file, e.g. a tabix region query that cuts off the parent line — is
returned as a top-level feature, after the features that appeared in the input,
rather than being dropped.

## API

### `parseStringSync(str: string): GffFeature[]`

Synchronously parse a GFF3 string and return an array of features. Comments,
directives, and `##FASTA` sections are ignored.

### `parseRecords<R>(records: readonly R[]): ParsedRecord<R>[]`

Parse an array of records wrapping raw GFF3 lines. Useful when managing raw line
data directly (e.g. from a tabix-indexed file). Each top-level feature is
returned paired with the record it came from, so a caller can attach its own
stable id (a byte offset, a hash, …) without the parser stamping anything onto
the feature. Records may carry extra fields (`R` is inferred), which pass
through untouched on `record`.

```ts
const features = parseRecords(lines.map(line => ({ line, offset }))).map(
  ({ feature, record }) => ({ ...feature, id: record.offset }),
)
```

### `extractType(line: string): string`

Extract the feature type (GFF3 column 3) from a raw line without fully splitting
it.

### `LineRecord` / `ParsedRecord`

```ts
interface LineRecord {
  line: string
}

interface ParsedRecord<R extends LineRecord = LineRecord> {
  feature: GffFeature
  record: R // the input record this top-level feature was parsed from
}
```

## Publishing

[Trusted publishing](https://docs.npmjs.com/about-trusted-publishing) via GitHub
Actions.

```bash
pnpm version patch  # or minor/major
```
