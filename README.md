# gff-nostream

[![NPM version](https://img.shields.io/npm/v/gff-nostream.svg?style=flat-square)](https://npmjs.org/package/gff-nostream)
[![Build Status](https://img.shields.io/github/actions/workflow/status/GMOD/gff-nostream/publish.yml?branch=main)](https://github.com/GMOD/gff-nostream/actions?query=branch%3Amain+workflow%3APush+)

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

For browser or other non-Node environments, pass any GFF3 string directly — for
example from `fetch`:

```js
import { parseStringSyncJBrowse } from 'gff-nostream'

const text = await fetch('my_annotations.gff3').then(r => r.text())
const features = parseStringSyncJBrowse(text)
```

## Object format

### GFF3 format

Features are returned as arrays of all lines sharing the same ID (to represent
multi-location features). Values that are `.` in GFF3 are `null` in the output.

A simple feature located in one place:

```json
[
  {
    "seq_id": "ctg123",
    "source": null,
    "type": "gene",
    "start": 1000,
    "end": 9000,
    "score": null,
    "strand": "+",
    "phase": null,
    "attributes": {
      "ID": ["gene00001"],
      "Name": ["EDEN"]
    },
    "child_features": [],
    "derived_features": []
  }
]
```

A CDS called `cds00001` located in two places:

```json
[
  {
    "seq_id": "ctg123",
    "source": null,
    "type": "CDS",
    "start": 1201,
    "end": 1500,
    "score": null,
    "strand": "+",
    "phase": "0",
    "attributes": {
      "ID": ["cds00001"],
      "Parent": ["mRNA00001"]
    },
    "child_features": [],
    "derived_features": []
  },
  {
    "seq_id": "ctg123",
    "source": null,
    "type": "CDS",
    "start": 3000,
    "end": 3902,
    "score": null,
    "strand": "+",
    "phase": "0",
    "attributes": {
      "ID": ["cds00001"],
      "Parent": ["mRNA00001"]
    },
    "child_features": [],
    "derived_features": []
  }
]
```

### JBrowse format

The `JBrowse` variants return flat objects with coordinates converted to 0-based
half-open, `strand` as a number (`1`/`-1`/`0`), attributes spread as lowercase
top-level keys, and `subfeatures` instead of `child_features`.

The same gene feature in JBrowse format:

```json
{
  "refName": "ctg123",
  "source": null,
  "type": "gene",
  "start": 999,
  "end": 9000,
  "strand": 1,
  "subfeatures": [],
  "id": "gene00001",
  "name": "EDEN"
}
```

Note: multi-location features (same ID on multiple lines) are not merged in
JBrowse format — only the first occurrence is kept.

## API

### `parseStringSync(str: string): GFF3Feature[]`

Synchronously parse a GFF3 string and return an array of features. Comments,
directives, and `##FASTA` sections are ignored.

### `parseStringSyncJBrowse(str: string): JBrowseFeature[]`

Synchronously parse a GFF3 string and return features in JBrowse format.

### `parseRecords(records: LineRecord[]): GFF3Feature[]`

Parse an array of `LineRecord` objects. Useful when managing raw line data
directly (e.g. from a tabix-indexed file with byte offsets).

### `parseRecordsJBrowse(records: LineRecord[]): JBrowseFeature[]`

Same as `parseRecords` but returns JBrowse-format features.

### `LineRecord`

```ts
interface LineRecord {
  line: string
  hasEscapes: boolean // set true when line contains '%' to enable URL-decoding
  lineHash?: string | number // propagated to attributes._lineHash on the parsed feature
  start?: number // byte offset passthrough (not used by the parser)
  end?: number // byte offset passthrough (not used by the parser)
}
```

## Publishing

[Trusted publishing](https://docs.npmjs.com/about-trusted-publishing) via GitHub
Actions.

```bash
pnpm version patch  # or minor/major
```
