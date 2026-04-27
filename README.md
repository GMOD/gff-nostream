# gff-nostream

[![NPM version](https://img.shields.io/npm/v/gff-nostream.svg?style=flat-square)](https://npmjs.org/package/gff-nostream)
[![Build Status](https://img.shields.io/github/actions/workflow/status/GMOD/gff-nostream/push.yml?branch=main)](https://github.com/GMOD/gff-nostream/actions?query=branch%3Amain+workflow%3APush+)

Parse GFF3 data. A simplified version of [@gmod/gff](https://github.com/GMOD/gff-js) with no Node.js stream dependency.

## Install

    $ npm install gff-nostream

## Usage

```js
import { parseStringSync } from 'gff-nostream'
import fs from 'fs'

const stringOfGFF3 = fs.readFileSync('my_annotations.gff3', 'utf8')
const features = parseStringSync(stringOfGFF3)
```

## Object format

In GFF3, features can have more than one location. Features are returned as arrays of all lines sharing the same ID. Values that are `.` in GFF3 are `null` in the output.

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

## API

### `parseStringSync(str: string): GFF3Feature[]`

Synchronously parse a GFF3 string and return an array of features.

### `parseStringSyncJBrowse(str: string): JBrowseFeature[]`

Synchronously parse a GFF3 string and return features in JBrowse format (flat objects with `subfeatures` instead of `child_features`).

### `parseRecords(records: LineRecord[]): GFF3Feature[]`

Parse an array of `LineRecord` objects. Useful when managing raw line data directly (e.g. from an indexed file with byte offsets).

### `parseRecordsJBrowse(records: LineRecord[]): JBrowseFeature[]`

Same as `parseRecords` but returns JBrowse-format features.
