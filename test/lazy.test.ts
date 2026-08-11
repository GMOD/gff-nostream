import fs from 'node:fs'

import { describe, expect, it } from 'vitest'

import {
  getAttribute,
  getAttributes,
  parseLines,
  parseLinesLazy,
  parseRecords,
  parseRecordsLazy,
} from '../src/index.ts'

import type { GffFeature, LazyGffFeature } from '../src/index.ts'

const FILES = [
  'messy_protein_domains.gff3',
  'gff3_with_syncs.gff3',
  'au9_scaffold_subset.gff3',
  'tomato_chr4_head.gff3',
  'hybrid1.gff3',
  'hybrid2.gff3',
  'knownGene.gff3',
  'knownGene2.gff3',
  'knownGene_out_of_order.gff3',
  'tomato_test.gff3',
  'spec_eden.gff3',
  'spec_match.gff3',
  'quantitative.gff3',
  'refGene_excerpt.gff3',
  'tair10.gff3',
  'mm9_sample_ensembl.gff3',
  'Saccharomyces_cerevisiae_EF3_e64.gff3',
]

function featureLines(path: string) {
  return fs
    .readFileSync(`test/data/${path}`, 'utf8')
    .split(/\r?\n/)
    .filter(line => line.length !== 0 && !line.startsWith('#'))
}

/**
 * Rebuild the eager parser's output shape from a lazy feature, so the two can
 * be compared directly. This is exactly what a consumer that wants everything
 * has to do, so it doubles as a check that `getAttributes` really is the whole
 * of what the lazy path deferred.
 */
function materialize(f: LazyGffFeature): GffFeature {
  const { attributeString, subfeatures, ...columns } = f
  return {
    ...columns,
    ...getAttributes(f),
    subfeatures: subfeatures.map(materialize),
  }
}

describe('parseLinesLazy', () => {
  it.each(FILES)(
    'materializes to exactly what the eager parser produces: %s',
    file => {
      const lines = featureLines(file)
      expect(parseLinesLazy(lines).map(materialize)).toEqual(parseLines(lines))
    },
  )

  // the lazy path reads ID/Parent with a targeted scan rather than off parsed
  // properties, so the tree it builds is the thing most at risk of diverging
  it.each(FILES)('nests subfeatures identically: %s', file => {
    const lines = featureLines(file)
    const shape = (fs: { subfeatures: unknown[] }[]): unknown =>
      fs.map(f => shape(f.subfeatures as { subfeatures: unknown[] }[]))
    expect(shape(parseLinesLazy(lines))).toEqual(shape(parseLines(lines)))
  })
})

describe('parseRecordsLazy', () => {
  it.each(FILES)('pairs the same records as parseRecords: %s', file => {
    const records = featureLines(file).map((line, i) => ({ line, offset: i }))
    const lazy = parseRecordsLazy(records)
    const eager = parseRecords(records)
    expect(lazy.map(r => r.record)).toEqual(eager.map(r => r.record))
    expect(lazy.map(r => materialize(r.feature))).toEqual(
      eager.map(r => r.feature),
    )
  })
})

// Columns are read by scanning for tabs, with a split-based fallback for a
// line that has fewer than the nine columns GFF3 requires. The 17-file
// differential above only ever exercises the scanning path, since real files
// are well formed.
describe('column parsing against the eager parser', () => {
  const agrees = (line: string) => {
    const [lazy] = parseLinesLazy([line])
    const [eager] = parseLines([line])
    const { attributeString, subfeatures, ...columns } = lazy!
    const { subfeatures: _s, ...eagerAll } = eager!
    // every fixed column the eager parser produced must match
    for (const [k, v] of Object.entries(columns)) {
      expect([k, v]).toEqual([k, eagerAll[k]])
    }
    return lazy!
  }

  it.each([
    ['well formed', 'c\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1'],
    ['all columns "."', 'c\t.\t.\t.\t.\t.\t.\t.\t.'],
    ['negative strand and phase', 'c\tsrc\tCDS\t5\t50\t3.5\t-\t2\tID=x'],
    ['empty columns', 'c\t\t\t1\t100\t\t\t\tID=g1'],
    ['escaped columns', 'c%20x\tsrc%20y\tgene\t1\t100\t.\t+\t.\tID=g1'],
    ['unknown strand character', 'c\tsrc\tgene\t1\t100\t.\t?\t.\tID=g1'],
  ])('matches on a %s line', (_name, line) => {
    agrees(line)
  })

  // fewer than eight tabs takes the fallback; the eager parser reaches the
  // same absent-column handling through split, so the two must still agree
  it.each([
    ['one column', 'ctgA'],
    ['three columns', 'ctgA\tsrc\tgene'],
    ['eight columns', 'c\tsrc\tgene\t1\t100\t.\t+\t.'],
  ])('matches on a truncated line: %s', (_name, line) => {
    const lazy = agrees(line)
    expect(lazy.attributeString).toBe('')
  })

  // a tab inside column 9 is invalid, but `split` cuts the attributes at it and
  // the lazy path must cut them in the same place rather than keeping more
  it('bounds the attribute string exactly where split would', () => {
    const line = 'c\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1\tstray'
    expect(parseLinesLazy([line])[0]!.attributeString).toBe(line.split('\t')[8])
    agrees(line)
  })
})

describe('getAttribute', () => {
  const one = (attrs: string) =>
    parseLinesLazy([`c\ts\tt\t1\t2\t.\t+\t.\t${attrs}`])[0]!
  const eagerOne = (attrs: string) =>
    parseLines([`c\ts\tt\t1\t2\t.\t+\t.\t${attrs}`])[0]!

  // every lookup below is also checked against the eager parser's key, so the
  // two normalizations cannot drift
  const agrees = (attrs: string, key: string) => {
    expect(getAttribute(one(attrs), key)).toEqual(eagerOne(attrs)[key])
  }

  it('reads a plain attribute', () => {
    agrees('gbkey=Src;Name=foo', 'gbkey')
    agrees('gbkey=Src;Name=foo', 'name')
  })

  it('lowercases the tag', () => {
    agrees('GBKEY=Src', 'gbkey')
    agrees('Gene_Name=ABC', 'gene_name')
  })

  it('collapses a single value but keeps multiples as an array', () => {
    agrees('Parent=a', 'parent')
    agrees('Parent=a,b,c', 'parent')
  })

  it('unescapes values', () => {
    agrees('Name=a%20b', 'name')
    agrees('Name=%E3%81%82', 'name')
  })

  it('suffixes a tag colliding with a fixed column', () => {
    agrees('Start=5', 'start2')
    agrees('Type=x', 'type2')
    // the unsuffixed name belongs to the column, so no attribute can be there
    expect(getAttribute(one('Start=5'), 'start')).toBeUndefined()
  })

  // a tag literally named `start2` and a tag named `Start` both normalize to
  // `start2`; the eager parser lets the last one win and so must this
  it('resolves a literal collision with a suffixed name', () => {
    agrees('start2=literal', 'start2')
    agrees('Start=5;start2=literal', 'start2')
    agrees('start2=literal;Start=5', 'start2')
  })

  it('takes the last of a repeated tag, as the eager parser does', () => {
    agrees('Name=first;Name=second', 'name')
  })

  it('skips a tag with no value', () => {
    agrees('Name=;gbkey=Src', 'name')
    agrees('Name=,,;gbkey=Src', 'name')
  })

  it('returns undefined for an absent attribute', () => {
    expect(getAttribute(one('Name=foo'), 'nope')).toBeUndefined()
  })

  it('handles an empty and a "." attribute column', () => {
    expect(getAttribute(one(''), 'name')).toBeUndefined()
    expect(getAttribute(one('.'), 'name')).toBeUndefined()
  })

  it('does not match a tag that merely ends with the key', () => {
    expect(getAttribute(one('transcript_ID=x'), 'id')).toBeUndefined()
    expect(getAttribute(one('transcript_ID=x'), 'transcript_id')).toBe('x')
  })

  // '=' inside a value must not be mistaken for the next attribute's separator
  it('does not read an "=" belonging to a later attribute', () => {
    agrees('Note=a=b;Name=x', 'note')
    agrees('Note=a=b;Name=x', 'name')
  })
})
