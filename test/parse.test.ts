import fs from 'node:fs'

import { describe, expect, it } from 'vitest'

import {
  extractType,
  parseLines,
  parseRecords,
  parseStringSync,
} from '../src/index.ts'
import { unescape } from '../src/util.ts'

describe('GFF3 parser', () => {
  ;(
    [
      'messy_protein_domains.gff3',
      'gff3_with_syncs.gff3',
      'au9_scaffold_subset.gff3',
      'tomato_chr4_head.gff3',
      'directives.gff3',
      'hybrid1.gff3',
      'hybrid2.gff3',
      'knownGene.gff3',
      'knownGene2.gff3',
      'tomato_test.gff3',
      'spec_eden.gff3',
      'spec_match.gff3',
      'quantitative.gff3',
      'refGene_excerpt.gff3',
      'tair10.gff3',
    ] as const
  ).forEach(filename => {
    it(`can cursorily parse ${filename}`, () => {
      const stuff = parseStringSync(
        fs.readFileSync(`test/data/${filename}`, 'utf8'),
      )
      expect(stuff).toMatchSnapshot()
    })
  })

  it('can parse chr1 TAIR10 gff3', () => {
    parseStringSync(fs.readFileSync('test/data/tair10_chr1.gff', 'utf8'))
  })

  it('parses 0-based start, numeric strand, refName, and lowercased attributes', () => {
    const result = parseStringSync(
      'ctg123\ttest\tgene\t1000\t9000\t0.5\t+\t.\tID=gene00001;Name=TestGene',
    )
    expect(result.length).toBe(1)
    const feature = result[0]!
    expect(feature.refName).toBe('ctg123')
    expect(feature.start).toBe(999) // 0-based (1000 - 1)
    expect(feature.end).toBe(9000)
    expect(feature.strand).toBe(1) // numeric
    expect(feature.score).toBe(0.5)
    expect(feature.type).toBe('gene')
    expect(feature.source).toBe('test')
    expect(feature.id).toBe('gene00001') // lowercased, unpacked
    expect(feature.name).toBe('TestGene') // lowercased, unpacked
    expect(feature.subfeatures).toEqual([])
  })

  it('parses negative and unknown strand correctly', () => {
    const result = parseStringSync(
      `chr1\t.\tgene\t100\t200\t.\t-\t.\tID=g1
chr1\t.\tgene\t300\t400\t.\t.\t.\tID=g2`,
    )
    expect(result[0]!.strand).toBe(-1)
    expect(result[1]!.strand).toBe(0)
  })

  it('parses phase as number', () => {
    const result = parseStringSync('chr1\t.\tCDS\t100\t200\t.\t+\t2\tID=cds1')
    expect(result[0]!.phase).toBe(2)
  })

  it('builds subfeatures from parent/child relationships', () => {
    const result = parseStringSync(
      `ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001
ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00001;Parent=gene00001
ctg123\t.\texon\t1050\t1500\t.\t+\t.\tID=exon1;Parent=mRNA00001`,
    )
    expect(result.length).toBe(1)
    const gene = result[0]!
    expect(gene.id).toBe('gene00001')
    expect(gene.subfeatures.length).toBe(1)
    const mrna = gene.subfeatures[0]!
    expect(mrna.id).toBe('mRNA00001')
    expect(mrna.subfeatures.length).toBe(1)
    expect(mrna.subfeatures[0]!.id).toBe('exon1')
  })

  it('attaches every segment of a multi-location child to its parent', () => {
    const result = parseStringSync(
      `ctgA\t.\tgene\t1\t1000\t.\t+\t.\tID=gene1
ctgA\t.\tmRNA\t1\t1000\t.\t+\t.\tID=mRNA1;Parent=gene1
ctgA\t.\tCDS\t1\t100\t.\t+\t0\tID=cds1;Parent=mRNA1
ctgA\t.\tCDS\t200\t300\t.\t+\t0\tID=cds1;Parent=mRNA1
ctgA\t.\tCDS\t400\t500\t.\t+\t0\tID=cds1;Parent=mRNA1`,
    )
    const mrna = result[0]!.subfeatures[0]!
    const cds = mrna.subfeatures.filter(f => f.type === 'CDS')
    expect(cds.length).toBe(3)
    expect(cds.map(f => f.start)).toEqual([0, 199, 399])
  })

  it('keeps every segment of a top-level discontinuous feature', () => {
    const result = parseStringSync(
      `ctgA\t.\tcDNA_match\t1050\t1500\t5.8e-42\t+\t.\tID=match1
ctgA\t.\tcDNA_match\t5000\t5500\t8.1e-43\t+\t.\tID=match1
ctgA\t.\tcDNA_match\t7000\t9000\t1.4e-40\t+\t.\tID=match1`,
    )
    expect(result.length).toBe(3)
    expect(result.map(f => [f.start, f.end])).toEqual([
      [1049, 1500],
      [4999, 5500],
      [6999, 9000],
    ])
  })

  it('keeps multi-value attributes as arrays', () => {
    const result = parseStringSync(
      'chr1\t.\tgene\t100\t200\t.\t+\t.\tID=g1;Dbxref=GO:123,GO:456',
    )
    expect(result[0]!.dbxref).toEqual(['GO:123', 'GO:456'])
  })

  it('adds suffix to attribute names that conflict with default fields', () => {
    const result = parseStringSync(
      'chr1\t.\tgene\t100\t200\t.\t+\t.\tID=g1;Start=custom_start;Type=custom_type',
    )
    expect(result[0]!.start).toBe(99) // actual start field
    expect(result[0]!.start2).toBe('custom_start') // attribute with suffix
    expect(result[0]!.type).toBe('gene') // actual type field
    expect(result[0]!.type2).toBe('custom_type') // attribute with suffix
  })

  it('takes the first value when ID is multi-valued', () => {
    const result = parseStringSync(
      `ctg\t.\tgene\t1\t10\t.\t+\t.\tID=a,b
ctg\t.\tmRNA\t1\t5\t.\t+\t.\tID=m;Parent=a`,
    )
    expect(result.length).toBe(1)
    expect(result[0]!.id).toEqual(['a', 'b'])
    // Parent=a should match the first ID value
    expect(result[0]!.subfeatures.length).toBe(1)
    expect(result[0]!.subfeatures[0]!.id).toBe('m')
  })

  it('handles escaped characters', () => {
    const result = parseStringSync(
      'SL2.40%25ch01\tIT%25AG\tgene\t100\t200\t.\t+\t.\tID=gene%3B1;Name=Test%20Gene',
    )
    expect(result[0]!.refName).toBe('SL2.40%ch01')
    expect(result[0]!.source).toBe('IT%AG')
    expect(result[0]!.id).toBe('gene;1')
    expect(result[0]!.name).toBe('Test Gene')
  })

  it('parseRecords pairs each top-level feature with its originating record', () => {
    const records = [
      {
        line: 'ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00001;Parent=gene00001',
        offset: 456,
      },
      {
        line: 'ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001',
        offset: 123,
      },
    ]
    const result = parseRecords(records)

    // only the parentless gene is top-level; the mRNA (seen first) is an orphan
    // that gets nested once its parent appears
    expect(result.length).toBe(1)
    expect(result[0]!.record.offset).toBe(123)
    expect(result[0]!.feature.type).toBe('gene')
    expect(result[0]!.feature.subfeatures[0]!.type).toBe('mRNA')
  })

  it('keeps features whose Parent never appears in the input', () => {
    const result = parseStringSync(
      `ctgA\t.\tmRNA\t1\t100\t.\t+\t.\tID=m1;Parent=missing_gene
ctgA\t.\texon\t1\t50\t.\t+\t.\tID=e1;Parent=m1`,
    )
    expect(result.length).toBe(1)
    expect(result[0]!.id).toBe('m1')
    expect(result[0]!.subfeatures[0]!.id).toBe('e1')
  })

  it('does not duplicate a feature that has one resolved and one missing parent', () => {
    const result = parseStringSync(
      `ctgA\t.\tgene\t1\t1000\t.\t+\t.\tID=g1
ctgA\t.\texon\t1\t50\t.\t+\t.\tID=e1;Parent=g1,missing_gene`,
    )
    expect(result.length).toBe(1)
    expect(result[0]!.id).toBe('g1')
    expect(result[0]!.subfeatures.length).toBe(1)
  })

  it('parseLines nests children whose parent appears later', () => {
    const result = parseLines([
      'ctgA\t.\texon\t1\t50\t.\t+\t.\tID=e1;Parent=m1',
      'ctgA\t.\tmRNA\t1\t100\t.\t+\t.\tID=m1;Parent=g1',
      'ctgA\t.\tgene\t1\t1000\t.\t+\t.\tID=g1',
    ])
    expect(result.length).toBe(1)
    expect(result[0]!.id).toBe('g1')
    expect(result[0]!.subfeatures[0]!.id).toBe('m1')
    expect(result[0]!.subfeatures[0]!.subfeatures[0]!.id).toBe('e1')
  })

  it('parseLines returns a feature whose parent never appears', () => {
    const result = parseLines([
      'ctgA\t.\tgene\t1\t1000\t.\t+\t.\tID=g1',
      'ctgA\t.\texon\t1\t50\t.\t+\t.\tID=e1;Parent=missing',
    ])
    // dangling-parent features come after the genuine top-level ones
    expect(result.map(f => f.id)).toEqual(['g1', 'e1'])
  })

  it('parseLines agrees with parseStringSync on a real file', () => {
    const str = fs.readFileSync('test/data/spec_eden.gff3', 'utf8')
    const lines = str
      .split('\n')
      .filter(line => line.length !== 0 && !line.startsWith('#'))
    expect(parseLines(lines)).toEqual(parseStringSync(str))
  })

  it('parseRecords keeps orphans paired with their record', () => {
    const records = [
      { line: 'ctgA\t.\texon\t1\t50\t.\t+\t.\tParent=missing', offset: 7 },
    ]
    const result = parseRecords(records)
    expect(result.length).toBe(1)
    expect(result[0]!.record.offset).toBe(7)
    expect(result[0]!.feature.type).toBe('exon')
  })

  it('suffixes an attribute that would overwrite the subfeatures array', () => {
    const result = parseStringSync(
      `chr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1;Subfeatures=weird
chr1\t.\texon\t1\t50\t.\t+\t.\tParent=g1`,
    )
    expect(result[0]!.subfeatures2).toBe('weird')
    expect(result[0]!.subfeatures.length).toBe(1)
  })

  it('parses truncated lines without throwing', () => {
    const result = parseStringSync('chr1\t.\tgene\t100\t200\t.\t+\t.')
    expect(result.length).toBe(1)
    expect(result[0]!.type).toBe('gene')
    expect(result[0]!.start).toBe(99)
    expect(parseStringSync('chr1')[0]!.refName).toBe('chr1')
  })

  it('ignores an attribute whose value list is empty', () => {
    const result = parseStringSync(
      'chr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1;Foo=,,;Bar=',
    )
    expect(result[0]!.foo).toBeUndefined()
    expect(result[0]!.bar).toBeUndefined()
  })
})

describe('extractType', () => {
  it('reads column 3 of a full line', () => {
    expect(extractType('chr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1')).toBe('gene')
  })

  it('handles lines that end at or before the type column', () => {
    expect(extractType('chr1\t.\tgene')).toBe('gene')
    expect(extractType('chr1\t.')).toBe('')
    expect(extractType('chr1')).toBe('')
  })
})

describe('unescape', () => {
  it('decodes valid escapes and leaves invalid ones literal', () => {
    expect(unescape('SL2.40%25ch01')).toBe('SL2.40%ch01')
    expect(unescape('Test%20Gene')).toBe('Test Gene')
    expect(unescape('no escapes here')).toBe('no escapes here')
    expect(unescape('%2')).toBe('%2')
  })

  it('does not let an invalid escape swallow a following valid one', () => {
    expect(unescape('a%b%20c')).toBe('a%b c')
    expect(unescape('a%20%xy%21b')).toBe('a %xy!b')
  })

  it('decodes escaped multi-byte UTF-8 characters', () => {
    expect(unescape('%E3%81%82%E3%81%82')).toBe('ああ')
    expect(unescape('caf%C3%A9')).toBe('café')
    expect(unescape('a%20%E3%81%82%2Cb')).toBe('a あ,b')
    expect(unescape('%F0%9F%A7%AC')).toBe('🧬')
    const sample = 'あ漢字🧬 é ñ ü ; = , %'
    expect(unescape(encodeURIComponent(sample))).toBe(sample)
  })

  it('replaces bytes that are not valid UTF-8', () => {
    expect(unescape('caf%E9')).toBe('caf�')
  })

  it('decodes escapes regardless of hex digit case', () => {
    expect(unescape('a%2Fb')).toBe('a/b')
    expect(unescape('a%2fb')).toBe('a/b')
    expect(unescape('a%eFb')).toBe(unescape('a%EFb'))
    expect(unescape('a%Efb')).toBe(unescape('a%EFb'))
  })
})
