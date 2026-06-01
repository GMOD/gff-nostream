import fs from 'node:fs'

import { describe, expect, it } from 'vitest'

import {
  parseRecords,
  parseRecordsJBrowse,
  parseStringSync,
  parseStringSyncJBrowse,
} from '../src/index.ts'

describe('GFF3 parser', () => {
  it('can parse gff3_with_syncs.gff3', () => {
    const stuff = parseStringSync(
      fs.readFileSync('test/data/gff3_with_syncs.gff3', 'utf8'),
    )
    expect(stuff).toMatchSnapshot()
  })
  ;(
    [
      [1010, 'messy_protein_domains.gff3'],
      [4, 'gff3_with_syncs.gff3'],
      [51, 'au9_scaffold_subset.gff3'],
      [14, 'tomato_chr4_head.gff3'],
      [5, 'directives.gff3'],
      [5, 'hybrid1.gff3'],
      [3, 'hybrid2.gff3'],
      [6, 'knownGene.gff3'],
      [6, 'knownGene2.gff3'],
      [16, 'tomato_test.gff3'],
      [3, 'spec_eden.gff3'],
      [1, 'spec_match.gff3'],
      [8, 'quantitative.gff3'],
    ] as const
  ).forEach(([_count, filename]) => {
    it(`can cursorily parse ${filename}`, () => {
      const stuff = parseStringSync(
        fs.readFileSync(`test/data/${filename}`, 'utf8'),
      )
      expect(stuff).toMatchSnapshot()
    })
  })

  it('can parse the EDEN gene from the gff3 spec', () => {
    const stuff = parseStringSync(
      fs.readFileSync('test/data/spec_eden.gff3', 'utf8'),
    )
    expect(stuff).toMatchSnapshot()
  })

  it('can parse an excerpt of the refGene gff3', () => {
    const stuff = parseStringSync(
      fs.readFileSync('test/data/refGene_excerpt.gff3', 'utf8'),
    )
    expect(stuff).toMatchSnapshot()
  })

  it('can parse an excerpt of the TAIR10 gff3', () => {
    const stuff = parseStringSync(
      fs.readFileSync('test/data/tair10.gff3', 'utf8'),
    )
    expect(stuff).toMatchSnapshot()
  })

  it('can parse chr1 TAIR10 gff3', () => {
    parseStringSync(fs.readFileSync('test/data/tair10_chr1.gff', 'utf8'))
  })

  it('can parse a string synchronously', () => {
    const gff3 = fs.readFileSync('test/data/spec_eden.gff3').toString('utf8')
    const result = parseStringSync(gff3)
    expect(result).toMatchSnapshot()
  })

  it('can parse LineRecord objects with lineHash', () => {
    const records = [
      {
        line: 'ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001',
        lineHash: 'hash123',
        start: 1000,
        end: 9000,
        hasEscapes: false,
      },
      {
        line: 'ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00001;Parent=gene00001',
        lineHash: 456,
        start: 1050,
        end: 9000,
        hasEscapes: false,
      },
    ]
    const result = parseRecords(records)
    expect(result.length).toBe(1)
    expect(result[0]?.[0]?.attributes?._lineHash).toEqual(['hash123'])
    expect(
      result[0]?.[0]?.child_features[0]?.[0]?.attributes?._lineHash,
    ).toEqual(['456'])
  })

  it('can parse LineRecord objects without lineHash', () => {
    const records = [
      {
        line: 'ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001',
        start: 1000,
        end: 9000,
        hasEscapes: false,
      },
    ]
    const result = parseRecords(records)
    expect(result.length).toBe(1)
    expect(result[0]?.[0]?.attributes?._lineHash).toBeUndefined()
    expect(result[0]?.[0]?.attributes?.ID).toEqual(['gene00001'])
  })

  it('can parse LineRecord objects with empty attributes', () => {
    const records = [
      {
        line: 'ctg123\t.\tgene\t1000\t9000\t.\t+\t.\t.',
        lineHash: 'hashOnly',
        start: 1000,
        end: 9000,
        hasEscapes: false,
      },
    ]
    const result = parseRecords(records)
    expect(result.length).toBe(1)
    expect(result[0]?.[0]?.attributes?._lineHash).toEqual(['hashOnly'])
  })

  it('can parse some whitespace and escapes', () => {
    const gff3 = `
SL2.40%25ch01	IT%25AG eugene	g%25e;ne	80999140	81004317	.	+	.	multivalue=val1,val2,val3;testing=blah
`

    const result = parseStringSync(gff3)
    const referenceResult = [
      [
        {
          seq_id: 'SL2.40%ch01',
          source: 'IT%AG eugene',
          type: 'g%e;ne',
          start: 80999140,
          end: 81004317,
          score: null,
          strand: '+',
          phase: null,
          attributes: {
            multivalue: ['val1', 'val2', 'val3'],
            testing: ['blah'],
          },
          child_features: [],
          derived_features: [],
        },
      ],
    ]

    expect(result).toEqual(referenceResult)
  })

  it('shares child_features across all lines of a multi-location feature', () => {
    const records = [
      'ctg\t.\tgene\t1\t10\t.\t+\t.\tID=g',
      'ctg\t.\tmRNA\t1\t5\t.\t+\t.\tID=m1;Parent=g',
      'ctg\t.\tgene\t20\t30\t.\t+\t.\tID=g',
      'ctg\t.\tmRNA\t20\t25\t.\t+\t.\tID=m2;Parent=g',
    ].map(line => ({
      line,
      start: 0,
      end: 0,
      hasEscapes: false,
    }))

    const result = parseRecords(records)
    expect(result.length).toBe(1)
    const gene = result[0]!
    expect(gene.length).toBe(2)
    // both lines of the multi-location gene see all children, regardless of
    // whether the child arrived before or after the duplicate-ID line
    expect(gene[0]!.child_features).toHaveLength(2)
    expect(gene[1]!.child_features).toHaveLength(2)
    expect(gene[0]!.child_features).toBe(gene[1]!.child_features)
  })

  it('can parse another string synchronously', () => {
    const gff3 = `
SL2.40%25ch01	IT%25AG eugene	g%25e;ne	80999140	81004317	.	+	.	Alias=Solyc01g098840;ID=gene:Solyc01g098840.2;Name=Solyc01g098840.2;from_BOGAS=1;length=5178
`

    const result = parseStringSync(gff3)
    expect(result).toMatchSnapshot()
  })
})

describe('JBrowse format parser', () => {
  it('parses with 0-based start, numeric strand, refName, and lowercased attributes', () => {
    const gff3 = `ctg123\ttest\tgene\t1000\t9000\t0.5\t+\t.\tID=gene00001;Name=TestGene`

    const result = parseStringSyncJBrowse(gff3)
    expect(result.length).toBe(1)
    const feature = result[0]
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
    const gff3 = `chr1\t.\tgene\t100\t200\t.\t-\t.\tID=g1
chr1\t.\tgene\t300\t400\t.\t.\t.\tID=g2`

    const result = parseStringSyncJBrowse(gff3)
    expect(result[0].strand).toBe(-1)
    expect(result[1].strand).toBe(0)
  })

  it('parses phase as number', () => {
    const gff3 = `chr1\t.\tCDS\t100\t200\t.\t+\t2\tID=cds1`

    const result = parseStringSyncJBrowse(gff3)
    expect(result[0].phase).toBe(2)
  })

  it('builds subfeatures from parent/child relationships', () => {
    const gff3 = `ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001
ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00001;Parent=gene00001
ctg123\t.\texon\t1050\t1500\t.\t+\t.\tID=exon1;Parent=mRNA00001`

    const result = parseStringSyncJBrowse(gff3)
    expect(result.length).toBe(1)
    const gene = result[0]
    expect(gene.id).toBe('gene00001')
    expect(gene.subfeatures.length).toBe(1)
    const mrna = gene.subfeatures[0]
    expect(mrna.id).toBe('mRNA00001')
    expect(mrna.subfeatures.length).toBe(1)
    expect(mrna.subfeatures[0].id).toBe('exon1')
  })

  it('attaches every segment of a multi-location child to its parent', () => {
    const gff3 = `ctgA\t.\tgene\t1\t1000\t.\t+\t.\tID=gene1
ctgA\t.\tmRNA\t1\t1000\t.\t+\t.\tID=mRNA1;Parent=gene1
ctgA\t.\tCDS\t1\t100\t.\t+\t0\tID=cds1;Parent=mRNA1
ctgA\t.\tCDS\t200\t300\t.\t+\t0\tID=cds1;Parent=mRNA1
ctgA\t.\tCDS\t400\t500\t.\t+\t0\tID=cds1;Parent=mRNA1`

    const result = parseStringSyncJBrowse(gff3)
    const mrna = result[0].subfeatures[0]
    const cds = mrna.subfeatures.filter(f => f.type === 'CDS')
    expect(cds.length).toBe(3)
    expect(cds.map(f => f.start)).toEqual([0, 199, 399])
  })

  it('keeps multi-value attributes as arrays', () => {
    const gff3 = `chr1\t.\tgene\t100\t200\t.\t+\t.\tID=g1;Dbxref=GO:123,GO:456`

    const result = parseStringSyncJBrowse(gff3)
    expect(result[0].dbxref).toEqual(['GO:123', 'GO:456'])
  })

  it('adds suffix to attribute names that conflict with default fields', () => {
    const gff3 = `chr1\t.\tgene\t100\t200\t.\t+\t.\tID=g1;Start=custom_start;Type=custom_type`

    const result = parseStringSyncJBrowse(gff3)
    expect(result[0].start).toBe(99) // actual start field
    expect(result[0].start2).toBe('custom_start') // attribute with suffix
    expect(result[0].type).toBe('gene') // actual type field
    expect(result[0].type2).toBe('custom_type') // attribute with suffix
  })

  it('parseRecordsJBrowse includes _lineHash', () => {
    const records = [
      {
        line: 'ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001',
        lineHash: 'offset123',
        start: 1000,
        end: 9000,
        hasEscapes: false,
      },
      {
        line: 'ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00001;Parent=gene00001',
        lineHash: 456,
        start: 1050,
        end: 9000,
        hasEscapes: false,
      },
    ]
    const result = parseRecordsJBrowse(records)
    expect(result.length).toBe(1)
    expect(result[0]._lineHash).toBe('offset123')
    expect(result[0].subfeatures[0]._lineHash).toBe('456')
  })

  it('takes the first value when ID is multi-valued', () => {
    const gff3 = `ctg\t.\tgene\t1\t10\t.\t+\t.\tID=a,b
ctg\t.\tmRNA\t1\t5\t.\t+\t.\tID=m;Parent=a`

    const result = parseStringSyncJBrowse(gff3)
    expect(result.length).toBe(1)
    expect(result[0]!.id).toEqual(['a', 'b'])
    // Parent=a should match the first ID value
    expect(result[0]!.subfeatures.length).toBe(1)
    expect(result[0]!.subfeatures[0]!.id).toBe('m')
  })

  it('handles escaped characters', () => {
    const gff3 = `SL2.40%25ch01\tIT%25AG\tgene\t100\t200\t.\t+\t.\tID=gene%3B1;Name=Test%20Gene`

    const result = parseStringSyncJBrowse(gff3)
    expect(result[0].refName).toBe('SL2.40%ch01')
    expect(result[0].source).toBe('IT%AG')
    expect(result[0].id).toBe('gene;1')
    expect(result[0].name).toBe('Test Gene')
  })
})
