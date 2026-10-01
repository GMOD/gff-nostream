import fs from 'node:fs'

import { describe, expect, it } from 'vitest'

import {
  extractType,
  getAttribute,
  getLinkAttributes,
  hasIdAttribute,
  parseLines,
  parseLinesLazy,
  parseRecords,
  parseRecordsLazy,
  parseStringSync,
} from '../src/index.ts'
import { unescape } from '../src/util.ts'

import type { GffFeature } from '../src/index.ts'

const gagPolLines = [
  'NC_001802.1\tRefSeq\tgene\t336\t4642\t.\t+\t.\tID=gene-HIV1gp1;Name=gag-pol',
  'NC_001802.1\tRefSeq\tCDS\t336\t1637\t.\t+\t0\tID=cds-NP_057849.4;Parent=gene-HIV1gp1',
  'NC_001802.1\tRefSeq\tCDS\t1637\t4642\t.\t+\t0\tID=cds-NP_057849.4;Parent=gene-HIV1gp1',
  'NC_001802.1\tRefSeq\tmature_protein_region_of_CDS\t1799\t2095\t.\t+\t.\tID=id-NP_057849.4:489..587;Parent=cds-NP_057849.4;product=protease',
  'NC_001802.1\tRefSeq\tmature_protein_region_of_CDS\t3776\t4639\t.\t+\t.\tID=id-NP_057849.4:1148..1435;Parent=cds-NP_057849.4;product=integrase',
]

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

  it('folds a top-level discontinuous feature into one spanning its segments', () => {
    const result = parseStringSync(
      `ctgA\t.\tcDNA_match\t1050\t1500\t5.8e-42\t+\t.\tID=match1;Target=NM_1 1 451 +
ctgA\t.\tcDNA_match\t5000\t5500\t8.1e-43\t+\t.\tID=match1;Target=NM_1 452 952 +
ctgA\t.\tcDNA_match\t7000\t9000\t1.4e-40\t+\t.\tID=match1;Target=NM_1 953 2953 +`,
    )
    expect(result.length).toBe(1)
    const [aln] = result
    expect([aln!.start, aln!.end, aln!.type]).toEqual([
      1049,
      9000,
      'cDNA_match',
    ])
    expect(aln!.subfeatures.map(f => [f.start, f.end, f.type])).toEqual([
      [1049, 1500, 'cDNA_match'],
      [4999, 5500, 'cDNA_match'],
      [6999, 9000, 'cDNA_match'],
    ])
    expect(aln!.subfeatures.map(f => f.target)).toEqual([
      'NM_1 1 451 +',
      'NM_1 452 952 +',
      'NM_1 953 2953 +',
    ])
    expect(aln!.subfeatures.every(f => f.subfeatures.length === 0)).toBe(true)
  })

  it('gives a single-line top-level feature no segment of itself', () => {
    const result = parseStringSync(
      `ctgA\t.\tmatch\t1050\t1500\t.\t+\t.\tID=m1
ctgA\t.\tmatch\t5000\t5500\t.\t+\t.\tID=m2`,
    )
    expect(result.map(f => f.subfeatures.length)).toEqual([0, 0])
  })

  // a match carrying match_part children is already a container, so a second
  // match line under its ID is not a segment of it
  it('does not fold into a line that already has children', () => {
    const result = parseStringSync(
      `ctgA\t.\tmatch\t1000\t1500\t.\t+\t.\tID=m1
ctgA\t.\tmatch_part\t1000\t1100\t.\t+\t.\tParent=m1
ctgA\t.\tmatch\t5000\t5500\t.\t+\t.\tID=m1`,
    )
    expect(result.map(f => [f.start, f.subfeatures.map(s => s.type)])).toEqual([
      [999, ['match_part']],
      [4999, []],
    ])
  })

  it('parseRecords pairs a folded feature with its first record', () => {
    const records = [
      'ctgA\t.\tcDNA_match\t1050\t1500\t.\t+\t.\tID=match1',
      'ctgA\t.\tcDNA_match\t5000\t5500\t.\t+\t.\tID=match1',
    ].map((line, offset) => ({ line, offset }))
    const result = parseRecords(records)
    expect(result.map(r => r.record.offset)).toEqual([0])
    expect(result[0]!.feature.subfeatures.length).toBe(2)
  })

  // ngsutils names every gene by its symbol, so two genes share an ID; the
  // second is not a segment of the first, and its children still attach to
  // whichever line registered the ID, as before.
  it('leaves a repeated ID alone once the first line has children', () => {
    const result = parseStringSync(
      `chr1\t.\tgene\t100\t500\t.\t+\t.\tID=APITD1
chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=NM_1;Parent=APITD1
chr1\t.\tgene\t150\t500\t.\t+\t.\tID=APITD1
chr1\t.\tmRNA\t150\t500\t.\t+\t.\tID=NM_2;Parent=APITD1`,
    )
    expect(result.map(f => [f.start, f.subfeatures.length])).toEqual([
      [99, 2],
      [149, 0],
    ])
  })

  it('leaves a repeated ID alone across types', () => {
    const result = parseStringSync(
      `chr1\t.\tgene\t100\t500\t.\t+\t.\tID=x
chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=x`,
    )
    expect(result.map(f => [f.type, f.subfeatures.length])).toEqual([
      ['gene', 0],
      ['mRNA', 0],
    ])
  })

  // The four below guard the same code path as the two above, in the shapes
  // that path is easiest to get wrong. It regressed once: 3.0.6 through 3.0.9
  // registered a shared ID per line and dropped the continuation lines, so a
  // GENCODE transcript came out with one CDS of its four while all its exons
  // survived — exons carry unique IDs, CDS segments share one, which is what
  // made it look like a rendering quirk rather than a parse bug. Consumers
  // translate from the full CDS set, so the visible symptom was a silently
  // shortened protein (NRAS 190aa -> 40aa) with nothing reporting an error.
  // It reached no release: jbrowse's lockfile sat on a bad version for thirteen
  // days on main and never in a tag. Keep these even if they look redundant —
  // the two tests above pass under several ways of getting this wrong.

  it('attaches continuation lines that arrive before their parent', () => {
    const result = parseStringSync(
      `ctgA\t.\tCDS\t1\t100\t.\t+\t0\tID=cds1;Parent=mRNA1
ctgA\t.\tCDS\t200\t300\t.\t+\t0\tID=cds1;Parent=mRNA1
ctgA\t.\tCDS\t400\t500\t.\t+\t0\tID=cds1;Parent=mRNA1
ctgA\t.\tmRNA\t1\t500\t.\t+\t.\tID=mRNA1`,
    )
    const cds = result[0]!.subfeatures.filter(f => f.type === 'CDS')
    expect(cds.map(f => f.start)).toEqual([0, 199, 399])
  })

  it('attaches continuation lines interleaved with uniquely-identified siblings', () => {
    const result = parseStringSync(
      `ctgA\t.\tmRNA\t1\t500\t.\t+\t.\tID=mRNA1
ctgA\t.\tCDS\t1\t100\t.\t+\t0\tID=cds1;Parent=mRNA1
ctgA\t.\texon\t1\t100\t.\t+\t.\tID=exon1;Parent=mRNA1
ctgA\t.\tCDS\t200\t300\t.\t+\t0\tID=cds1;Parent=mRNA1
ctgA\t.\texon\t200\t300\t.\t+\t.\tID=exon2;Parent=mRNA1
ctgA\t.\tCDS\t400\t500\t.\t+\t0\tID=cds1;Parent=mRNA1`,
    )
    const subs = result[0]!.subfeatures
    expect(subs.filter(f => f.type === 'CDS').map(f => f.start)).toEqual([
      0, 199, 399,
    ])
    expect(subs.filter(f => f.type === 'exon').length).toBe(2)
  })

  it('gives every parent the full set of a multi-parent shared ID', () => {
    const result = parseStringSync(
      `ctgA\t.\tmRNA\t1\t500\t.\t+\t.\tID=mRNA1
ctgA\t.\tmRNA\t1\t500\t.\t+\t.\tID=mRNA2
ctgA\t.\tCDS\t1\t100\t.\t+\t0\tID=cds1;Parent=mRNA1,mRNA2
ctgA\t.\tCDS\t200\t300\t.\t+\t0\tID=cds1;Parent=mRNA1,mRNA2`,
    )
    for (const mrna of result) {
      expect(mrna.subfeatures.filter(f => f.type === 'CDS').length).toBe(2)
    }
  })

  it('keeps two transcripts’ shared-ID CDS sets apart', () => {
    const result = parseStringSync(
      `ctgA\t.\tmRNA\t1\t500\t.\t+\t.\tID=mRNA1
ctgA\t.\tCDS\t1\t100\t.\t+\t0\tID=cds1;Parent=mRNA1
ctgA\t.\tCDS\t200\t300\t.\t+\t0\tID=cds1;Parent=mRNA1
ctgA\t.\tmRNA\t1\t500\t.\t+\t.\tID=mRNA2
ctgA\t.\tCDS\t1\t100\t.\t+\t0\tID=cds2;Parent=mRNA2
ctgA\t.\tCDS\t200\t300\t.\t+\t0\tID=cds2;Parent=mRNA2`,
    )
    expect(result.map(m => m.subfeatures.length)).toEqual([2, 2])
  })

  describe('a parented ID over several lines', () => {
    const sars = fs.readFileSync('test/data/sars_cov2_NC_045512.2.gff3', 'utf8')
    const hiv = fs.readFileSync('test/data/hiv1_NC_001802.1.gff3', 'utf8')

    const span = (f: GffFeature) => [f.type, f.start, f.end]
    const byGene = (features: GffFeature[], id: string) =>
      features.find(f => f.id === id)!.subfeatures

    it('folds SARS-CoV-2 ORF1ab, a -1 frameshift with mature peptides', () => {
      const [pp1ab, pp1a, ...rest] = byGene(
        parseStringSync(sars),
        'gene-GU280_gp01',
      )
      expect(rest).toEqual([])
      expect([pp1ab!.id, ...span(pp1ab!)]).toEqual([
        'cds-YP_009724389.1',
        'CDS',
        265,
        21555,
      ])
      expect([pp1a!.id, ...span(pp1a!)]).toEqual([
        'cds-YP_009725295.1',
        'CDS',
        265,
        13483,
      ])
      const segments = pp1ab!.subfeatures.filter(f => f.type === 'CDS')
      expect(segments.map(span)).toEqual([
        ['CDS', 265, 13468],
        ['CDS', 13467, 21555],
      ])
      expect(segments.every(f => f.subfeatures.length === 0)).toBe(true)
      const peptides = pp1ab!.subfeatures.slice(2)
      expect(peptides.length).toBe(16)
      expect(
        peptides.every(f => f.type === 'mature_protein_region_of_CDS'),
      ).toBe(true)
      expect(peptides.at(-1)!.end).toBe(21552)
      // nsp12 straddles the frameshift on two childless lines, so it stays two
      expect(
        peptides
          .filter(f => f.product === 'RNA-dependent RNA polymerase')
          .map(span),
      ).toEqual([
        ['mature_protein_region_of_CDS', 13441, 13468],
        ['mature_protein_region_of_CDS', 13467, 16236],
      ])
    })

    it('folds HIV-1 gag-pol and leaves the spliced, childless tat, rev and vpr alone', () => {
      const features = parseStringSync(hiv)
      const [gagPol, ...rest] = byGene(features, 'gene-HIV1gp1')
      expect(rest).toEqual([])
      expect(span(gagPol!)).toEqual(['CDS', 335, 4642])
      expect(gagPol!.subfeatures.map(f => f.type)).toEqual([
        'CDS',
        'CDS',
        ...Array<string>(7).fill('mature_protein_region_of_CDS'),
      ])
      expect(gagPol!.subfeatures.slice(0, 2).map(span)).toEqual([
        ['CDS', 335, 1637],
        ['CDS', 1636, 4642],
      ])
      for (const [gene, blocks] of [
        ['gene-HIV1gp4', [5104, 5320]],
        ['gene-HIV1gp5', [5376, 7924]],
        ['gene-HIV1gp6', [5515, 7924]],
      ] as const) {
        const cds = byGene(features, gene)
        expect(cds.map(f => f.start)).toEqual(blocks)
        expect(cds.every(f => f.subfeatures.length === 0)).toBe(true)
      }
    })

    it('folds whichever order the lines arrive in', () => {
      const lines = sars
        .split('\n')
        .filter(line => line.length !== 0 && !line.startsWith('#'))
      const byStart = (a: string, b: string) =>
        +a.split('\t')[3]! - +b.split('\t')[3]!
      for (const order of [[...lines].sort(byStart), [...lines].reverse()]) {
        const [pp1ab] = byGene(parseLines(order), 'gene-GU280_gp01').filter(
          f => f.id === 'cds-YP_009724389.1',
        )
        expect(span(pp1ab!)).toEqual(['CDS', 265, 21555])
        expect(pp1ab!.subfeatures.length).toBe(18)
      }
    })

    interface Tree {
      subfeatures: Tree[]
    }
    const entryPoints: [string, (lines: string[]) => Tree[]][] = [
      ['parseLines', l => parseLines(l)],
      ['parseStringSync', l => parseStringSync(l.join('\n'))],
      [
        'parseRecords',
        l => parseRecords(l.map(line => ({ line }))).map(r => r.feature),
      ],
      ['parseLinesLazy', l => parseLinesLazy(l)],
      [
        'parseRecordsLazy',
        l => parseRecordsLazy(l.map(line => ({ line }))).map(r => r.feature),
      ],
    ]

    it.each(entryPoints)('folds through %s', (_name, parse) => {
      const [gene, ...rest] = parse(gagPolLines)
      expect(rest).toEqual([])
      expect(gene!.subfeatures.map(f => f.subfeatures.length)).toEqual([4])
    })

    // a slice that cuts off the gene leaves the CDS lines unparented
    it('folds lines whose parent is absent into one top-level feature', () => {
      const records = gagPolLines.slice(1).map((line, offset) => ({
        line,
        offset,
      }))
      const result = parseRecords(records)
      expect(result.map(r => r.record.offset)).toEqual([0])
      expect(span(result[0]!.feature)).toEqual(['CDS', 335, 4642])
    })

    it('leaves an NCBI mRNA’s childless shared-ID CDS lines one per line', () => {
      const [gene] = parseStringSync(
        `NC_000001.11\tBestRefSeq\tgene\t100\t900\t.\t+\t.\tID=gene-X
NC_000001.11\tBestRefSeq\tmRNA\t100\t900\t.\t+\t.\tID=rna-NM_1;Parent=gene-X
NC_000001.11\tBestRefSeq\texon\t100\t300\t.\t+\t.\tID=exon-NM_1-1;Parent=rna-NM_1
NC_000001.11\tBestRefSeq\texon\t500\t900\t.\t+\t.\tID=exon-NM_1-2;Parent=rna-NM_1
NC_000001.11\tBestRefSeq\tCDS\t150\t300\t.\t+\t0\tID=cds-NP_1;Parent=rna-NM_1
NC_000001.11\tBestRefSeq\tCDS\t500\t800\t.\t+\t2\tID=cds-NP_1;Parent=rna-NM_1`,
      )
      const mrna = gene!.subfeatures[0]!
      expect(mrna.subfeatures.map(span)).toEqual([
        ['exon', 99, 300],
        ['exon', 499, 900],
        ['CDS', 149, 300],
        ['CDS', 499, 800],
      ])
    })

    it('leaves a further line under other parents where it is', () => {
      const result = parseStringSync(
        `ctgA\t.\tgene\t1\t500\t.\t+\t.\tID=g1
ctgA\t.\tgene\t600\t900\t.\t+\t.\tID=g2
ctgA\t.\tmRNA\t1\t500\t.\t+\t.\tID=t1;Parent=g1
ctgA\t.\tmRNA\t600\t900\t.\t+\t.\tID=t1;Parent=g2
ctgA\t.\texon\t1\t500\t.\t+\t.\tParent=t1`,
      )
      expect(result.map(g => g.subfeatures.map(span))).toEqual([
        [['mRNA', 0, 500]],
        [['mRNA', 599, 900]],
      ])
    })
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

describe('whitespace around an attribute tag', () => {
  const gene = (attrs: string) => `c\ts\tgene\t1\t100\t.\t+\t.\t${attrs}`
  const mrna = (attrs: string) => `c\ts\tmRNA\t1\t100\t.\t+\t.\t${attrs}`

  it('links a child to a parent whose ID follows "; "', () => {
    const [eager] = parseLines([gene('Name=A; ID=x'), mrna('Parent=x')])
    expect(eager?.id).toBe('x')
    expect(eager?.subfeatures).toHaveLength(1)
    const [lazy] = parseLinesLazy([gene('Name=A; ID=x'), mrna(' Parent=x')])
    expect(getLinkAttributes(lazy!).id).toBe('x')
    expect(lazy?.subfeatures).toHaveLength(1)
  })

  it('trims space before the "=" as well', () => {
    expect(parseLines([gene('ID =x')])[0]?.id).toBe('x')
    expect(getAttribute(parseLinesLazy([gene('ID =x')])[0]!, 'id')).toBe('x')
  })
})

describe('hasIdAttribute', () => {
  const line = (attrs: string) => `chr1\tRefSeq\tgene\t1\t2\t.\t+\t.\t${attrs}`

  it('is true for a record the linker would register', () => {
    expect(hasIdAttribute(line('ID=gene-A;Name=A'))).toBe(true)
    expect(hasIdAttribute(line('Name=A;ID=gene-A'))).toBe(true)
    expect(hasIdAttribute(line('Name=A; ID=gene-A'))).toBe(true)
    expect(hasIdAttribute(line('id=gene-A'))).toBe(true)
  })

  it('is false for a record nothing can reference', () => {
    expect(
      hasIdAttribute(
        'chr1\tRefSeq\tmatch\t585989\t121976459\t.\t+\t.\t' +
          'Target=chr1 585989 121976459 +;gap_count=0;pct_identity_gap=100',
      ),
    ).toBe(false)
    expect(hasIdAttribute(line('Parent=gene-A'))).toBe(false)
    expect(hasIdAttribute(line('.'))).toBe(false)
    expect(hasIdAttribute(line('ID=;Name=A'))).toBe(false)
    expect(hasIdAttribute('chr1\tRefSeq\tgene')).toBe(false)
  })

  it('reads the tag, not the text', () => {
    expect(hasIdAttribute(line('geneID=7157;Name=TP53'))).toBe(false)
    expect(hasIdAttribute(line('Note=see ID=other'))).toBe(false)
    expect(hasIdAttribute(line('ID=a\tParent=b'))).toBe(true)
    expect(hasIdAttribute(line('Name=a\tID=b'))).toBe(false)
  })
})
