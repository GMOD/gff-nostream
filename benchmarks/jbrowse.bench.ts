import fs from 'fs'
import { bench, describe } from 'vitest'

import { parseStringSync } from '../src/api.ts'
import { parseFeature, parseAttributes, unescape } from '../src/util.ts'
import type { GffFeature } from '../src/util.ts'

// Local record shape for the optimization-variant harness below; the shipped
// `parseRecords` now derives the unescape flag itself from the line.
interface LineRecord {
  line: string
  lineHash?: string | number
  start: number
  end: number
  hasEscapes: boolean
}

// The only optimization that showed improvement: avoid toLowerCase for common GFF3 attribute names
const JBROWSE_DEFAULT_FIELDS = new Set([
  'start',
  'end',
  'seq_id',
  'score',
  'type',
  'source',
  'phase',
  'strand',
])

// Common GFF3 attribute names pre-lowercased (these are defined in the GFF3 spec)
const COMMON_ATTRS: Record<string, string> = {
  ID: 'id',
  Name: 'name',
  Parent: 'parent',
  Note: 'note',
  Dbxref: 'dbxref',
  Ontology_term: 'ontology_term',
  Is_circular: 'is_circular',
  Alias: 'alias',
  Target: 'target',
  Gap: 'gap',
  Derives_from: 'derives_from',
  // Add lowercase versions too for case-insensitive matching
  id: 'id',
  name: 'name',
  parent: 'parent',
  note: 'note',
  dbxref: 'dbxref',
  alias: 'alias',
  target: 'target',
  gap: 'gap',
}

function parseAttributesJBrowseOptimized(
  attrString: string,
  result: Record<string, unknown>,
) {
  if (attrString.length === 0 || attrString === '.') {
    return
  }

  let len = attrString.length
  if (attrString[len - 1] === '\n') {
    len = attrString[len - 2] === '\r' ? len - 2 : len - 1
    attrString = attrString.slice(0, len)
  }

  let start = 0
  while (start < len) {
    let semiIdx = attrString.indexOf(';', start)
    if (semiIdx === -1) {
      semiIdx = len
    }

    if (semiIdx > start) {
      const eqIdx = attrString.indexOf('=', start)
      if (eqIdx !== -1 && eqIdx < semiIdx && eqIdx + 1 < semiIdx) {
        const tag = attrString.slice(start, eqIdx)
        if (tag === '_lineHash') {
          start = semiIdx + 1
          continue
        }

        // Use pre-computed lowercase for common attributes, fall back to toLowerCase
        let key = COMMON_ATTRS[tag]
        if (!key) {
          key = tag.toLowerCase()
          if (JBROWSE_DEFAULT_FIELDS.has(key)) {
            key += '2'
          }
        }

        const values: string[] = []
        let valStart = eqIdx + 1
        while (valStart < semiIdx) {
          let commaIdx = attrString.indexOf(',', valStart)
          if (commaIdx === -1 || commaIdx > semiIdx) {
            commaIdx = semiIdx
          }
          if (commaIdx > valStart) {
            const val = attrString.slice(valStart, commaIdx)
            values.push(unescape(val))
          }
          valStart = commaIdx + 1
        }

        result[key] = values.length === 1 ? values[0] : values
      }
    }
    start = semiIdx + 1
  }
}

function parseAttributesJBrowseNoUnescapeOptimized(
  attrString: string,
  result: Record<string, unknown>,
) {
  if (attrString.length === 0 || attrString === '.') {
    return
  }

  let len = attrString.length
  if (attrString[len - 1] === '\n') {
    len = attrString[len - 2] === '\r' ? len - 2 : len - 1
    attrString = attrString.slice(0, len)
  }

  let start = 0
  while (start < len) {
    let semiIdx = attrString.indexOf(';', start)
    if (semiIdx === -1) {
      semiIdx = len
    }

    if (semiIdx > start) {
      const eqIdx = attrString.indexOf('=', start)
      if (eqIdx !== -1 && eqIdx < semiIdx && eqIdx + 1 < semiIdx) {
        const tag = attrString.slice(start, eqIdx)
        if (tag === '_lineHash') {
          start = semiIdx + 1
          continue
        }

        let key = COMMON_ATTRS[tag]
        if (!key) {
          key = tag.toLowerCase()
          if (JBROWSE_DEFAULT_FIELDS.has(key)) {
            key += '2'
          }
        }

        const values: string[] = []
        let valStart = eqIdx + 1
        while (valStart < semiIdx) {
          let commaIdx = attrString.indexOf(',', valStart)
          if (commaIdx === -1 || commaIdx > semiIdx) {
            commaIdx = semiIdx
          }
          if (commaIdx > valStart) {
            values.push(attrString.slice(valStart, commaIdx))
          }
          valStart = commaIdx + 1
        }

        result[key] = values.length === 1 ? values[0] : values
      }
    }
    start = semiIdx + 1
  }
}

function parseFeatureJBrowseOptimized(line: string): GffFeature {
  const f = line.split('\t')
  const seq_id = f[0]!
  const source = f[1]!
  const type = f[2]!
  const startStr = f[3]!
  const endStr = f[4]!
  const scoreStr = f[5]!
  const strand = f[6]!
  const phase = f[7]!
  const attrString = f[8]!

  const strandVal =
    strand === '+' ? 1 : strand === '-' ? -1 : strand === '.' ? 0 : undefined

  const result: GffFeature = {
    refName: seq_id.length === 0 || seq_id === '.' ? '' : unescape(seq_id),
    source: source.length === 0 || source === '.' ? null : unescape(source),
    type: type.length === 0 || type === '.' ? null : unescape(type),
    start: startStr.length === 0 || startStr === '.' ? 0 : +startStr - 1,
    end: endStr.length === 0 || endStr === '.' ? 0 : +endStr,
    score: scoreStr.length === 0 || scoreStr === '.' ? undefined : +scoreStr,
    strand: strandVal,
    phase: phase.length === 0 || phase === '.' ? undefined : +phase,
    subfeatures: [],
  }

  parseAttributesJBrowseOptimized(attrString, result)
  return result
}

function parseFeatureJBrowseNoUnescapeOptimized(line: string): GffFeature {
  const f = line.split('\t')
  const seq_id = f[0]!
  const source = f[1]!
  const type = f[2]!
  const startStr = f[3]!
  const endStr = f[4]!
  const scoreStr = f[5]!
  const strand = f[6]!
  const phase = f[7]!
  const attrString = f[8]!

  const strandVal =
    strand === '+' ? 1 : strand === '-' ? -1 : strand === '.' ? 0 : undefined

  const result: GffFeature = {
    refName: seq_id.length === 0 || seq_id === '.' ? '' : seq_id,
    source: source.length === 0 || source === '.' ? null : source,
    type: type.length === 0 || type === '.' ? null : type,
    start: startStr.length === 0 || startStr === '.' ? 0 : +startStr - 1,
    end: endStr.length === 0 || endStr === '.' ? 0 : +endStr,
    score: scoreStr.length === 0 || scoreStr === '.' ? undefined : +scoreStr,
    strand: strandVal,
    phase: phase.length === 0 || phase === '.' ? undefined : +phase,
    subfeatures: [],
  }

  parseAttributesJBrowseNoUnescapeOptimized(attrString, result)
  return result
}

function parseRecordsJBrowseOptimized(records: LineRecord[]): GffFeature[] {
  const items: GffFeature[] = []
  const byId = new Map<string, GffFeature>()
  const orphans = new Map<string, GffFeature[]>()

  for (let i = 0; i < records.length; i++) {
    const record = records[i]!
    const feature = record.hasEscapes
      ? parseFeatureJBrowseOptimized(record.line)
      : parseFeatureJBrowseNoUnescapeOptimized(record.line)

    if (record.lineHash !== undefined) {
      feature._lineHash = String(record.lineHash)
    }

    const id = feature.id as string | undefined
    const parent = feature.parent as string | string[] | undefined

    if (!id && !parent) {
      items.push(feature)
      continue
    }

    if (id) {
      const existing = byId.get(id)
      if (!existing) {
        if (!parent) {
          items.push(feature)
        }
        byId.set(id, feature)
        const waiting = orphans.get(id)
        if (waiting) {
          for (let j = 0; j < waiting.length; j++) {
            feature.subfeatures.push(waiting[j]!)
          }
          orphans.delete(id)
        }
      }
    }

    if (parent) {
      const parents = Array.isArray(parent) ? parent : [parent]
      for (let j = 0; j < parents.length; j++) {
        const parentId = parents[j]!
        const parentFeature = byId.get(parentId)
        if (parentFeature) {
          parentFeature.subfeatures.push(feature)
        } else {
          let arr = orphans.get(parentId)
          if (!arr) {
            arr = []
            orphans.set(parentId, arr)
          }
          arr.push(feature)
        }
      }
    }
  }

  return items
}

function stringToRecords(str: string) {
  const lines = str.split(/\r?\n/)
  const records: LineRecord[] = []
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i]!
    if (line.length === 0 || line[0] === '#') {
      if (line.startsWith('##FASTA')) {
        break
      }
      continue
    }
    if (line[0] === '>') {
      break
    }
    records.push({
      line,
      start: 0,
      end: 0,
      hasEscapes: line.includes('%'),
    })
  }
  return records
}

function parseStringSyncJBrowseOptimized(str: string): GffFeature[] {
  return parseRecordsJBrowseOptimized(stringToRecords(str))
}

// Benchmark attributes parsing with common attributes (isolated test)
describe('parseAttributesJBrowse - common attrs', () => {
  const attrs = 'ID=AT1G01010;Name=NAC001;Parent=gene1;Note=something'

  bench('current', () => {
    const result: Record<string, unknown> = {}
    parseAttributes(attrs, result, true)
  })

  bench('common attr lookup', () => {
    const result: Record<string, unknown> = {}
    parseAttributesJBrowseOptimized(attrs, result)
  })
})

// Benchmark attributes with uncommon attrs
describe('parseAttributesJBrowse - uncommon attrs', () => {
  const attrs =
    'gene_id=ENSG00000223972;transcript_id=ENST00000456328;biotype=lncRNA'

  bench('current', () => {
    const result: Record<string, unknown> = {}
    parseAttributes(attrs, result, true)
  })

  bench('common attr lookup', () => {
    const result: Record<string, unknown> = {}
    parseAttributesJBrowseOptimized(attrs, result)
  })
})

// Benchmark single feature parsing
describe('parseFeatureJBrowse - single line', () => {
  const line =
    'chr1\tAraport11\tgene\t3631\t5899\t.\t+\t.\tID=AT1G01010;Name=NAC001;Note=NAC domain containing protein 1'

  bench('current', () => {
    parseFeature(line)
  })

  bench('optimized', () => {
    parseFeatureJBrowseOptimized(line)
  })
})

// Benchmark full file parsing - large file
describe('parseStringSyncJBrowse - large file (tair10_chr1)', () => {
  const data = fs.readFileSync('test/data/tair10_chr1.gff', 'utf8')

  bench(
    'current',
    () => {
      parseStringSync(data)
    },
    { iterations: 50, warmupIterations: 5 },
  )

  bench(
    'optimized (common attr lookup)',
    () => {
      parseStringSyncJBrowseOptimized(data)
    },
    { iterations: 50, warmupIterations: 5 },
  )
})

// Benchmark full file parsing - medium file
describe('parseStringSyncJBrowse - medium file (au9)', () => {
  const data = fs.readFileSync('test/data/au9_scaffold_subset.gff3', 'utf8')

  bench('current', () => {
    parseStringSync(data)
  })

  bench('optimized (common attr lookup)', () => {
    parseStringSyncJBrowseOptimized(data)
  })
})

// Benchmark messy protein domains (lots of escapes and uncommon attrs)
describe('parseStringSyncJBrowse - messy file', () => {
  const data = fs.readFileSync('test/data/messy_protein_domains.gff3', 'utf8')

  bench('current', () => {
    parseStringSync(data)
  })

  bench('optimized (common attr lookup)', () => {
    parseStringSyncJBrowseOptimized(data)
  })
})
