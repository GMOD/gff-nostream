import { parseFeature } from './util.ts'

import type { GffFeature } from './util.ts'

interface ParseInput {
  line: string
  lineHash?: string | number
  hasEscapes: boolean
}

export interface LineRecord extends ParseInput {
  /** Genomic start coordinate from the tabix index (1-based) */
  start: number
  /** Genomic end coordinate from the tabix index */
  end: number
  /** GFF3 feature type (column 3) */
  type: string
}

/** Extract the GFF3 feature type (column 3) from a raw line without a full split. */
export function extractType(line: string): string {
  const t1 = line.indexOf('\t')
  const t2 = line.indexOf('\t', t1 + 1)
  const t3 = line.indexOf('\t', t2 + 1)
  return line.slice(t2 + 1, t3)
}

/** Append a value to the array stored under key, creating the array if absent. */
function appendOrphan<T>(orphans: Map<string, T[]>, key: string, value: T) {
  const arr = orphans.get(key)
  if (arr) {
    arr.push(value)
  } else {
    orphans.set(key, [value])
  }
}

/**
 * The parser collapses single-element attribute arrays to scalars, so a raw
 * ID/Parent value can be a string, a string array, or absent. These coerce
 * those `unknown` values without typecasts.
 */
function firstString(value: unknown): string | undefined {
  const v: unknown = Array.isArray(value) ? value[0] : value
  return typeof v === 'string' ? v : undefined
}

function toStringArray(value: unknown): string[] {
  if (Array.isArray(value)) {
    return value.filter((v): v is string => typeof v === 'string')
  }
  return typeof value === 'string' ? [value] : []
}

/**
 * Synchronously parse a string containing GFF3 and return an array of the
 * parsed features. Comments, directives, and `##FASTA` sections are ignored.
 *
 * @param str - GFF3 string
 * @returns array of parsed features
 */
export function parseStringSync(str: string): GffFeature[] {
  return parseRecords(stringToRecords(str))
}

function stringToRecords(str: string) {
  const lines = str.split(/\r?\n/)
  const records: ParseInput[] = []
  for (const line of lines) {
    if (line.startsWith('##FASTA') || line.startsWith('>')) {
      break
    }
    if (line.length === 0 || line.startsWith('#')) {
      continue
    }
    records.push({
      line,
      hasEscapes: line.includes('%'),
    })
  }
  return records
}

/**
 * Parse an array of LineRecord objects containing raw GFF3 lines.
 * Supports parent/child relationships via subfeatures.
 *
 * @param records - Array of LineRecord objects with raw line and metadata
 * @returns array of parsed features
 */
export function parseRecords(records: ParseInput[]): GffFeature[] {
  const items: GffFeature[] = []
  const byId = new Map<string, GffFeature>()
  const orphans = new Map<string, GffFeature[]>()

  for (const record of records) {
    const feature = parseFeature(record.line, record.hasEscapes)

    if (record.lineHash !== undefined) {
      feature._lineHash = String(record.lineHash)
    }

    const id = firstString(feature.id)
    const parents = toStringArray(feature.parent)

    // A parentless line is a top-level item. Every line of a top-level
    // discontinuous feature (e.g. cDNA_match/EST_match spanning several
    // segments under one shared ID, with no Parent) is its own top-level
    // item, so push regardless of whether the id is already registered.
    if (parents.length === 0) {
      items.push(feature)
    }

    // Register the id only the first time it is seen. Continuation lines
    // (multi-location features such as a CDS spanning several segments share
    // one ID across lines) skip registration but must still be attached to
    // their parent below, so this is independent of the parent handling.
    if (id && !byId.has(id)) {
      byId.set(id, feature)
      const waiting = orphans.get(id)
      if (waiting) {
        for (const w of waiting) {
          feature.subfeatures.push(w)
        }
        orphans.delete(id)
      }
    }

    for (const parentId of parents) {
      const parentFeature = byId.get(parentId)
      if (parentFeature) {
        parentFeature.subfeatures.push(feature)
      } else {
        appendOrphan(orphans, parentId, feature)
      }
    }
  }

  return items
}

export type { GffFeature } from './util.ts'
