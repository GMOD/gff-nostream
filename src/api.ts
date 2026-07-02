import { parseFeature } from './util.ts'

import type { GffFeature } from './util.ts'

export interface LineRecord {
  /** Raw GFF3 feature line */
  line: string
}

/**
 * A top-level parsed feature paired with the input record it came from. The
 * parser stamps no identity onto the feature itself; callers that need a stable
 * per-feature id (e.g. from a tabix byte offset) read it off their own `record`.
 */
export interface ParsedRecord<R extends LineRecord = LineRecord> {
  feature: GffFeature
  record: R
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
 * Register a feature's ID and attach it to its parent(s), building the
 * subfeature tree in `byId`/`orphans`. Returns true when the feature is
 * top-level (has no Parent) and the caller should collect it.
 */
function linkFeature(
  feature: GffFeature,
  byId: Map<string, GffFeature>,
  orphans: Map<string, GffFeature[]>,
): boolean {
  const id = firstString(feature.id)
  const parents = toStringArray(feature.parent)

  // Register the id only the first time it is seen. Continuation lines
  // (multi-location features such as a CDS spanning several segments share one
  // ID across lines) skip registration but must still be attached to their
  // parent below, so this is independent of the parent handling.
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

  // Every line of a top-level discontinuous feature (e.g. cDNA_match spanning
  // several segments under one shared ID, with no Parent) is its own top-level
  // item, so this is independent of whether the id was just registered.
  return parents.length === 0
}

/**
 * Synchronously parse a string containing GFF3 and return an array of the
 * parsed features. Comments, directives, and `##FASTA` sections are ignored.
 *
 * @param str - GFF3 string
 * @returns array of parsed features
 */
export function parseStringSync(str: string): GffFeature[] {
  const items: GffFeature[] = []
  const byId = new Map<string, GffFeature>()
  const orphans = new Map<string, GffFeature[]>()

  for (const line of str.split(/\r?\n/)) {
    if (line.startsWith('##FASTA') || line.startsWith('>')) {
      break
    }
    if (line.length !== 0 && !line.startsWith('#')) {
      const feature = parseFeature(line)
      if (linkFeature(feature, byId, orphans)) {
        items.push(feature)
      }
    }
  }

  return items
}

/**
 * Parse an array of records wrapping raw GFF3 lines, resolving parent/child
 * relationships into `subfeatures`. Returns each top-level feature paired with
 * the record it came from, so callers can attach their own identity (e.g. a
 * byte offset) without the parser stamping anything onto the feature.
 *
 * @param records - Array of records, each carrying a raw GFF3 `line`
 * @returns top-level features, each paired with its originating record
 */
export function parseRecords<R extends LineRecord>(
  records: readonly R[],
): ParsedRecord<R>[] {
  const items: ParsedRecord<R>[] = []
  const byId = new Map<string, GffFeature>()
  const orphans = new Map<string, GffFeature[]>()

  for (const record of records) {
    const feature = parseFeature(record.line)
    if (linkFeature(feature, byId, orphans)) {
      items.push({ feature, record })
    }
  }

  return items
}

export type { GffFeature } from './util.ts'
