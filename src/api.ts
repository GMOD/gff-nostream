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

/**
 * Extract the GFF3 feature type (column 3) from a raw line without a full
 * split. Returns '' for a line with fewer than two tabs, where there is no
 * third column to read.
 */
export function extractType(line: string): string {
  const t1 = line.indexOf('\t')
  const t2 = t1 === -1 ? -1 : line.indexOf('\t', t1 + 1)
  if (t2 === -1) {
    return ''
  } else {
    const t3 = line.indexOf('\t', t2 + 1)
    return line.slice(t2 + 1, t3 === -1 ? line.length : t3)
  }
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
 * - `top-level`: no Parent, collect it now
 * - `attached`: nested under at least one parent seen so far
 * - `orphaned`: every Parent is still unseen; collect it at the end unless a
 *   later line turns out to define one of them
 */
type LinkStatus = 'top-level' | 'attached' | 'orphaned'

/**
 * Register a feature's ID and attach it to its parent(s), building the
 * subfeature tree in `byId`/`orphans`.
 */
function linkFeature(
  feature: GffFeature,
  byId: Map<string, GffFeature>,
  orphans: Map<string, GffFeature[]>,
): LinkStatus {
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

  let attached = false
  for (const parentId of parents) {
    const parentFeature = byId.get(parentId)
    if (parentFeature) {
      parentFeature.subfeatures.push(feature)
      attached = true
    } else {
      appendOrphan(orphans, parentId, feature)
    }
  }

  // Every line of a top-level discontinuous feature (e.g. cDNA_match spanning
  // several segments under one shared ID, with no Parent) is its own top-level
  // item, so this is independent of whether the id was just registered.
  return parents.length === 0 ? 'top-level' : attached ? 'attached' : 'orphaned'
}

/**
 * True when none of a feature's Parent ids were ever defined in the input, so
 * it was never nested anywhere. Registering an id adopts everything waiting on
 * it, so presence in `byId` after the full pass means the feature was attached.
 */
function isUnparented(feature: GffFeature, byId: Map<string, GffFeature>) {
  return !toStringArray(feature.parent).some(parentId => byId.has(parentId))
}

/**
 * Synchronously parse a string containing GFF3 and return an array of the
 * parsed features. Comments, directives, and `##FASTA` sections are ignored.
 * Features whose Parent is never defined in the input (common when parsing a
 * slice of a file, e.g. a tabix region query) are returned at the end as
 * top-level items rather than dropped.
 *
 * @param str - GFF3 string
 * @returns array of parsed features
 */
export function parseStringSync(str: string): GffFeature[] {
  const items: GffFeature[] = []
  const byId = new Map<string, GffFeature>()
  const orphans = new Map<string, GffFeature[]>()
  const pending: GffFeature[] = []

  for (const line of str.split(/\r?\n/)) {
    if (line.startsWith('##FASTA') || line.startsWith('>')) {
      break
    }
    if (line.length !== 0 && !line.startsWith('#')) {
      const feature = parseFeature(line)
      const status = linkFeature(feature, byId, orphans)
      if (status === 'top-level') {
        items.push(feature)
      } else if (status === 'orphaned') {
        pending.push(feature)
      }
    }
  }

  for (const feature of pending) {
    if (isUnparented(feature, byId)) {
      items.push(feature)
    }
  }

  return items
}

/**
 * Parse an array of records wrapping raw GFF3 lines, resolving parent/child
 * relationships into `subfeatures`. Returns each top-level feature paired with
 * the record it came from, so callers can attach their own identity (e.g. a
 * byte offset) without the parser stamping anything onto the feature.
 * Features whose Parent is never defined in `records` (common for a tabix
 * region query that cuts off the parent line) are returned at the end as
 * top-level items rather than dropped.
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
  const pending: ParsedRecord<R>[] = []

  for (const record of records) {
    const feature = parseFeature(record.line)
    const status = linkFeature(feature, byId, orphans)
    if (status === 'top-level') {
      items.push({ feature, record })
    } else if (status === 'orphaned') {
      pending.push({ feature, record })
    }
  }

  for (const parsed of pending) {
    if (isUnparented(parsed.feature, byId)) {
      items.push(parsed)
    }
  }

  return items
}

export type { GffFeature } from './util.ts'
