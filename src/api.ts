import {
  parseFeature,
  parseFeatureJBrowse,
  parseFeatureJBrowseNoUnescape,
  parseFeatureNoUnescape,
} from './util.ts'

import type {
  GFF3Feature,
  GFF3FeatureLineWithRefs,
  JBrowseFeature,
} from './util.ts'

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
 * The JBrowse parser collapses single-element attribute arrays to scalars, so a
 * raw ID/Parent value can be a string, a string array, or absent. These coerce
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
 * parsed items.
 *
 * @param str - GFF3 string
 * @returns array of parsed features
 */
export function parseStringSync(str: string): GFF3Feature[] {
  return parseRecords(stringToRecords(str))
}

/**
 * Synchronously parse a string containing GFF3 directly into JBrowse format.
 *
 * @param str - GFF3 string
 * @returns array of JBrowse-format features
 */
export function parseStringSyncJBrowse(str: string): JBrowseFeature[] {
  return parseRecordsJBrowse(stringToRecords(str))
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
 * Supports parent/child relationships.
 *
 * @param records - Array of LineRecord objects with raw line and metadata
 * @returns array of parsed features
 */
export function parseRecords(records: ParseInput[]): GFF3Feature[] {
  const items: GFF3Feature[] = []
  const byId = new Map<string, GFF3Feature>()
  const orphans = new Map<string, GFF3Feature[]>()

  for (const record of records) {
    const parsed = record.hasEscapes
      ? parseFeature(record.line)
      : parseFeatureNoUnescape(record.line)
    const featureLine: GFF3FeatureLineWithRefs = {
      ...parsed,
      child_features: [],
      derived_features: [],
    }

    if (record.lineHash !== undefined) {
      featureLine.attributes ??= {}
      featureLine.attributes._lineHash = [String(record.lineHash)]
    }

    const attrs = featureLine.attributes
    const ids = attrs?.ID
    const parents = attrs?.Parent

    if (!ids && !parents) {
      items.push([featureLine])
    } else {
      let feature: GFF3Feature
      if (ids) {
        const id = ids[0]!
        const existing = byId.get(id)
        if (existing) {
          // Multi-location continuation: share child_features/derived_features
          // with the first line so children remain visible across all lines
          // regardless of arrival order.
          featureLine.child_features = existing[0]!.child_features
          featureLine.derived_features = existing[0]!.derived_features
          existing.push(featureLine)
          feature = existing
        } else {
          feature = [featureLine]
          if (!parents) {
            items.push(feature)
          }
          byId.set(id, feature)
          const waiting = orphans.get(id)
          if (waiting) {
            for (const w of waiting) {
              featureLine.child_features.push(w)
            }
            orphans.delete(id)
          }
        }
      } else {
        feature = [featureLine]
      }

      if (parents) {
        for (const parentId of parents) {
          const parent = byId.get(parentId)
          if (parent) {
            // child_features is shared across all parent feature lines,
            // so push once via the first line.
            parent[0]!.child_features.push(feature)
          } else {
            appendOrphan(orphans, parentId, feature)
          }
        }
      }
    }
  }

  return items
}

/**
 * Parse an array of LineRecord objects directly into JBrowse feature format.
 * Supports parent/child relationships via subfeatures.
 *
 * @param records - Array of LineRecord objects with raw line and metadata
 * @returns array of JBrowse-format features
 */
export function parseRecordsJBrowse(records: ParseInput[]): JBrowseFeature[] {
  const items: JBrowseFeature[] = []
  const byId = new Map<string, JBrowseFeature>()
  const orphans = new Map<string, JBrowseFeature[]>()

  for (const record of records) {
    const feature = record.hasEscapes
      ? parseFeatureJBrowse(record.line)
      : parseFeatureJBrowseNoUnescape(record.line)

    if (record.lineHash !== undefined) {
      feature._lineHash = String(record.lineHash)
    }

    const id = firstString(feature.id)
    const parents = toStringArray(feature.parent)

    if (!id && parents.length === 0) {
      items.push(feature)
    } else {
      // A parentless line is a top-level item. Every line of a top-level
      // discontinuous feature (e.g. cDNA_match/EST_match spanning several
      // segments under one shared ID, with no Parent) is its own top-level
      // item, so push here regardless of whether the id is already registered.
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
  }

  return items
}

export type {
  GFF3Comment,
  GFF3Directive,
  GFF3Feature,
  GFF3FeatureLine,
  GFF3FeatureLineWithRefs,
  GFF3Item,
  GFF3Sequence,
  JBrowseFeature,
} from './util.ts'
