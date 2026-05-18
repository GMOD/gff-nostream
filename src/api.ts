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

export interface LineRecord {
  line: string
  lineHash?: string | number
  /** Optional passthrough byte offsets — not used by the parser */
  start?: number
  /** Optional passthrough byte offsets — not used by the parser */
  end?: number
  hasEscapes: boolean
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
  const records: LineRecord[] = []
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
export function parseRecords(records: LineRecord[]): GFF3Feature[] {
  const items: GFF3Feature[] = []
  const byId = new Map<string, GFF3Feature>()
  const orphans = new Map<string, GFF3Feature[]>()

  for (const record of records) {
    const featureLine = (
      record.hasEscapes
        ? parseFeature(record.line)
        : parseFeatureNoUnescape(record.line)
    ) as GFF3FeatureLineWithRefs
    featureLine.child_features = []
    featureLine.derived_features = []

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
export function parseRecordsJBrowse(records: LineRecord[]): JBrowseFeature[] {
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

    // attribute parsing collapses single-element arrays to scalars, so id can
    // be string | string[]; defensively take the first if multi-valued.
    const rawId = feature.id as string | string[] | undefined
    const id = Array.isArray(rawId) ? rawId[0] : rawId
    const parent = feature.parent as string | string[] | undefined

    if (!id && !parent) {
      items.push(feature)
    } else if (!id || !byId.has(id)) {
      if (id) {
        if (!parent) {
          items.push(feature)
        }
        byId.set(id, feature)
        const waiting = orphans.get(id)
        if (waiting) {
          for (const w of waiting) {
            feature.subfeatures.push(w)
          }
          orphans.delete(id)
        }
      }

      if (parent) {
        const parents = Array.isArray(parent) ? parent : [parent]
        for (const parentId of parents) {
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
