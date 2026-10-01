import {
  attributeStringHasId,
  getLinkAttributes,
  parseFeature,
  parseFeatureLazy,
} from './util.ts'

import type { GffFeature, LazyGffFeature } from './util.ts'

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

/** {@link ParsedRecord} whose feature still carries column 9 as raw text. */
export interface ParsedLazyRecord<R extends LineRecord = LineRecord> {
  feature: LazyGffFeature
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

/**
 * Whether a raw GFF3 line carries an `ID` attribute, which is the exact test for
 * whether it can have children: a child names its parent with `Parent=<ID>`, so
 * a record with no `ID` can be referenced by nothing. Decided by the same scan
 * the parser links with, so a line this admits is one whose children the tree
 * will actually attach, and one it rejects is one they never would be.
 *
 * Column 9 is found by counting tabs rather than from the last tab, and is
 * bounded at a stray tab inside it, exactly as {@link parseFeatureLazy} bounds
 * it. A line with fewer than nine columns has no attributes to carry an ID.
 */
export function hasIdAttribute(line: string) {
  let p = 0
  for (let i = 0; i < 8; i++) {
    const t = line.indexOf('\t', p)
    if (t === -1) {
      return false
    }
    p = t + 1
  }
  const attrEnd = line.indexOf('\t', p)
  return attributeStringHasId(
    line.slice(p, attrEnd === -1 ? undefined : attrEnd),
  )
}

/** Append a value to the array stored under key, creating the array if absent. */
function appendTo<K, T>(map: Map<K, T[]>, key: K, value: T) {
  const arr = map.get(key)
  if (arr) {
    arr.push(value)
  } else {
    map.set(key, [value])
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
 * - `folded`: a further line of a parentless feature already collected, now a
 *   segment under it
 */
type LinkStatus = 'top-level' | 'attached' | 'orphaned' | 'folded'

interface Linkable<F> {
  type: string | null
  start: number
  end: number
  subfeatures: F[]
}

/**
 * One parse's tree under construction. `split` holds the parentless features
 * already folded into containers; `continuations` holds the further same-type
 * lines of each parented ID, keyed by its first line.
 */
interface LinkState<F> {
  byId: Map<string, F>
  orphans: Map<string, F[]>
  split: Set<F>
  continuations: Map<F, F[]>
}

function newLinkState<F>(): LinkState<F> {
  return {
    byId: new Map(),
    orphans: new Map(),
    split: new Set(),
    continuations: new Map(),
  }
}

function parentsOf(feature: GffFeature) {
  return toStringArray(feature.parent)
}

function parentsOfLazy(feature: LazyGffFeature) {
  return toStringArray(getLinkAttributes(feature).parent)
}

/**
 * Register a feature's ID and attach it to its parent(s), building the
 * subfeature tree in `state`.
 */
function linkFeature(
  feature: GffFeature,
  state: LinkState<GffFeature>,
): LinkStatus {
  return linkResolved(
    feature,
    firstString(feature.id),
    parentsOf(feature),
    state,
  )
}

/**
 * {@link linkFeature} for a feature whose attributes are still raw text: the
 * tree-building logic is identical, only the two attributes it needs are read
 * with a targeted scan rather than off already-parsed properties.
 */
function linkFeatureLazy(
  feature: LazyGffFeature,
  state: LinkState<LazyGffFeature>,
): LinkStatus {
  const { id, parent } = getLinkAttributes(feature)
  return linkResolved(feature, firstString(id), toStringArray(parent), state)
}

/**
 * A further line of a parentless discontinuous feature (a cDNA_match written
 * one line per aligned block, all under one ID) folds into the first: that line
 * becomes the feature, spanning every segment, and each segment — its own
 * included — hangs under it as a subfeature of the same type.
 */
function foldSegment<F extends Linkable<F>>(
  container: F,
  segment: F,
  split: Set<F>,
) {
  if (!split.has(container)) {
    split.add(container)
    container.subfeatures = [{ ...container, subfeatures: [] }]
  }
  container.start = Math.min(container.start, segment.start)
  container.end = Math.max(container.end, segment.end)
  container.subfeatures.push(segment)
}

/**
 * The link-and-attach step itself, over an id and parent list the caller has
 * already extracted — the one copy shared by the eager and lazy parsers, which
 * differ only in where those two values come from.
 */
function linkResolved<F extends Linkable<F>>(
  feature: F,
  id: string | undefined,
  parents: string[],
  { byId, orphans, split, continuations }: LinkState<F>,
): LinkStatus {
  // Only the first line of an ID registers it. A further line with no parent
  // has nowhere else to go, so it folds into that first line — unless the first
  // line already has children, or is of another type: then the ID is a
  // duplicate (two genes both named by their symbol), not one feature over
  // several lines, and the line stands on its own. A further line with parents
  // attaches to them as its own segment, and foldContinuations decides at the
  // end whether its ID turned out to be one feature.
  if (id) {
    const first = byId.get(id)
    if (first === undefined) {
      byId.set(id, feature)
      const waiting = orphans.get(id)
      if (waiting) {
        for (const w of waiting) {
          feature.subfeatures.push(w)
        }
        orphans.delete(id)
      }
    } else if (first.type === feature.type) {
      if (parents.length > 0) {
        appendTo(continuations, first, feature)
      } else if (split.has(first) || first.subfeatures.length === 0) {
        foldSegment(first, feature, split)
        return 'folded'
      }
    }
  }

  let attached = false
  for (const parentId of parents) {
    const parentFeature = byId.get(parentId)
    if (parentFeature) {
      parentFeature.subfeatures.push(feature)
      attached = true
    } else {
      appendTo(orphans, parentId, feature)
    }
  }

  return parents.length === 0 ? 'top-level' : attached ? 'attached' : 'orphaned'
}

function sameParents(a: string[], b: string[]) {
  return a.length === b.length && a.every((p, i) => p === b[i])
}

/**
 * A parented ID written over several lines is one discontinuous feature when
 * its first line has children: NCBI writes a ribosomal-frameshift polyprotein's
 * CDS (SARS-CoV-2 ORF1ab, HIV-1 gag-pol) as one line per reading frame and
 * hangs its mature peptides off that ID. The first line becomes the feature,
 * spanning every line, with its segments (a childless copy of itself first)
 * followed by its children; the further lines leave their parents. A childless
 * repeated ID — every CDS of an NCBI or Ensembl transcript — keeps one sibling
 * per line, as does a further line naming other parents than the first.
 *
 * @returns the lines folded away, which an entry point must not collect as
 *   unparented
 */
function foldContinuations<F extends Linkable<F>>(
  { byId, continuations }: LinkState<F>,
  getParents: (feature: F) => string[],
) {
  const folded = new Set<F>()
  for (const [first, lines] of continuations) {
    if (first.subfeatures.length > 0) {
      const parents = getParents(first)
      const segments = lines.filter(line =>
        sameParents(getParents(line), parents),
      )
      if (segments.length > 0) {
        first.subfeatures = [
          { ...first, subfeatures: [] },
          ...segments,
          ...first.subfeatures,
        ]
        for (const segment of segments) {
          first.start = Math.min(first.start, segment.start)
          first.end = Math.max(first.end, segment.end)
          folded.add(segment)
        }
        for (const parentId of parents) {
          const parent = byId.get(parentId)
          if (parent) {
            parent.subfeatures = parent.subfeatures.filter(f => !folded.has(f))
          }
        }
      }
    }
  }
  return folded
}

/**
 * True when none of a feature's Parent ids were ever defined in the input, so
 * it was never nested anywhere. Registering an id adopts everything waiting on
 * it, so presence in `byId` after the full pass means the feature was attached.
 */
function isUnparented(feature: GffFeature, byId: Map<string, GffFeature>) {
  return !parentsOf(feature).some(parentId => byId.has(parentId))
}

/**
 * {@link isUnparented} for a lazily-parsed feature. Rescanning the attribute
 * string is fine here: this runs only over the features still orphaned after
 * the whole input has been read, not over every line.
 */
function isUnparentedLazy(
  feature: LazyGffFeature,
  byId: Map<string, LazyGffFeature>,
) {
  return !parentsOfLazy(feature).some(parentId => byId.has(parentId))
}

/*
 * The five entry points below run the same link-and-collect loop over
 * differently-shaped input. Routing them through one core that reads lines via
 * a callback was tried and reverted — it put every caller on a shared
 * polymorphic call site for no measured gain. The duplication is deliberate;
 * keep the loops in sync by hand. (The link step and the closing fold are *not*
 * duplicated: those are `linkResolved` and `foldContinuations`, which the eager
 * and lazy paths share.)
 */

/**
 * Parse an array of raw GFF3 feature lines, resolving parent/child
 * relationships into `subfeatures`. The lines must already be free of blanks,
 * comments, and any `##FASTA` section — this is the entry point for a caller
 * that has split and filtered the file itself (a tabix region query, or a
 * whole-file scan grouping lines by reference sequence) and has no per-line
 * identity to carry through.
 * Features whose Parent is never defined in `lines` are returned at the end as
 * top-level items rather than dropped.
 *
 * @param lines - raw GFF3 feature lines
 * @returns top-level features
 */
export function parseLines(lines: readonly string[]): GffFeature[] {
  const items: GffFeature[] = []
  const pending: GffFeature[] = []
  const state = newLinkState<GffFeature>()

  for (const line of lines) {
    const feature = parseFeature(line)
    const status = linkFeature(feature, state)
    if (status === 'top-level') {
      items.push(feature)
    } else if (status === 'orphaned') {
      pending.push(feature)
    }
  }

  const folded = foldContinuations(state, parentsOf)
  for (const feature of pending) {
    if (!folded.has(feature) && isUnparented(feature, state.byId)) {
      items.push(feature)
    }
  }

  return items
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
  const state = newLinkState<GffFeature>()
  const pending: GffFeature[] = []

  // filters and parses in one pass rather than collecting the kept lines and
  // handing them to parseLines, which would hold a second array the size of the
  // file alongside the split
  for (const line of str.split(/\r?\n/)) {
    if (line.startsWith('##FASTA') || line.startsWith('>')) {
      break
    }
    if (line.length !== 0 && !line.startsWith('#')) {
      const feature = parseFeature(line)
      const status = linkFeature(feature, state)
      if (status === 'top-level') {
        items.push(feature)
      } else if (status === 'orphaned') {
        pending.push(feature)
      }
    }
  }

  const folded = foldContinuations(state, parentsOf)
  for (const feature of pending) {
    if (!folded.has(feature) && isUnparented(feature, state.byId)) {
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
  const state = newLinkState<GffFeature>()
  const pending: ParsedRecord<R>[] = []

  for (const record of records) {
    const feature = parseFeature(record.line)
    const status = linkFeature(feature, state)
    if (status === 'top-level') {
      items.push({ feature, record })
    } else if (status === 'orphaned') {
      pending.push({ feature, record })
    }
  }

  const folded = foldContinuations(state, parentsOf)
  for (const parsed of pending) {
    if (
      !folded.has(parsed.feature) &&
      isUnparented(parsed.feature, state.byId)
    ) {
      items.push(parsed)
    }
  }

  return items
}

/**
 * {@link parseLines}, leaving each feature's attributes as raw column-9 text.
 * Only ID and Parent are read, because the parent/child tree cannot be built
 * without them; everything else is materialized on demand — see
 * {@link LazyGffFeature} for when that is the right trade.
 *
 * @param lines - raw GFF3 feature lines
 * @returns top-level features, attributes unparsed
 */
export function parseLinesLazy(lines: readonly string[]): LazyGffFeature[] {
  const items: LazyGffFeature[] = []
  const pending: LazyGffFeature[] = []
  const state = newLinkState<LazyGffFeature>()

  for (const line of lines) {
    const feature = parseFeatureLazy(line)
    const status = linkFeatureLazy(feature, state)
    if (status === 'top-level') {
      items.push(feature)
    } else if (status === 'orphaned') {
      pending.push(feature)
    }
  }

  const folded = foldContinuations(state, parentsOfLazy)
  for (const feature of pending) {
    if (!folded.has(feature) && isUnparentedLazy(feature, state.byId)) {
      items.push(feature)
    }
  }

  return items
}

/**
 * {@link parseRecords}, leaving each feature's attributes as raw column-9 text.
 * The lazy counterpart for a caller that also carries its own per-line identity
 * (a tabix byte offset, say) — see {@link LazyGffFeature}.
 *
 * @param records - Array of records, each carrying a raw GFF3 `line`
 * @returns top-level features with attributes unparsed, each paired with its
 *   originating record
 */
export function parseRecordsLazy<R extends LineRecord>(
  records: readonly R[],
): ParsedLazyRecord<R>[] {
  const items: ParsedLazyRecord<R>[] = []
  const state = newLinkState<LazyGffFeature>()
  const pending: ParsedLazyRecord<R>[] = []

  for (const record of records) {
    const feature = parseFeatureLazy(record.line)
    const status = linkFeatureLazy(feature, state)
    if (status === 'top-level') {
      items.push({ feature, record })
    } else if (status === 'orphaned') {
      pending.push({ feature, record })
    }
  }

  const folded = foldContinuations(state, parentsOfLazy)
  for (const parsed of pending) {
    if (
      !folded.has(parsed.feature) &&
      isUnparentedLazy(parsed.feature, state.byId)
    ) {
      items.push(parsed)
    }
  }

  return items
}

export {
  getAttribute,
  getAttributes,
  getLinkAttributes,
  parseFeatureLazy,
} from './util.ts'

export type { GffFeature, LazyGffFeature } from './util.ts'
