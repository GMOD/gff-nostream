// Fast, low-level functions for parsing and formatting GFF3.
// JavaScript port of Robert Buels's Bio::GFF3::LowLevel Perl module.

const directiveRegex = /^\s*##\s*(\S+)\s*(.*)/
const whitespaceRegex = /\s+/
const nonDigitRegex = /\D/g

const HEX_LOOKUP: Record<string, string | undefined> = {}
for (let i = 0; i < 256; i++) {
  const hex = i.toString(16).toUpperCase().padStart(2, '0')
  HEX_LOOKUP[hex] = String.fromCharCode(i)
  HEX_LOOKUP[hex.toLowerCase()] = String.fromCharCode(i)
}

/**
 * Unescape a string value used in a GFF3 attribute.
 *
 * @param stringVal - Escaped GFF3 string value
 * @returns An unescaped string value
 */

export function unescape(stringVal: string) {
  const idx = stringVal.indexOf('%')
  if (idx === -1) {
    return stringVal
  }

  let result = ''
  let lastIdx = 0
  let i = idx

  while (i < stringVal.length) {
    if (stringVal[i] === '%' && i + 2 < stringVal.length) {
      result += stringVal.slice(lastIdx, i)
      const hex = stringVal.slice(i + 1, i + 3)
      const char = HEX_LOOKUP[hex]
      if (char !== undefined) {
        result += char
      } else {
        result += stringVal.slice(i, i + 3)
      }
      i += 3
      lastIdx = i
    } else {
      i++
    }
  }

  return result + stringVal.slice(lastIdx)
}

function parseAttributesImpl(
  attrString: string,
  shouldUnescape: boolean,
): GFF3Attributes {
  if (attrString.length === 0 || attrString === '.') {
    return {}
  }

  const attrs: GFF3Attributes = {}
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
        let arec = attrs[tag]
        if (!arec) {
          arec = []
          attrs[tag] = arec
        }

        let valStart = eqIdx + 1
        while (valStart < semiIdx) {
          let commaIdx = attrString.indexOf(',', valStart)
          if (commaIdx === -1 || commaIdx > semiIdx) {
            commaIdx = semiIdx
          }
          if (commaIdx > valStart) {
            const val = attrString.slice(valStart, commaIdx)
            arec.push(shouldUnescape ? unescape(val) : val)
          }
          valStart = commaIdx + 1
        }
      }
    }
    start = semiIdx + 1
  }
  return attrs
}

/**
 * Parse the 9th column (attributes) of a GFF3 feature line.
 *
 * @param attrString - String of GFF3 9th column
 * @returns Parsed attributes
 */
export function parseAttributes(attrString: string): GFF3Attributes {
  return parseAttributesImpl(attrString, true)
}

/**
 * Parse the 9th column (attributes) of a GFF3 feature line without unescaping.
 * Fast path for data known to contain no escaped characters.
 *
 * @param attrString - String of GFF3 9th column
 * @returns Parsed attributes
 */
export function parseAttributesNoUnescape(attrString: string): GFF3Attributes {
  return parseAttributesImpl(attrString, false)
}

function normImpl(s: string, shouldUnescape: boolean) {
  if (s.length === 0 || s === '.') {
    return null
  }
  return shouldUnescape ? unescape(s) : s
}

function parseFeatureImpl(
  line: string,
  shouldUnescape: boolean,
): GFF3FeatureLine {
  const f = line.split('\t')
  const startStr = f[3]!
  const endStr = f[4]!
  const scoreStr = f[5]!
  const strand = f[6]!
  const phase = f[7]!
  const attrString = f[8]!

  return {
    seq_id: normImpl(f[0]!, shouldUnescape),
    source: normImpl(f[1]!, shouldUnescape),
    type: normImpl(f[2]!, shouldUnescape),
    start: startStr.length === 0 || startStr === '.' ? null : +startStr,
    end: endStr.length === 0 || endStr === '.' ? null : +endStr,
    score: scoreStr.length === 0 || scoreStr === '.' ? null : +scoreStr,
    strand: normImpl(strand, false),
    phase: normImpl(phase, false),
    attributes:
      attrString.length === 0 || attrString === '.'
        ? null
        : parseAttributesImpl(attrString, shouldUnescape),
  }
}

/**
 * Parse a GFF3 feature line
 *
 * @param line - GFF3 feature line
 * @returns The parsed feature
 */
export function parseFeature(line: string): GFF3FeatureLine {
  return parseFeatureImpl(line, true)
}

/**
 * Parse a GFF3 feature line without unescaping.
 * Fast path for data known to contain no escaped characters.
 *
 * @param line - GFF3 feature line
 * @returns The parsed feature
 */
export function parseFeatureNoUnescape(line: string): GFF3FeatureLine {
  return parseFeatureImpl(line, false)
}

function parseFieldsArrayImpl(
  f: (string | null | undefined)[],
  shouldUnescape: boolean,
): GFF3FeatureLine {
  const seq_id = f[0]
  const source = f[1]
  const type = f[2]
  const startStr = f[3]
  const endStr = f[4]
  const scoreStr = f[5]
  const strand = f[6]
  const phase = f[7]
  const attrString = f[8]

  return {
    seq_id: seq_id ? normImpl(seq_id, shouldUnescape) : null,
    source: source ? normImpl(source, shouldUnescape) : null,
    type: type ? normImpl(type, shouldUnescape) : null,
    start: !startStr || startStr === '.' ? null : +startStr,
    end: !endStr || endStr === '.' ? null : +endStr,
    score: !scoreStr || scoreStr === '.' ? null : +scoreStr,
    strand: strand && strand !== '.' ? strand : null,
    phase: phase && phase !== '.' ? phase : null,
    attributes:
      !attrString || attrString === '.'
        ? null
        : parseAttributesImpl(attrString, shouldUnescape),
  }
}

/**
 * Parse a GFF3 feature from a pre-split fields array
 *
 * @param f - Array of 9 GFF3 column values (use null or '.' for empty values)
 * @returns The parsed feature
 */
export function parseFieldsArray(
  f: (string | null | undefined)[],
): GFF3FeatureLine {
  return parseFieldsArrayImpl(f, true)
}

/**
 * Parse a GFF3 feature from a pre-split fields array without unescaping.
 * Fast path for data known to contain no escaped characters.
 *
 * @param f - Array of 9 GFF3 column values (use null or '.' for empty values)
 * @returns The parsed feature
 */
export function parseFieldsArrayNoUnescape(
  f: (string | null | undefined)[],
): GFF3FeatureLine {
  return parseFieldsArrayImpl(f, false)
}

/**
 * Parse a GFF3 directive line.
 *
 * @param line - GFF3 directive line
 * @returns The parsed directive
 */
export function parseDirective(
  line: string,
):
  | GFF3Directive
  | GFF3SequenceRegionDirective
  | GFF3GenomeBuildDirective
  | null {
  const match = directiveRegex.exec(line)
  if (!match) {
    return null
  }

  const name = match[1]!
  let contents = match[2]!

  const parsed: GFF3Directive = { directive: name }
  if (contents.length) {
    parsed.value = contents.trimEnd()
  }

  if (name === 'sequence-region') {
    const c = contents.split(whitespaceRegex, 3)
    return {
      ...parsed,
      seq_id: c[0]!,
      start: c[1]!.replaceAll(nonDigitRegex, ''),
      end: c[2]!.replaceAll(nonDigitRegex, ''),
    }
  } else if (name === 'genome-build') {
    const [source, buildName] = contents.split(whitespaceRegex, 2)
    return {
      ...parsed,
      source: source!,
      buildName: buildName!,
    }
  }

  return parsed
}

/** A record of GFF3 attribute identifiers and the values of those identifiers */
export type GFF3Attributes = Record<string, string[] | undefined>

/** A representation of a single line of a GFF3 file */
export interface GFF3FeatureLine {
  /** The ID of the landmark used to establish the coordinate system for the current feature */
  seq_id: string | null
  /** A free text qualifier intended to describe the algorithm or operating procedure that generated this feature */
  source: string | null
  /** The type of the feature */
  type: string | null
  /** The start coordinates of the feature */
  start: number | null
  /** The end coordinates of the feature */
  end: number | null
  /** The score of the feature */
  score: number | null
  /** The strand of the feature */
  strand: string | null
  /** For features of type "CDS", the phase indicates where the next codon begins relative to the 5' end of the current CDS feature */
  phase: string | null
  /** Feature attributes */
  attributes: GFF3Attributes | null
}

/**
 * A GFF3 Feature line that includes references to other features defined in
 * their "Parent" or "Derives_from" attributes
 */
export interface GFF3FeatureLineWithRefs extends GFF3FeatureLine {
  /** An array of child features */
  child_features: GFF3Feature[]
  /** An array of features derived from this feature */
  derived_features: GFF3Feature[]
}

/**
 * A GFF3 feature, which may include multiple individual feature lines
 */
export type GFF3Feature = GFF3FeatureLineWithRefs[]

/** A GFF3 directive */
export interface GFF3Directive {
  /** The name of the directive */
  directive: string
  /** The string value of the directive */
  value?: string
}

/** A GFF3 sequence-region directive */
export interface GFF3SequenceRegionDirective extends GFF3Directive {
  /** The string value of the directive */
  value: string
  /** The sequence ID parsed from the directive */
  seq_id: string
  /** The sequence start parsed from the directive */
  start: string
  /** The sequence end parsed from the directive */
  end: string
}

/** A GFF3 genome-build directive */
export interface GFF3GenomeBuildDirective extends GFF3Directive {
  /** The string value of the directive */
  value: string
  /** The genome build source parsed from the directive */
  source: string
  /** The genome build name parsed from the directive */
  buildName: string
}

/** A GFF3 comment */
export interface GFF3Comment {
  /** The text of the comment */
  comment: string
}

/** A GFF3 FASTA single sequence */
export interface GFF3Sequence {
  /** The ID of the sequence */
  id: string
  /** The description of the sequence */
  description?: string
  /** The sequence */
  sequence: string
}

export type GFF3Item = GFF3Feature | GFF3Directive | GFF3Comment | GFF3Sequence

// JBrowse format types and parsing functions

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

// Pre-computed lowercase for common GFF3 spec attribute names to avoid
// toLowerCase() calls in the hot path
const COMMON_ATTRS: Record<string, string | undefined> = {
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
  id: 'id',
  name: 'name',
  parent: 'parent',
  note: 'note',
  dbxref: 'dbxref',
  alias: 'alias',
  target: 'target',
  gap: 'gap',
}

export interface JBrowseFeature {
  start: number
  end: number
  strand?: number
  type: string | null
  source: string | null
  refName: string
  phase?: number
  score?: number
  subfeatures: JBrowseFeature[]
  [key: string]: unknown
}

function parseAttributesJBrowseImpl(
  attrString: string,
  result: Record<string, unknown>,
  shouldUnescape: boolean,
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
        let key = COMMON_ATTRS[tag]
        if (key === undefined) {
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
            values.push(shouldUnescape ? unescape(val) : val)
          }
          valStart = commaIdx + 1
        }

        result[key] = values.length === 1 ? values[0] : values
      }
    }
    start = semiIdx + 1
  }
}

export function parseAttributesJBrowse(
  attrString: string,
  result: Record<string, unknown>,
) {
  parseAttributesJBrowseImpl(attrString, result, true)
}

export function parseAttributesJBrowseNoUnescape(
  attrString: string,
  result: Record<string, unknown>,
) {
  parseAttributesJBrowseImpl(attrString, result, false)
}

function parseFeatureJBrowseImpl(
  line: string,
  shouldUnescape: boolean,
): JBrowseFeature {
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

  const result: JBrowseFeature = {
    refName:
      seq_id.length === 0 || seq_id === '.'
        ? ''
        : shouldUnescape
          ? unescape(seq_id)
          : seq_id,
    source:
      source.length === 0 || source === '.'
        ? null
        : shouldUnescape
          ? unescape(source)
          : source,
    type:
      type.length === 0 || type === '.'
        ? null
        : shouldUnescape
          ? unescape(type)
          : type,
    start: startStr.length === 0 || startStr === '.' ? 0 : +startStr - 1,
    end: endStr.length === 0 || endStr === '.' ? 0 : +endStr,
    score: scoreStr.length === 0 || scoreStr === '.' ? undefined : +scoreStr,
    strand:
      strand === '+' ? 1 : strand === '-' ? -1 : strand === '.' ? 0 : undefined,
    phase: phase.length === 0 || phase === '.' ? undefined : +phase,
    subfeatures: [],
  }

  parseAttributesJBrowseImpl(attrString, result, shouldUnescape)
  return result
}

export function parseFeatureJBrowse(line: string): JBrowseFeature {
  return parseFeatureJBrowseImpl(line, true)
}

export function parseFeatureJBrowseNoUnescape(line: string): JBrowseFeature {
  return parseFeatureJBrowseImpl(line, false)
}
