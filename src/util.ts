// Fast, low-level functions for parsing GFF3.
// JavaScript port of Robert Buels's Bio::GFF3::LowLevel Perl module.

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
    const char =
      stringVal[i] === '%' && i + 2 < stringVal.length
        ? HEX_LOOKUP[stringVal.slice(i + 1, i + 3)]
        : undefined
    if (char !== undefined) {
      result += stringVal.slice(lastIdx, i) + char
      i += 3
      lastIdx = i
    } else {
      // Not a valid escape: advance one char so a '%' that begins a real
      // escape immediately after isn't swallowed (e.g. the %20 in "a%b%20c").
      i++
    }
  }

  return result + stringVal.slice(lastIdx)
}

function isEmpty(s: string) {
  return s.length === 0 || s === '.'
}

function strField<E extends null | ''>(
  s: string,
  shouldUnescape: boolean,
  empty: E,
) {
  return isEmpty(s) ? empty : shouldUnescape ? unescape(s) : s
}

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

const STRAND_MAP: Record<string, number | undefined> = {
  '+': 1,
  '-': -1,
  '.': 0,
}

/**
 * A parsed GFF3 feature: a flat object with 0-based half-open coordinates,
 * numeric strand (`1`/`-1`/`0`), attributes spread as lowercase top-level keys,
 * and child features nested under `subfeatures`.
 */
export interface GffFeature {
  start: number
  end: number
  strand?: number
  type: string | null
  source: string | null
  refName: string
  phase?: number
  score?: number
  subfeatures: GffFeature[]
  [key: string]: unknown
}

/**
 * Parse the 9th column (attributes) of a GFF3 feature line into `result`,
 * lowercasing keys and suffixing any that collide with a default field name.
 * Pass shouldUnescape=false as a fast path for data with no escaped characters.
 */
export function parseAttributes(
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

/**
 * Parse a GFF3 feature line. Unescaping is skipped entirely for lines with no
 * '%' character, which is the common case.
 *
 * @param line - GFF3 feature line
 * @returns The parsed feature
 */
export function parseFeature(line: string): GffFeature {
  const f = line.split('\t')
  const shouldUnescape = line.includes('%')
  const startStr = f[3]!
  const endStr = f[4]!
  const scoreStr = f[5]!
  const phase = f[7]!
  const attrString = f[8]!

  const result: GffFeature = {
    refName: strField(f[0]!, shouldUnescape, ''),
    source: strField(f[1]!, shouldUnescape, null),
    type: strField(f[2]!, shouldUnescape, null),
    start: isEmpty(startStr) ? 0 : +startStr - 1,
    end: isEmpty(endStr) ? 0 : +endStr,
    score: isEmpty(scoreStr) ? undefined : +scoreStr,
    strand: STRAND_MAP[f[6]!],
    phase: isEmpty(phase) ? undefined : +phase,
    subfeatures: [],
  }

  parseAttributes(attrString, result, shouldUnescape)
  return result
}
