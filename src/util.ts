// Fast, low-level functions for parsing GFF3.
// JavaScript port of Robert Buels's Bio::GFF3::LowLevel Perl module.

const PERCENT = '%'.charCodeAt(0)
const NOT_HEX = -1

// Value of each hex digit, indexed by char code, so that escapes are read
// without slicing out the digits and in any digit case (%eF is %EF).
const HEX_DIGIT = new Int8Array(128).fill(NOT_HEX)
for (let i = 0; i < 16; i++) {
  const digit = i.toString(16)
  HEX_DIGIT[digit.charCodeAt(0)] = i
  HEX_DIGIT[digit.toUpperCase().charCodeAt(0)] = i
}

const utf8Decoder = new TextDecoder()

function hexDigit(charCode: number) {
  return charCode < 128 ? HEX_DIGIT[charCode]! : NOT_HEX
}

/** Byte encoded by the escape at `i`, or NOT_HEX if `i` starts no valid escape. */
function escapeByteAt(s: string, i: number) {
  const hi =
    i + 2 < s.length && s.charCodeAt(i) === PERCENT
      ? hexDigit(s.charCodeAt(i + 1))
      : NOT_HEX
  const lo = hi === NOT_HEX ? NOT_HEX : hexDigit(s.charCodeAt(i + 2))
  return lo === NOT_HEX ? NOT_HEX : hi * 16 + lo
}

/** Index just past the run of consecutive escapes starting at `i`. */
function escapeRunEnd(s: string, i: number) {
  let end = i
  while (escapeByteAt(s, end) !== NOT_HEX) {
    end += 3
  }
  return end
}

/**
 * Decode the escapes in `s` between `start` and `end` as UTF-8, the encoding
 * GFF3 files are written in. Bytes that aren't valid UTF-8 become U+FFFD.
 */
function decodeEscapeRun(s: string, start: number, end: number) {
  const bytes = new Uint8Array((end - start) / 3)
  for (let i = start, b = 0; i < end; i += 3, b++) {
    bytes[b] = escapeByteAt(s, i)
  }
  return utf8Decoder.decode(bytes)
}

/**
 * Unescape a string value used in a GFF3 attribute.
 *
 * @param stringVal - Escaped GFF3 string value
 * @returns An unescaped string value
 */
export function unescape(stringVal: string) {
  let i = stringVal.indexOf('%')
  if (i === -1) {
    return stringVal
  }

  let result = ''
  let lastIdx = 0

  while (i !== -1) {
    const byte = escapeByteAt(stringVal, i)
    if (byte === NOT_HEX) {
      // Not a valid escape: keep looking from the next '%' so one that begins
      // a real escape isn't swallowed (e.g. the %20 in "a%b%20c").
      i = stringVal.indexOf('%', i + 1)
    } else if (byte < 0x80) {
      // An escaped byte below 0x80 stands for itself, no UTF-8 decoding needed
      result += stringVal.slice(lastIdx, i) + String.fromCharCode(byte)
      lastIdx = i + 3
      i = stringVal.indexOf('%', lastIdx)
    } else {
      // A byte >= 0x80 is the lead byte of a multi-byte UTF-8 character spread
      // over consecutive escapes (%E3%81%82 is あ), so decode the whole run
      // together rather than one byte at a time.
      const end = escapeRunEnd(stringVal, i)
      result += stringVal.slice(lastIdx, i) + decodeEscapeRun(stringVal, i, end)
      lastIdx = end
      i = stringVal.indexOf('%', end)
    }
  }

  return result + stringVal.slice(lastIdx)
}

// Columns past the end of a truncated line are treated the same as '.'
function isEmpty(s: string | undefined): s is undefined | '' | '.' {
  return s === undefined || s.length === 0 || s === '.'
}

function strField<E extends null | ''>(
  s: string | undefined,
  shouldUnescape: boolean,
  empty: E,
) {
  return isEmpty(s) ? empty : shouldUnescape ? unescape(s) : s
}

// Attribute names that collide with a field of the parsed feature are suffixed
// with '2'. 'subfeatures' is essential: an attribute of that name would replace
// the subfeature array with a string and break tree building. 'refname' does
// not collide with the camelCase 'refName' field, but is reserved so the two
// don't sit side by side.
const JBROWSE_DEFAULT_FIELDS = new Set([
  'start',
  'end',
  'seq_id',
  'refname',
  'score',
  'type',
  'source',
  'phase',
  'strand',
  'subfeatures',
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
  attrString: string | undefined,
  result: Record<string, unknown>,
  shouldUnescape: boolean,
) {
  if (
    attrString === undefined ||
    attrString.length === 0 ||
    attrString === '.'
  ) {
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

        // `Foo=,,` has no values at all, same as the `Foo=` that the eqIdx
        // check above already rejects, so don't leave an empty array behind
        if (values.length !== 0) {
          result[key] = values.length === 1 ? values[0] : values
        }
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
  const startStr = f[3]
  const endStr = f[4]
  const scoreStr = f[5]
  const strandStr = f[6]
  const phase = f[7]
  const attrString = f[8]

  const result: GffFeature = {
    refName: strField(f[0], shouldUnescape, ''),
    source: strField(f[1], shouldUnescape, null),
    type: strField(f[2], shouldUnescape, null),
    start: isEmpty(startStr) ? 0 : +startStr - 1,
    end: isEmpty(endStr) ? 0 : +endStr,
    score: isEmpty(scoreStr) ? undefined : +scoreStr,
    strand: strandStr === undefined ? undefined : STRAND_MAP[strandStr],
    phase: isEmpty(phase) ? undefined : +phase,
    subfeatures: [],
  }

  parseAttributes(attrString, result, shouldUnescape)
  return result
}
