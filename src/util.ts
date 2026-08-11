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

/** Index of `char` at or after `from`, or `end` if it does not occur before it. */
function indexOrEnd(s: string, char: string, from: number, end: number) {
  const idx = s.indexOf(char, from)
  return idx === -1 || idx > end ? end : idx
}

/** Raw lines handed to parseRecords can still carry their line terminator. */
function trimLineEnd(s: string) {
  let end = s.length
  if (s[end - 1] === '\n') {
    end -= s[end - 2] === '\r' ? 2 : 1
  }
  return end === s.length ? s : s.slice(0, end)
}

/** Feature key an attribute tag is stored under: lowercased, and suffixed if reserved. */
function attributeKey(tag: string) {
  const common = COMMON_ATTRS[tag]
  if (common === undefined) {
    const key = tag.toLowerCase()
    return JBROWSE_DEFAULT_FIELDS.has(key) ? `${key}2` : key
  } else {
    return common
  }
}

/** Comma-separated attribute values between `from` and `end`. */
function parseValues(
  s: string,
  from: number,
  end: number,
  shouldUnescape: boolean,
) {
  const values: string[] = []
  let valStart = from
  while (valStart < end) {
    const commaIdx = indexOrEnd(s, ',', valStart, end)
    if (commaIdx > valStart) {
      const val = s.slice(valStart, commaIdx)
      values.push(shouldUnescape ? unescape(val) : val)
    }
    valStart = commaIdx + 1
  }
  return values
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

  const attrs = trimLineEnd(attrString)
  const len = attrs.length

  let start = 0
  while (start < len) {
    const semiIdx = indexOrEnd(attrs, ';', start, len)
    const eqIdx = attrs.indexOf('=', start)

    // Nothing to store for a tag whose '=' belongs to a later attribute, or
    // that has no value before the next ';' — including `Foo=` and `Foo=,,`
    if (eqIdx !== -1 && eqIdx + 1 < semiIdx) {
      const values = parseValues(attrs, eqIdx + 1, semiIdx, shouldUnescape)
      if (values.length !== 0) {
        const key = attributeKey(attrs.slice(start, eqIdx))
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
/**
 * A parsed GFF3 feature whose attributes are still the raw text of column 9.
 *
 * The eight fixed columns are parsed exactly as {@link GffFeature}'s are; what
 * is deferred is turning column 9 into object keys, which on annotation-grade
 * input (GENCODE, NCBI, Ensembl — 15-20 attributes a line) is about two thirds
 * of the cost of parsing a line. Read attributes with {@link getAttribute} for
 * one key, or {@link getAttributes} to materialize them all.
 *
 * Worth deferring only when most attributes are never read. A caller that will
 * read every attribute of every feature should use {@link parseFeature}, which
 * parses each attribute string once instead of once per lookup.
 */
export interface LazyGffFeature {
  start: number
  end: number
  strand?: number
  type: string | null
  source: string | null
  refName: string
  phase?: number
  score?: number
  subfeatures: LazyGffFeature[]
  /**
   * Raw, unparsed column 9. Present so attributes can be materialized later;
   * it is a normal enumerable property, so a caller that spreads or
   * JSON-stringifies the feature gets this string and *not* the attributes.
   */
  attributeString: string
}

/**
 * The raw tags (lowercased) that {@link attributeKey} maps onto `key`, as a
 * primary and an optional alternate — two `string | undefined` locals rather
 * than an array, because this is resolved on every lookup and an array here
 * costs more than the scan it feeds.
 *
 * Every COMMON_ATTRS entry equals its tag's `toLowerCase()`, so that table is
 * purely a fast path and the mapping is exactly "lowercase, then suffix '2' if
 * the result is reserved". Inverting it: a reserved key is unreachable (any tag
 * lowercasing to it is suffixed away), and a key ending in '2' whose stem is
 * reserved is reachable two ways — from the stem (`Start` -> `start2`) and from
 * a tag literally named `start2`. Which of those wins is decided by position,
 * not by preference: the scan keeps the last match either way, as the eager
 * parser's last-wins assignment does.
 */
function primaryTag(key: string) {
  return JBROWSE_DEFAULT_FIELDS.has(key) ? undefined : key
}

function alternateTag(key: string) {
  const stem = key.endsWith('2') ? key.slice(0, -1) : undefined
  return stem !== undefined && JBROWSE_DEFAULT_FIELDS.has(stem)
    ? stem
    : undefined
}

/**
 * Case-insensitive compare of `s[start,end)` against an already-lowercased
 * `lower`, without slicing the tag out first. ASCII is folded inline; a tag
 * containing anything else falls back to `toLowerCase`, which is the only way
 * to be right about non-ASCII case folding (and about the rare fold that
 * changes length, hence the second loop).
 */
function tagMatches(s: string, start: number, end: number, lower: string) {
  if (end - start === lower.length) {
    for (let i = 0; i < lower.length; i++) {
      const c = s.charCodeAt(start + i)
      if (c >= 128) {
        return s.slice(start, end).toLowerCase() === lower
      }
      const lc = c >= 65 && c <= 90 ? c + 32 : c
      if (lc !== lower.charCodeAt(i)) {
        return false
      }
    }
    return true
  }
  for (let i = start; i < end; i++) {
    if (s.charCodeAt(i) >= 128) {
      return s.slice(start, end).toLowerCase() === lower
    }
  }
  return false
}

/**
 * Values for several attribute keys in one pass over the attribute string,
 * returned positionally. Mirrors {@link parseAttributes}'s scan exactly —
 * same delimiter handling, same "skip a tag with no value" rule, same
 * single-value collapse — but slices and unescapes a value only for a tag that
 * matches one of the wanted keys. Later occurrences overwrite earlier ones,
 * matching the eager parser's last-wins assignment.
 *
 * One pass for N keys because the linking loop wants ID and Parent together,
 * and scanning twice for them would undo much of the point.
 */
/** True when the tag at `[start,end)` is `tag`, or `alt` if one is given. */
function tagIs(
  s: string,
  start: number,
  end: number,
  tag: string | undefined,
  alt: string | undefined,
) {
  return (
    (tag !== undefined && tagMatches(s, start, end, tag)) ||
    (alt !== undefined && tagMatches(s, start, end, alt))
  )
}

/**
 * Values for up to two attribute slots in one pass. Mirrors
 * {@link parseAttributes}'s scan exactly — same delimiter handling, same "skip
 * a tag with no value" rule, same single-value collapse — but slices and
 * unescapes a value only for a tag that matches a wanted slot, and keeps the
 * last match, as the eager parser's last-wins assignment does.
 *
 * Two slots rather than a list because the only multi-key caller is the linking
 * loop, which wants ID and Parent together and would otherwise scan twice. It
 * is written to allocate nothing per attribute: an earlier version took a
 * `string[]` of keys and tested them with `.some()`, which cost an array per
 * call and a closure per attribute per key, and measured *slower* than parsing
 * every attribute eagerly.
 *
 * Values are always unescaped rather than gated on the attribute string
 * containing a '%'. `unescape` already returns its input untouched when the
 * value has none, so the gate would only add a scan of the whole string to
 * every lookup.
 */
function scanAttributes(
  attrString: string,
  tag0: string | undefined,
  alt0: string | undefined,
  tag1: string | undefined,
  alt1: string | undefined,
) {
  let value0: unknown
  let value1: unknown
  if (
    attrString.length === 0 ||
    attrString === '.' ||
    (tag0 === undefined && alt0 === undefined && tag1 === undefined)
  ) {
    return { value0, value1 }
  }

  const attrs = trimLineEnd(attrString)
  const len = attrs.length

  let start = 0
  while (start < len) {
    const semiIdx = indexOrEnd(attrs, ';', start, len)
    const eqIdx = attrs.indexOf('=', start)
    if (eqIdx !== -1 && eqIdx + 1 < semiIdx) {
      const slot0 = tagIs(attrs, start, eqIdx, tag0, alt0)
      if (slot0 || tagIs(attrs, start, eqIdx, tag1, alt1)) {
        const values = parseValues(attrs, eqIdx + 1, semiIdx, true)
        if (values.length !== 0) {
          const value = values.length === 1 ? values[0] : values
          if (slot0) {
            value0 = value
          } else {
            value1 = value
          }
        }
      }
    }
    start = semiIdx + 1
  }
  return { value0, value1 }
}

/** The value of one attribute, parsed on demand. See {@link LazyGffFeature}. */
export function getAttribute(feature: LazyGffFeature, key: string) {
  return scanAttributes(
    feature.attributeString,
    primaryTag(key),
    alternateTag(key),
    undefined,
    undefined,
  ).value0
}

/**
 * Every attribute of a lazily-parsed feature, as the eager parser would have
 * spread them onto it. This is the full column-9 parse the lazy path defers,
 * so call it when a caller genuinely needs all of them — serializing a feature
 * for a detail view, say.
 */
export function getAttributes(feature: LazyGffFeature) {
  const result: Record<string, unknown> = {}
  // The eager path derives shouldUnescape from the whole line rather than just
  // column 9, but `unescape` returns its input untouched when there is no '%',
  // so narrowing it to the attribute string cannot change any value.
  parseAttributes(
    feature.attributeString,
    result,
    feature.attributeString.includes('%'),
  )
  return result
}

/**
 * ID and Parent together, the two attributes tree-building needs, in one pass.
 * Neither name is reserved and neither ends in '2', so both resolve to a plain
 * tag with no alternate.
 */
export function getLinkAttributes(feature: LazyGffFeature) {
  const { value0, value1 } = scanAttributes(
    feature.attributeString,
    'id',
    undefined,
    'parent',
    undefined,
  )
  return { id: value0, parent: value1 }
}

/**
 * Parse a GFF3 feature line, leaving column 9 as raw text. The
 * attribute-deferring counterpart to {@link parseFeature}; see
 * {@link LazyGffFeature} for when that pays.
 *
 * @param line - GFF3 feature line
 * @returns The parsed feature, attributes unparsed
 */
export function parseFeatureLazy(line: string): LazyGffFeature {
  const f = line.split('\t')
  const shouldUnescape = line.includes('%')
  const startStr = f[3]
  const endStr = f[4]
  const scoreStr = f[5]
  const strandStr = f[6]
  const phase = f[7]

  return {
    refName: strField(f[0], shouldUnescape, ''),
    source: strField(f[1], shouldUnescape, null),
    type: strField(f[2], shouldUnescape, null),
    start: isEmpty(startStr) ? 0 : +startStr - 1,
    end: isEmpty(endStr) ? 0 : +endStr,
    score: isEmpty(scoreStr) ? undefined : +scoreStr,
    strand: strandStr === undefined ? undefined : STRAND_MAP[strandStr],
    phase: isEmpty(phase) ? undefined : +phase,
    subfeatures: [],
    attributeString: f[8] ?? '',
  }
}

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
