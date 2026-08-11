## [5.2.1](https://github.com/GMOD/gff-nostream/compare/v5.2.0...v5.2.1) (2026-08-11)

### Performance Improvements

- Read the lazy parser's columns by scanning for tabs ([4d65c65](https://github.com/GMOD/gff-nostream/commit/4d65c655902f93ec62f0ef0c1cf1061f5af0d33a))

## [5.2.0](https://github.com/GMOD/gff-nostream/compare/v5.1.3...v5.2.0) (2026-08-11)

### Chores

- Render only the commit subject, and link the commit ([1a790ca](https://github.com/GMOD/gff-nostream/commit/1a790ca7ed1a75febc894e3ea42b67b92f1f4527))
- Create a GitHub release for each published tag ([0550173](https://github.com/GMOD/gff-nostream/commit/05501738cfd7e163b84566b0d080f870786dcf31))
- Enforce type strippability in tsconfig ([dc1c5f3](https://github.com/GMOD/gff-nostream/commit/dc1c5f31941f3003c9afea6c7bd21cddc9fec809))

### Features

- Parse GFF3 attributes on demand ([7459c2a](https://github.com/GMOD/gff-nostream/commit/7459c2ad178ee64afaf86feb6bbaa1d371cd9579))

## [5.1.3](https://github.com/GMOD/gff-nostream/compare/v5.1.2...v5.1.3) (2026-08-10)

### Chores

- Drop prepublishOnly and the redundant pull_request workflow
- Add git-cliff for changelog generation
- Type-check the tests and enforce prettier, as @gmod/bam does
- Let npm publish stop auto-correcting repository.url
- Exempt our own packages from the release quarantine
- Bump pnpm/action-setup to v6.0.10
- Run the test suite as `pnpm test --run`
- Gate preversion on format:check, as CI does
- Gate preversion on typecheck too, as CI does
- Converge package.json on the shape its siblings use

### Documentation

- Mark breaking changes in the generated changelog

### Other Changes

- Backfill CHANGELOG.md for v0.0.1 through v5.1.2 ([5e18083](https://github.com/GMOD/gff-nostream/commit/5e1808354a158ead51ab08dc97d990f0a3e706b4))
- Revert "chore: converge package.json" — the CHANGELOG prettier step ([6aeae4c](https://github.com/GMOD/gff-nostream/commit/6aeae4c017b63958b770b36d7655fa663895c94d))

### Tests

- Pin the shared-ID linking path in the shapes it broke in
- Correct the note -- the regression never reached a release

# Unreleased

- Drop the redundant `prepublishOnly` script and `pull_request.yml` workflow: `preversion` already gates local publishes and the push workflow's `test` job already gates CI, so both were a third redundant run. `pull_request.yml` only added coverage for fork PRs, which this project's CI convention deliberately excludes. No user-facing change.

# v5.1.2 (2026-07-26)

- Add `parseLines(lines: readonly string[]): GffFeature[]`, a new entry point for callers that already have raw feature lines split out and have no per-line identity to carry through `parseRecords`' `{line}` wrapper (e.g. JBrowse's plain-GFF3 adapter, which groups lines by reference sequence before parsing)

# v5.1.1 (2026-07-25)

- Internal maintenance only: pin pnpm version, declare `sideEffects`, sha-pin CI actions, move CI to Node 24

# v5.1.0 (2026-07-25)

- Features whose `Parent` is never defined in the input (common when parsing a slice of a file, e.g. a tabix region query that cuts off the parent line) are now returned as top-level items instead of being silently dropped
- Percent-escaped multi-byte UTF-8 (e.g. `%E3%81%82`) now decodes correctly via `TextDecoder` instead of one Latin-1 character per byte; mixed-case escapes like `%eF` also decode correctly
- Lines with fewer than 9 columns no longer throw; missing trailing columns are treated the same as `.`
- An attribute named `Subfeatures` no longer replaces the `subfeatures` child array and breaks tree building — reserved field names are now suffixed with `2` on collision
- `extractType` returns the full type on a line with no fourth column, and `''` when there's no third column, instead of dropping the last character
- An attribute with an empty value list (`Foo=,,`) no longer leaves an empty array on the parsed feature

# v5.0.0 (2026-07-01)

- BREAKING: `parseRecords` now returns `ParsedRecord<R>[]` (`{ feature, record }` pairs) instead of `GffFeature[]`. Callers derive their own stable per-feature id (e.g. a tabix byte offset) from the record they passed in, instead of the parser stamping an opaque `_lineHash` onto the feature
- BREAKING: `LineRecord` is now `{ line: string }` only — `lineHash`, `start`, `end`, `type`, and `hasEscapes` are gone, and `GffFeature` no longer carries `_lineHash`
- `parseStringSync` gets a direct, allocation-free parse loop instead of building an intermediate records array

# v4.0.0 (2026-06-26)

- BREAKING: `parseStringSync`/`parseRecords` now return the flat `GffFeature[]` format (0-based half-open coordinates, numeric `strand`, lowercased attributes, `subfeatures`) that was previously only available via the `*JBrowse` variants. `parseStringSyncJBrowse`, `parseRecordsJBrowse`, and `JBrowseFeature` are removed — use `parseStringSync`/`parseRecords`/`GffFeature` instead
- Fix: a multi-location child (e.g. a discontinuous CDS/exon spanning several segments) is now attached to its parent's `subfeatures` once, instead of once per line duplicating the feature
- Fix: an invalid `%`-escape no longer swallows a following valid escape (e.g. `a%b%20c` now correctly unescapes to `a%b c`)
- Drop dead code: `parseDirective`, the `GFF3FeatureLine` type vocabulary, and the unescape/no-unescape wrapper function pairs

# v3.0.11 (2026-06-22)

- Fix: `parseRecordsJBrowse` dropped every line but the first of a parentless, multi-location top-level feature (e.g. `cDNA_match`/`EST_match` sharing one ID with no `Parent`) — every segment is now kept as its own top-level item

# v3.0.10 (2026-06-01)

- Fix: `parseRecordsJBrowse` dropped continuation lines of a multi-location child feature (e.g. a CDS split across several segments sharing one ID) from its parent's `subfeatures` — all segments are now attached

# v3.0.9 (2026-05-19)

- `LineRecord` now requires `start`, `end`, and `type` fields, all available from a tabix callback, so consumers no longer need intersection-type workarounds to carry them
- Add `extractType(line: string): string`, exported for fast column-3 (feature type) extraction without splitting the whole line

# v3.0.8 (2026-05-18)

- Internal only: rename the merged CI workflow back to `publish.yml` so npm OIDC trusted publishing still matches the expected workflow filename

# v3.0.7 (2026-05-18)

- Internal only: merge the publish workflow into the push workflow, gated on the test job succeeding; README CI badge updates

# v3.0.6 (2026-05-18)

- Fix: `parseRecordsJBrowse` no longer processes parent relationships on a discarded duplicate-ID feature, which had been creating dangling subfeature references
- Fix: `parseRecords` now shares `child_features`/`derived_features` across all lines of a multi-location feature (same ID on multiple lines), so children stay visible on every line regardless of arrival order — previously children added before the duplicate-ID line were visible only on the first line, while those added after were duplicated across both

# v3.0.5 (2026-04-27)

- Internal maintenance only: enable `noUncheckedIndexedAccess` in `tsconfig.json`, switch `eslint-plugin-import` to `eslint-plugin-import-x`, standardize the build pipeline on pnpm scripts and `rimraf`

# v3.0.4 (2026-04-27)

- Repo moved to the GMOD org (`GMOD/gff-nostream`); README and package.json updated accordingly
- Simplified the `exports` map in package.json to flat `import`/`require` entries instead of redundant nested conditions
- Dropped the unused `@jbrowse/quick-lru` dependency and the `documentation` doc-generation tooling
- Internal maintenance only otherwise: eslint bumped to v10, dep bumps, README cleanup

# v3.0.3 (2026-03-28)

- Migrated build tooling from yarn to pnpm and upgraded to TypeScript 6 with `nodenext` module resolution
- Added a GitHub Actions publish workflow using npm trusted publishing (OIDC, no stored token)
- Internal refactor: deduplicated the escaped/non-escaped variants of `parseFeature`, `parseFieldsArray`, and their JBrowse counterparts into shared implementations taking a `shouldUnescape` flag; no behavior change
- Added a JBrowse-format parsing benchmark suite; README rewritten for ESM usage, added CONTRIBUTING.md

# v3.0.2 (2025-12-24)

- Fixed `unescape()` to tolerate malformed or truncated `%XX` percent-escapes by decoding via a manual hex lookup table instead of `decodeURIComponent`, which throws on invalid sequences
- Added test coverage for the JBrowse-format parsing path (introduced in v3.0.1)

# v3.0.1 (2025-12-24)

- Added `parseStringSyncJBrowse` and `parseRecordsJBrowse`, which parse GFF3 directly into a flattened JBrowse feature shape (0-based `start`, numeric `strand`/`phase`, lowercased attribute keys unpacked onto the feature, `subfeatures` array for parent/child nesting) instead of the GFF3 array-of-lines format
- New `JBrowseFeature` type exported alongside the existing GFF3 types

# v3.0.0 (2025-12-24)

- BREAKING: `LineRecord` (used by `parseRecords`) changed shape from `{ fields: string[], lineHash? }` to `{ line: string, start, end, hasEscapes, lineHash? }` — callers now pass the raw tab-delimited line and an `hasEscapes` flag instead of pre-split fields
- BREAKING: `parseArraySync` and `parseRecordsSyncFast` removed; consolidated into a single `parseRecords(records: LineRecord[])`, with the escaped/non-escaped choice now made per-record via `hasEscapes`
- BREAKING: replaced the stateful `Parser` class (`src/parse.ts`) with a flat `Map`-based implementation in `api.ts`; as a result, multi-line features whose lines disagree on `type` are no longer rejected, and a feature with a `Parent` that never resolves is now silently dropped from the output instead of throwing a parse error
- Removed the formatting/escaping helpers (`formatFeature`, `formatAttributes`, `formatDirective`, `escape`, `escapeColumn`, etc.) from `src/util.ts` — the library is parse-only; these were already unreachable from the package's public entry point
- Performance: attribute and feature-line parsing rewritten to index-based scanning instead of `split`/`map`, reducing intermediate allocations

# v2.0.1 (2025-12-23)

- Added `parseRecordsSyncFast(records, hasEscapes)`, a fast path that skips per-value unescaping entirely when the caller knows the input has no percent-encoded characters

# v2.0.0 (2025-11-25)

- BREAKING: `##FASTA` sections are no longer parsed into sequence objects; the parser now stops at the FASTA marker instead of emitting `GFF3Sequence` items (FASTA support is fully removed, not just undocumented)
- BREAKING: package.json switched to `"type": "module"` with a proper conditional `exports` map, replacing the old `main`/`module` fields
- Added `parseArraySync` (parse an array of GFF3 line strings) and `parseRecordsSync` (parse an array of pre-tokenized `LineRecord` objects with fields already split and an optional `lineHash` to attach), for callers that already have GFF3 data split into lines or fields
- Performance: precompiled regexes, `unescape` fast-paths strings with no `%`, attribute/feature parsing avoids redundant work

# v1.3.9 (2025-09-29)

- Internal maintenance only: build scripts switched from npm to yarn, no source changes

# v1.3.7 (2025-09-29)

- Fixed a correctness bug for large files: the parser's internal buffer-eviction limit (previously defaulting to 1000 lines) is now unbounded by default, so features referenced by a `Parent`/`ID` more than 1000 lines apart are no longer prematurely evicted and mismatched

# v1.3.6 (2025-05-13)

- Internal maintenance only: dependency bumps, no source changes

# v1.3.5 (2025-05-13)

- Added a postbuild step that writes `dist/package.json` with `{"type": "commonjs"}`, fixing CJS consumers that were otherwise resolving the CJS build's `.js` files as ESM
- Dependency bumps (vitest 2 to 3, etc.)

# v1.3.4 (2024-09-04)

- Internal cleanup: removed the now-dead `ParseOptions`/`ParseOptionsProcessed` types and their unused imports left over from the v1.3.3 simplification; no behavior change

# v1.3.3 (2024-09-04)

- BREAKING: `parseStringSync` no longer accepts an options argument — it always parses only features and always resolves with `disableDerivesFromReferences: true`. Parsing of directives, comments, and sequences via options is no longer possible

# v1.3.1 (2024-09-04)

- BREAKING: package renamed from `@gmod/gff` to `gff-nostream`. The entire Node.js stream-based API (`parseStream`, `formatStream`, `formatFile`, `formatSync`, the `GFFTransform` class) and the `gff-to-json` CLI binary are removed; only synchronous string parsing remains
- BREAKING: `parseStringSync` signature simplified to a single overload returning `GFF3Item[]`, replacing the large set of conditional-type overloads keyed off `parseFeatures`/`parseDirectives`/`parseComments`/`parseSequences`. The `bufferSize` option is also removed from the public API
- BREAKING: default export removed; only named exports (`parseStringSync`, types) remain

# v1.3.0 (2022-12-06)

- Added a `stream-browserify` polyfill dependency so the package builds under webpack 5, which no longer auto-polyfills Node's `stream` module in the browser

# v1.2.3 (2022-10-18)

- Added `disableDerivesFromReferences` option to `parseStringSync`, allowing callers to skip resolving `Derives_from` references

# v1.2.2 (2022-10-18)

- Add `disableDerivesFromReferences` parse option to skip resolving `Derives_from` links, for files where that reference graph is large or unwanted
- Shorten the "features reference other features that do not exist" parse error to list the missing IDs instead of dumping the full orphan-tracking structure via `JSON.stringify`
- Switch the release script from `postpublish` to `postversion`; no runtime change

# v1.2.1 (2022-03-30)

- Publish the `src` directory in the npm package so source maps in `dist`/`esm` resolve correctly for consumers
- Bump dependencies

# v1.2.0 (2021-12-15)

- Rewrite the codebase in TypeScript; the published package now ships type declarations plus an ESM build (`esm/`) alongside the existing CommonJS `dist/`
- BREAKING: remove `parseFile` — it wrapped `fs.createReadStream(...).pipe(gff.parseStream(options))`; call that directly instead
- Fix `##sequence-region` directive parsing throwing when the line is missing start/end coordinates
- Move CI from Travis to GitHub Actions
- Internal maintenance otherwise: large batch of dependabot dependency bumps accumulated since v1.1.2 (2018)

# v1.1.2 (2018-11-12)

- Trim whitespace from attribute names and values when parsing the 9th column, so lines with stray spaces around `key=value` pairs parse cleanly instead of preserving the whitespace
- `escape`/`unescape` in `util.js` now coerce non-string input via `String(s)` instead of throwing when passed a non-string value

# v1.1.0 (2018-04-05)

- No user-facing API change: broadened the Babel build target to also compile for common browsers (previously targeted Node 6 only), fixing consumers who bundle this package for browser use

# v1.0.0 (2018-04-05)

- Column-value escaping now uses a narrower character set (`escapeColumn`) than attribute-value escaping (`escape`), so `;`, `=`, `&`, and `,` are no longer escaped in plain column values (only in the 9th-column attributes) — fixes over-escaping of `seqid`, `source`, `type`, etc.
- The `type` column (column 3) is now unescaped like `seqid` and `source` are
- `parseStringSync('')` now returns `[]` instead of throwing/erroring on empty input
- Avoids `require('fs')` when running under webpack, so browser bundles no longer fail on a missing `fs` module
- `src/index.js` now re-exports the individual named functions (`parseStream`, `parseFile`, `parseStringSync`, `formatSync`, `formatStream`, `formatFile`, `util`) as the default export instead of spreading the whole `api` module — cosmetic internally, but makes the exported shape explicit

# v0.0.1 (2018-04-03)

- Initial release of `@gmod/gff`, a streaming GFF3 parser/formatter for Node and the browser
- Core API: `parseStream`/`parseFile` (streaming parse, optionally including directives and comments via `parseAll`), `parseStringSync` (synchronous parse of a string), and `formatSync`/`formatStream`/`formatFile` for writing GFF3 back out
- Handles multi-location features and multiple parents, reconstructs feature hierarchies via `Parent` and `Derived_from`, and supports both implicit and explicit FASTA sections
- Proper escaping/unescaping of attribute values per the GFF3 spec
