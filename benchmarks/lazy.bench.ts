import fs from 'fs'
import { bench, describe } from 'vitest'

import { getAttribute, parseLines, parseLinesLazy } from '../src/index.ts'

/*
 * Eager vs lazy attribute parsing.
 *
 * Read the numbers here with care: vitest bench runs all iterations of one
 * `bench` before the next, so whichever runs second inherits a heap the first
 * one filled, and the eager parser fills a lot of it. Measured that way the
 * lazy path looked 9.6x faster on attribute-heavy input; interleaving the two
 * and taking per-run medians put the real figure at 2.7x. Treat a difference
 * under ~1.5x here as noise, and interleave before believing anything.
 *
 * Fairly measured (medians of 25 interleaved runs, node 24):
 *
 *   tair10_chr1   1.8 attrs/line   eager 446ms  lazy 469ms  0.95x (break-even)
 *   gencode-like 16.6 attrs/line   eager 529ms  lazy 195ms  2.71x
 *
 * Retained heap after parsing, which for a caller holding a whole file resident
 * is the bigger effect:
 *
 *   tair10_chr1                    eager 60.7MB  lazy 42.0MB  1.4x smaller
 *   gencode-like                   eager 58.4MB  lazy  6.9MB  8.5x smaller
 *
 * So the trade is: no faster when column 9 is nearly empty, materially faster
 * and much smaller when it is not. Annotation-grade GFF3 is the latter.
 */

function featureLines(path: string) {
  return fs
    .readFileSync(path, 'utf8')
    .split(/\r?\n/)
    .filter(line => line.length !== 0 && !line.startsWith('#'))
}

const opts = { iterations: 30, warmupIterations: 5 }

function benchShape(name: string, lines: string[]) {
  const attrs =
    lines.reduce(
      (a, l) => a + l.slice(l.lastIndexOf('\t') + 1).split(';').length,
      0,
    ) / lines.length

  describe(`${name} — ${lines.length} lines, ${attrs.toFixed(1)} attrs/line`, () => {
    bench('parseLines (eager)', () => void parseLines(lines), opts)
    bench('parseLinesLazy', () => void parseLinesLazy(lines), opts)

    // What a render actually does: the tree, plus the handful of attributes a
    // default-configured JBrowse feature track reads per top-level feature —
    // `gbkey` for the stock NCBI source-record filter and `name` for the label.
    // If the lazy win did not survive these it would not be worth having.
    bench(
      'parseLinesLazy + gbkey/name per feature',
      () => {
        for (const f of parseLinesLazy(lines)) {
          getAttribute(f, 'gbkey')
          getAttribute(f, 'name')
        }
      },
      opts,
    )
  })
}

benchShape('tair10_chr1', featureLines('test/data/tair10_chr1.gff'))
benchShape(
  'mm9_ensembl x40',
  Array.from({ length: 40 }, () =>
    featureLines('test/data/mm9_sample_ensembl.gff3'),
  ).flat(),
)
