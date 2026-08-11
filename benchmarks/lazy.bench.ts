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
 * The fixture matters more than the timing method, and an earlier version of
 * this comment was wrong about the size of the effect because of it. Do NOT
 * size a GFF3 fixture up by concatenating a real file: ids repeat across the
 * copies, so every duplicate child attaches to the first parent carrying its
 * id. Twelve copies of a GENCODE excerpt gave 233 subfeatures per top-level
 * feature, against 12 in real tair10_chr1. That deepens the tree rather than
 * lengthening the file, and depth flatters the lazy side specifically, because
 * consumers that wrap features tend to copy each subfeature's attributes on the
 * way. It reported 3.03x end-to-end where a corpus of distinct genes says
 * 1.29x, and 8.5x on retained heap where the honest figure is 2.3x.
 *
 * Measured against a generated corpus of distinct genes (jb2bench
 * `ecosystem/results/gff3-lazy.md`, one process per arm, which is what this
 * vitest file cannot do):
 *
 *   parse only          ~1.2 attrs/line  1.2-1.5x   ~10 attrs/line  1.9-2.2x
 *   parse + the reads a render performs  1.0-1.2x                   1.0-1.2x
 *
 * The second row is the one to keep in mind: at the point of use the lazy side
 * gives most of the parse win back, because every getAttribute is a scan and a
 * render does several per feature. Deferring work only pays if the consumer
 * does not then pay for it on the way out.
 *
 * A consumer that also stops materializing attributes elsewhere does better —
 * JBrowse measures 1.29x/1.36x end to end, but that comes from its Feature
 * wrapper no longer spreading every subfeature, not from this library.
 *
 * Retained heap after parsing, which for a caller holding a whole file resident
 * is the more durable effect: on ~10 attrs/line, 16.1MB against 7.0MB (2.3x);
 * on ~1.2, 8.0MB against 7.1MB.
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
