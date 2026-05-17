# Deferred items

Tracked separately from the main commit log so they don't get lost.
None of these block existing functionality.

## LociSSD output (`--lociss-output`)

- **Typed BED4/6/12 columns** — Tail passthrough is implemented as a
  single catch-all `Tail` string column (raw tab-separated post-End
  bytes; preserved across all four sort paths). This is the
  FORMAT_SPEC §3.3-compliant fallback. A nicer experience for
  downstream tools would be to detect BED4 / BED5 / BED6 / BED12 on
  the first record and split into typed Arrow columns: Name (string),
  Score (string — BED scores can be non-numeric like `.`), Strand
  (dictionary&lt;string&gt;), ThickStart/ThickEnd (int32), ItemRgb
  (string), BlockCount (int32), BlockSizes / BlockStarts (string).
  Downstream tools could then push predicates on Strand etc.; today
  they have to `.str.split("\t")` the Tail column themselves. Same
  parser work in the writer regardless; mostly a question of when
  to do it.

- **`--collapse` + `--lociss-output`** — currently rejected at CLI
  parse. Could emit the collapsed weight as a `Score` column
  (`double` or `int32`). Needs design clarity on what column name
  the spec expects.

- **`int64` coord mode** — we always write `int32` Start/End/MaxEndSoFar.
  Spec allows int64 for assemblies with >2 Gbp chromosomes (axolotl,
  lily, onion). Add a `--lociss-int64` flag when a real >2 Gbp dataset
  comes up.

## General

- **External-merge: parallelise zstd compression in the writer** —
  currently zstd MT was disabled because it regressed at 1 MiB block
  sizes. Bigger blocks (4–8 MiB) might re-enable it usefully.
