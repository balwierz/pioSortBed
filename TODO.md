# Deferred items

Tracked separately from the main commit log so they don't get lost.
None of these block existing functionality.

## LociSSD output (`--lociss-output`)

- **Tail passthrough** — currently we drop BED columns 4+ (Name / Score
  / Strand / etc.) with a one-time stderr warning. Spec allows
  arbitrary user columns; we should:
  - Detect column count once on the first record.
  - For BED4: add `Name` (string).
  - For BED6: add `Score` (string — BED scores can be non-numeric like
    `.`) and `Strand` (`dictionary<string>` per spec).
  - For BED12: add the six remaining columns with sensible Arrow types.
  - Anything past the recognised columns: catch-all `Tail` string
    column.

- **§6.5 optional interval index** — Arrow IPC stream embedded under
  footer KV `lociSSD_interval_index`. Optional in the spec; readers
  fall back to predicate pushdown if absent. Worth adding if/when a
  reader actually needs row-precision pruning instead of row-group
  precision.

- **`--collapse` + `--lociss-output`** — currently rejected at CLI
  parse. Could emit the collapsed weight as a `Score` column
  (`double` or `int32`). Needs design clarity on what column name
  the spec expects.

- **`int64` coord mode** — we always write `int32` Start/End/MaxEndSoFar.
  Spec allows int64 for assemblies with >2 Gbp chromosomes (axolotl,
  lily, onion). Add a `--lociss-int64` flag when a real >2 Gbp dataset
  comes up.

- **Parquet `SortingColumn` hint** — spec §4.3 calls it advisory.
  Currently not emitted. ~5 lines to add.

## General

- **External-merge: parallelise zstd compression in the writer** —
  currently zstd MT was disabled because it regressed at 1 MiB block
  sizes. Bigger blocks (4–8 MiB) might re-enable it usefully.
