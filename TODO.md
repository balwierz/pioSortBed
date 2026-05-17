# Deferred items

Tracked separately from the main commit log so they don't get lost.
None of these block existing functionality.

## LociSSD output (`--lociss-output`)

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
