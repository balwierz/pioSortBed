# LociSSD file format specification

**Format version:** 2
**Status:** Stable
**Last updated:** 2026-05-06
**Reference implementation:** [lociSSD](https://github.com/) Python package, v0.2.0+

This document specifies the on-disk format of a LociSSD dataset precisely
enough that independent implementations (in any language with Apache
Parquet bindings) can write files the reference reader accepts and read
files the reference writer produces.

---

## 1. Conformance terminology

The keywords **MUST**, **MUST NOT**, **SHOULD**, **SHOULD NOT**, **MAY**,
and **REQUIRED** are to be interpreted as described in
[RFC 2119](https://www.rfc-editor.org/rfc/rfc2119).

A *conforming writer* is software that produces files satisfying every
**MUST** in this document. A *conforming reader* is software that
correctly interprets every file produced by a conforming writer.

---

## 2. Container

A LociSSD dataset **MUST** be a single Apache Parquet file.

* The file **MUST** conform to Parquet format `2.6` or later
  ([Parquet spec](https://github.com/apache/parquet-format)).
* Writers **SHOULD** use `data_page_version = 2.0`.
* Writers **SHOULD** enable per-row-group statistics
  (`write_statistics = true`).
* Writers **SHOULD** enable the page index
  (`write_page_index = true`).

The file extension is conventionally `.lociss` but is not enforced; any
extension is acceptable.

---

## 3. Schema

### 3.1 Required columns

Every conforming dataset **MUST** contain the following columns. Names
are case-sensitive.

| Column | Arrow type | Notes |
|---|---|---|
| `Chromosome` | `string` or `large_string` | Stored as plain string on disk. Dictionary encoding via Parquet `use_dictionary=true` is permitted and recommended. |
| `Start` | `int32` (default) or `int64` | 0-based, inclusive — BED convention. The default is `int32`, sufficient for any per-chromosome position up to `INT32_MAX = 2,147,483,647`. Use `int64` only when the assembly actually needs it (e.g. axolotl, lily, onion — chromosomes whose length exceeds INT32_MAX). |
| `End` | matches `Start`'s integer width | 0-based, exclusive. **MUST** satisfy `End > Start` for every row. |
| `MaxEndSoFar` | matches `Start`'s integer width | Derived; see §4 |

The integer width chosen by the writer is recorded in the manifest's
``coord_dtype`` field (§6.3). Conforming readers **MUST** accept either
width and **SHOULD** check that ``Start`` / ``End`` / ``MaxEndSoFar`` all
share a single width per file.

### 3.2 Optional columns

A dataset **MAY** contain any number of additional columns of any Arrow
type. Common cases:

* `Strand` — `dictionary<string>` or `string`. Conventionally one of
  `"+"`, `"-"`, or `"."`.
* User payload columns of any nested Arrow type:
  `List`, `LargeList`, `FixedSizeList`, `Map`, `Struct`, `LargeString`,
  `LargeBinary`, `Decimal128`, `Decimal256`, `Timestamp`, etc.
* Pickle-encoded columns; see §3.4.

### 3.3 Column ordering

Column order on disk is not normative — readers **MUST** locate columns
by name. However, conforming writers **SHOULD** order columns as:
required loci columns first (`Chromosome`, `Start`, `End`), then
`Strand` if present, then user columns in user-supplied order, then
`MaxEndSoFar` last.

### 3.4 Pickle-encoded columns

Columns whose values cannot be represented as a native Arrow type **MAY**
be stored as Python-pickled bytes in a `LargeBinary` column. Such columns
**MUST** be tagged in the manifest with `encoding = "pickle"` (see §6).

Non-Python implementations cannot interpret the contents of pickle
columns and **SHOULD** expose them as opaque bytes. They **MUST NOT**
silently coerce or discard the bytes.

---

## 4. Sort order and `MaxEndSoFar`

### 4.1 Canonical sort

Rows **MUST** be sorted ascending by the lexicographic tuple
`(Chromosome, Start, End)`.

* Chromosome ordering uses the underlying string's byte ordering. The
  natural-sort convention common in genomics (`chr2 < chr11`) is **NOT**
  required by this spec; both `"chr10" < "chr2"` (lexicographic) and a
  natural-sort ordering produced by some other tool are conformant as
  long as the chosen ordering is consistent within the file and matches
  the ordering recorded in the manifest's `chromosomes` list (§6).
* Within a chromosome, rows are ordered by `Start` ascending; ties are
  broken by `End` ascending.
* Strand-aware sorting (descending `Start` for `-` strand, as used by
  some PyRanges tooling) is **NOT** the storage order. Implementations
  that need strand-aware order **MUST** re-sort on read.

A consequence: all rows for a given chromosome occupy a contiguous run
in row order.

### 4.2 `MaxEndSoFar` definition

`MaxEndSoFar` is a derived `int64` column that **MUST** satisfy the
following for every row index `i`:

```
MaxEndSoFar[i] = max(End[j] for all j in [run_start(i), i])
```

where `run_start(i)` is the index of the first row in the contiguous
chromosome run containing row `i`.

Equivalently: `MaxEndSoFar` is `cumulative_max(End)` computed
*independently within each chromosome run*, resetting to the run's first
`End` value at every chromosome boundary.

Within a chromosome run, `MaxEndSoFar` is therefore non-decreasing.
Across the chromosome boundary, `MaxEndSoFar` **MUST NOT** carry over
from the previous chromosome.

### 4.3 Sort hint in Parquet metadata

Conforming writers **SHOULD** emit Parquet `SortingColumn` hints for
`Chromosome`, `Start`, and `End` (all ascending, nulls last). These
hints are advisory and a missing or incorrect hint does not invalidate
an otherwise-conforming file.

---

## 5. Recommended Parquet writer settings

These settings are **RECOMMENDED**, not required. Any Parquet-conformant
combination of options yields a valid LociSSD file as long as the schema,
sort order, manifest, and `MaxEndSoFar` invariants hold.

| Setting | Recommended value |
|---|---|
| `version` | `"2.6"` |
| `data_page_version` | `"2.0"` |
| `use_dictionary` | `true` |
| `write_statistics` | `true` |
| `write_page_index` | `true` |
| Default compression | `zstd` level 3 |
| `row_group_size` | `65536` |

Per-column compression overrides are encouraged for heavy columns
(sequences, embeddings, signal vectors) — see §6 for how the manifest
records them.

---

## 6. The manifest

The manifest is a UTF-8 JSON document embedded in the Parquet file's
footer key-value metadata.

### 6.1 Location

The manifest **MUST** be stored under key `lociSSD_manifest` (raw bytes:
`b"lociSSD_manifest"`) in the Parquet file-level key-value metadata
(`FileMetaData.key_value_metadata` in the
[Parquet thrift schema](https://github.com/apache/parquet-format/blob/master/src/main/thrift/parquet.thrift)).

Conforming writers **MAY** also mirror the manifest into the Arrow
schema metadata (under the same key) — this happens automatically when
the writer attaches the metadata to the Arrow schema before opening
`ParquetWriter`. Conforming readers **MUST** read from the file-level KV
metadata location, since that location is always populated regardless of
write path.

### 6.2 Encoding

The value associated with the `lociSSD_manifest` key **MUST** be a UTF-8
encoded JSON document.

### 6.3 Manifest JSON schema

```jsonc
{
  // REQUIRED. Currently 2. Reader behavior on unknown versions: §9.
  "format_version": 2,

  // REQUIRED. Free-form identifier of the producing software.
  "writer_version": "lociSSD 0.2.0",

  // REQUIRED. ISO 8601 timestamp with timezone offset.
  "created_utc": "2026-05-06T12:34:56+00:00",

  // REQUIRED. Total number of rows. MUST equal the sum of `chromosomes[].rows`
  // and MUST equal the Parquet file's row count.
  "row_count": 12345678,

  // REQUIRED. One entry per chromosome present in the file, in storage order.
  // The storage order MUST match the order in which chromosomes first appear
  // when reading the Chromosome column from row 0.
  "chromosomes": [
    {
      "name": "chr1",            // REQUIRED. String value matching Chromosome column.
      "rows": 234567,            // REQUIRED. Row count for this chromosome.
      "row_offset": 0,           // REQUIRED. Global row index of this chromosome's first row.
      "min_start": 0,            // REQUIRED. min(Start) over rows in this chromosome.
      "max_start": 248956421,    // REQUIRED. max(Start) over rows in this chromosome.
      "max_end": 248956422       // REQUIRED. max(End) over rows in this chromosome.
    }
  ],

  // REQUIRED. One entry per column in the Parquet file's schema. Order is
  // not normative; readers MUST look up by `name`.
  "schema_columns": [
    {
      "name": "Chromosome",      // REQUIRED. Matches Parquet column name exactly.
      "arrow_type": "string",    // REQUIRED. Free-form Arrow type repr.
      "compression": ["zstd", 3],// REQUIRED. [codec_name, level_or_null].
      "encoding": "native",      // REQUIRED. "native" or "pickle". §3.4.
      "derived": false           // OPTIONAL. true for MaxEndSoFar, false otherwise.
    },
    {
      "name": "MaxEndSoFar",
      "arrow_type": "int64",
      "compression": ["zstd", 3],
      "encoding": "native",
      "derived": true
    }
  ],

  // REQUIRED. The columns the data is sorted by, in priority order.
  // For format_version 2, this MUST be exactly ["Chromosome", "Start", "End"].
  "sort_keys": ["Chromosome", "Start", "End"],

  // REQUIRED. Target row group size used when writing.
  "row_group_size": 65536,

  // REQUIRED. Parquet data page version used.
  "data_page_version": "2.0",

  // REQUIRED. Default codec applied where no per-column override exists.
  // Format: [codec_name, level_or_null]. Codec name MUST be a Parquet
  // compression codec ("zstd", "snappy", "gzip", "lz4", "brotli", "none").
  "default_compression": ["zstd", 3],

  // OPTIONAL reference-genome / provenance metadata. Added in lociSSD
  // 0.3.0 as additive fields; readers MUST tolerate either explicit null
  // or a missing key on older files.

  // OPTIONAL. Reference assembly identifier (e.g. "GRCh38", "hg19", "mm39").
  "assembly": "GRCh38",

  // OPTIONAL. Species name (e.g. "Homo sapiens").
  "species": "Homo sapiens",

  // OPTIONAL. Provenance URL of the source data.
  "source_url": "https://encode.org/.../regions.bed.gz",

  // OPTIONAL. Free-form user metadata. Keys MUST be strings; values
  // MUST be JSON-serializable (string, number, boolean, null, array,
  // or object). String-valued keys are recommended for cross-tool
  // portability.
  "user_metadata": {
    "cell_type": "K562",
    "experiment": "ATAC-seq"
  },

  // OPTIONAL. Integer width chosen for Start / End / MaxEndSoFar:
  // "int32" or "int64". Conforming writers SHOULD set this so a
  // reader can tell that an int64 file was a deliberate choice rather
  // than an oversight by another tool that didn't think about width.
  // Conforming readers MUST cross-check this against the actual
  // schema and reject mismatches (e.g. coord_dtype="int32" with
  // schema int64 columns is invalid). Legacy datasets without this
  // key infer the width from the schema.
  "coord_dtype": "int32"
}
```

### 6.4 Manifest validation rules (normative)

A conforming reader **MUST**:

1. Reject files whose `format_version` is unknown to it.
2. Cross-reference `chromosomes[*].name` with distinct values in the
   Chromosome column; mismatches are corruption.
3. Treat `chromosomes[*].rows` and `chromosomes[*].row_offset` as
   authoritative for chromosome locations within the file.

A conforming reader **SHOULD**:

4. Verify that `sum(chromosomes[*].rows) == row_count == Parquet file's
   row count`.
5. Recompute `MaxEndSoFar` on demand and compare against the stored
   column when validating files (this is a `O(N)` check; it's optional
   for normal reads but recommended for a `verify` operation).

A conforming reader **MAY**:

6. Use the manifest's per-chromosome stats to short-circuit queries
   without opening row groups.

---

## 6.5 Optional interval index (footer KV)

A conforming writer **MAY** embed a per-chromosome **augmented sorted
interval index** in the Parquet file's footer key-value metadata under
the key ``b"lociSSD_interval_index"``. (This is *not* NCList — see the
note at the end of this section.) The presence of this key is
**OPTIONAL**; conforming readers **MUST** fall back to the §7
predicate-pushdown algorithm when the key is absent or the embedded
blob fails to parse.

When present, the index value is a **base64-encoded** Arrow IPC stream
(itself zstd-compressed) encoding a single Arrow Table whose schema is:

```
chromosome:        string
starts:            large_list<int64>
ends:              large_list<int64>
max_end_running:   large_list<int64>
row_id:            large_list<int64>
```

with one row per chromosome, in storage order. Each list field has the
same length as the chromosome's row count; element ``i`` of the lists
describes the ``i``-th row of that chromosome (in canonical sort order).

Semantics:

* ``starts[i]`` and ``ends[i]`` mirror the ``Start`` and ``End`` columns
  of the underlying Parquet for that chromosome.
* ``max_end_running[i] = max(ends[0..i])`` per chromosome — the same
  invariant as the ``MaxEndSoFar`` data column, kept here so a reader
  doing range queries doesn't have to fetch the column from Parquet
  just to prune.
* ``row_id[i]`` is the dataset's global row index of that chromosome's
  ``i``-th row (i.e. ``ChromosomeSpec.row_offset + i``).

The Arrow IPC stream's schema-level metadata carries one extra key,
``b"lociSSD_index_format_version"``, holding the index format version
("2" today). Readers **MUST** reject blobs whose declared version is
not understood.

**Why base64?** Parquet KV-metadata values are nominally thrift
``string``. Several mainstream Parquet readers (notably polars 1.38)
strictly validate them as UTF-8 and reject arbitrary byte sequences
with "Invalid thrift: bad data." The base64 wrapper keeps the value
pure ASCII so any conforming Parquet reader can still parse the file
even if it does not interpret the index. Index format version 1 stored
the IPC stream raw; v2 wraps it in base64. Conforming readers
**MUST** support v2 (current); supporting v1 is optional since v1 was
only ever produced by a development build that was never released.

The reference implementation answers a region query
``[q_start, q_end)`` on chromosome ``c`` by binary-searching ``starts``
for ``q_end`` and ``max_end_running`` for ``q_start``, then walking the
narrow candidate slice testing ``ends[i] > q_start``. This is
``O(log N + k)`` (where k is the size of the candidate window) and
gives row-precision rather than the row-group-precision of the
predicate path.

**Note on naming.** Earlier drafts of this document called the structure
"NCList-style." That was misleading: NCList (Alekseyenko & Lee 2007) is
a hierarchical forest of nested lists with parent/child pointers, while
this format is a flat sorted array augmented with a cumulative-max-end
column. Implementers porting an existing real NCList implementation
should still produce the four arrays described above; the on-disk format
is independent of how a writer/reader chooses to query in memory.

### 6.5.1 Scaling limit and future relocations

The index lives in the Parquet **file-footer KV metadata**. pyarrow (and
every Parquet reader I've checked) parses the entire FileMetaData
thrift structure when opening a file — there's no API to fetch one KV
key without reading them all. At ~7 bytes per row in the compressed
index, this becomes a problem at large scales:

| Rows | Index size | Footer behavior |
|---:|---:|---|
| 1M | ~7 MB | trivial |
| 10M | ~70 MB | fine local, awkward over HTTP |
| 100M | ~700 MB | painful — even `info` opens are slow |
| 1B | ~7 GB | broken — opening reads multi-GB |

The format spec keeps the index **OPTIONAL**, so writers producing
TB-scale datasets **SHOULD** omit the index (write with
``interval_index=False`` in the reference implementation) and let
readers fall back to the §7 predicate-pushdown path. Predicate
pushdown's cost scales with the row-group count, not row count, so
properly tuned ``row_group_size`` keeps queries fast at TB scale.

A future format version **MAY** relocate the index to one of:

* **Dedicated Parquet row groups inside the same file**, classified
  via a new manifest field (e.g. ``index_row_groups: list[int]``).
  Readers use Parquet's standard column-projection to load only the
  chromosomes they query. Preserves the single-file property.

* **An external sidecar file** (e.g. ``regions.lociss.idx``) holding
  the same Arrow IPC payload. Loses the single-file operational win
  format v2 was designed around, but is the simplest fix.

Per-key isolation in KV metadata is **not** a viable mitigation —
pyarrow loads all KV bytes when parsing the footer, so splitting one
big blob into per-chromosome keys does not reduce open cost.

Both relocations are spec-additive (new manifest fields, new file
position) and would not change the index payload format described
above. For now, conforming readers **MUST** support reading the index
from footer KV metadata (when present) and **MUST** fall back to the
predicate path when absent.

## 7. Range query predicate (advisory)

Tools that want to take advantage of LociSSD's row-group pruning for
BED-overlap queries `(query_chromosome, query_start, query_end)` **SHOULD**
push down the following predicate to Parquet:

```
   Chromosome == query_chromosome
   AND Start < query_end
   AND End > query_start
   AND MaxEndSoFar > query_start
```

The fourth clause is logically redundant with the third (because
`MaxEndSoFar[i] >= End[i]` for any single row), but it lets a
statistics-aware Parquet reader prune row groups via min/max stats on
the *cumulative* `MaxEndSoFar` even when a few long intervals would
inflate `End`'s row-group max statistic and defeat pruning of the
underlying `End` column directly.

This predicate is the recommended canonical form. Other equivalent forms
are acceptable; the on-disk format does not depend on the choice.

---

## 8. Recommended atomic write

Conforming writers **SHOULD** write atomically: write the Parquet file
to a temporary path in the same directory as the target, `fsync` the
file, then `rename(tmp, target)`. POSIX `rename(2)` is atomic across the
target name.

This recommendation prevents partial files from appearing at the target
path if the writer crashes. Tools that cannot satisfy it (e.g. cloud
object stores without atomic rename) **SHOULD** document the resulting
visibility semantics.

---

## 9. Versioning

The current format version is `2`. The manifest's `format_version` field
**MUST** be `2` for files produced today.

* **Format version 1** used a chromosome-partitioned Hive directory
  layout with an external `_manifest.json` file. Format version 1 is
  **NOT** supported by current readers.
* **Future format versions** (`3+`) will be backwards-incompatible.
  Conforming readers **MUST** reject unknown `format_version` values
  rather than guessing.

A future minor version bump within `format_version: 2` is reserved for
*additive* manifest fields. Conforming readers **MUST** ignore unknown
keys at the top level of the manifest object and within
`chromosomes[*]` and `schema_columns[*]` entries, so that producers can
add metadata fields (e.g. `assembly`, `species`, custom tags) without
breaking existing readers.

---

## 10. Minimal example

A conforming dataset with two intervals on `chr1` and one on `chrX`,
just the required loci columns plus `Score`.

### Parquet schema (Arrow notation)

```
Chromosome:   string
Start:        int64
End:          int64
Score:        double
MaxEndSoFar:  int64
```

### Sample row data (in storage order)

| Chromosome | Start | End | Score | MaxEndSoFar |
|---|---:|---:|---:|---:|
| chr1 | 100 | 200 | 0.5 | 200 |
| chr1 | 150 | 1000 | 0.7 | 1000 |
| chrX | 50 | 80 | 0.9 | 80 |

Note `MaxEndSoFar[2] == 80`, **not** `1000` — it resets at the chr1→chrX
boundary.

### Manifest (UTF-8 JSON, stored under file-level KV `lociSSD_manifest`)

```json
{
  "format_version": 2,
  "writer_version": "minimal-writer 0.0.1",
  "created_utc": "2026-05-06T00:00:00+00:00",
  "row_count": 3,
  "chromosomes": [
    {"name": "chr1", "rows": 2, "row_offset": 0, "min_start": 100, "max_start": 150, "max_end": 1000},
    {"name": "chrX", "rows": 1, "row_offset": 2, "min_start": 50,  "max_start": 50,  "max_end": 80}
  ],
  "schema_columns": [
    {"name": "Chromosome",  "arrow_type": "string", "compression": ["zstd", 3], "encoding": "native", "derived": false},
    {"name": "Start",       "arrow_type": "int64",  "compression": ["zstd", 3], "encoding": "native", "derived": false},
    {"name": "End",         "arrow_type": "int64",  "compression": ["zstd", 3], "encoding": "native", "derived": false},
    {"name": "Score",       "arrow_type": "double", "compression": ["zstd", 3], "encoding": "native", "derived": false},
    {"name": "MaxEndSoFar", "arrow_type": "int64",  "compression": ["zstd", 3], "encoding": "native", "derived": true}
  ],
  "sort_keys": ["Chromosome", "Start", "End"],
  "row_group_size": 65536,
  "data_page_version": "2.0",
  "default_compression": ["zstd", 3]
}
```

---

## 11. Conformance checklist

A file is a conforming LociSSD v2 dataset if and only if:

- [ ] It is a valid Apache Parquet file (any version ≥ 2.6 recommended).
- [ ] Its schema contains `Chromosome`, `Start`, `End`, and `MaxEndSoFar`
      columns with the types specified in §3.1.
- [ ] Every row satisfies `End > Start` (BED half-open invariant).
- [ ] Rows are sorted ascending by `(Chromosome, Start, End)` (§4.1).
- [ ] All rows for each chromosome occupy a contiguous run.
- [ ] `MaxEndSoFar` satisfies the per-chromosome cumulative-max
      invariant of §4.2.
- [ ] The Parquet file-level KV metadata contains a `lociSSD_manifest`
      key whose value is a UTF-8 JSON document conforming to §6.3.
- [ ] The manifest's `format_version` is `2`.
- [ ] The manifest's `chromosomes` list, `row_count`, and `sort_keys`
      are consistent with the actual data (§6.4).

A file that satisfies all of the above is interoperable with any
conforming reader.

---

## 12. Appendix: differences from format version 1 (historical)

Format version 1 used a directory layout with hive-style chromosome
partitions and an external `_manifest.json` sidecar. It was superseded
in May 2026 by format version 2 (single-file Parquet) for operational
reasons — single-file storage is easier to copy, `scp`, and upload to
object stores, and removes the partial-write recovery question that
arises when only some partition files exist on disk.

The current reference implementation does not include a v1 reader.
Tools that need to migrate v1 datasets should re-emit them via the
reference writer, which produces v2 output.
