#!/usr/bin/env python3
"""Region-query benchmark drivers for the bench_query.sh sweep.

Reads a query workload (chrom, start, end) from stdin (one query per line,
tab-separated) and a (reader, file) tuple from argv. Prints per-query
latencies in microseconds on stdout, one per line.

This is invoked once per (reader × cache state × region size) cell from
bench_query.sh; that script aggregates the median / p99 / throughput and
emits the CSV. Keeping all four readers in one process script means
Python startup cost is paid once per cell (not 1000x), and we don't
double-count import time as latency.

Readers covered:
  - tabix         (htslib CLI, via pysam)
  - polars        (pl.scan_parquet + predicate pushdown)
  - duckdb        (DuckDB SQL on the .lociss file)
  - loci          (loci.open_lociss + .region(); uses interval index if
                   --use-index is set, else predicate pushdown)
"""
from __future__ import annotations
import argparse
import sys
import time


def run_queries(reader: str, path: str, queries, use_index: bool = True):
    """Run all queries and yield per-query latency in microseconds."""
    if reader == "tabix":
        import pysam
        tbx = pysam.TabixFile(path)
        for chrom, start, end in queries:
            t0 = time.perf_counter_ns()
            n = 0
            for _ in tbx.fetch(chrom, start, end):
                n += 1
            t1 = time.perf_counter_ns()
            yield (t1 - t0) // 1000

    elif reader == "polars":
        import polars as pl
        # Use lazy frame + filter so predicate pushdown applies. The
        # MaxEndSoFar clause is logically redundant with End > q_start
        # but lets row-group stats prune more aggressively (see paper).
        for chrom, start, end in queries:
            t0 = time.perf_counter_ns()
            df = (
                pl.scan_parquet(path)
                .filter(
                    (pl.col("Chromosome") == chrom)
                    & (pl.col("Start") < end)
                    & (pl.col("End") > start)
                    & (pl.col("MaxEndSoFar") > start)
                )
                .collect()
            )
            t1 = time.perf_counter_ns()
            _ = df.height  # touch result so polars doesn't elide work
            yield (t1 - t0) // 1000

    elif reader == "duckdb":
        import duckdb
        # Re-open the connection per call would be unrealistic; in
        # practice you keep one connection. Match that.
        con = duckdb.connect(":memory:")
        # End is a SQL reserved word; quote it everywhere.
        # Use read_parquet() explicitly: .lociss is just our extension
        # for a Parquet file, but DuckDB's auto-format-detection only
        # looks at .parquet.
        sql = (
            f"SELECT Chromosome, Start, \"End\" FROM read_parquet('{path}') "
            "WHERE Chromosome = ? AND Start < ? AND \"End\" > ? "
            "AND MaxEndSoFar > ?"
        )
        for chrom, start, end in queries:
            t0 = time.perf_counter_ns()
            rows = con.execute(sql, [chrom, end, start, start]).fetchall()
            t1 = time.perf_counter_ns()
            _ = len(rows)
            yield (t1 - t0) // 1000

    elif reader == "loci":
        import loci
        ds = loci.open_lociss(path)
        for chrom, start, end in queries:
            t0 = time.perf_counter_ns()
            df = ds.region(chrom, start, end, use_index=use_index)
            t1 = time.perf_counter_ns()
            _ = len(df)
            yield (t1 - t0) // 1000

    else:
        raise SystemExit(f"unknown reader: {reader!r}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("reader", choices=["tabix", "polars", "duckdb", "loci"])
    ap.add_argument("path", help="indexed file path (.bed.gz for tabix, "
                                 ".lociss for the rest)")
    ap.add_argument("--use-index", action="store_true",
                    help="loci only: use embedded interval index (default off)")
    args = ap.parse_args()

    queries = []
    for line in sys.stdin:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        queries.append((parts[0], int(parts[1]), int(parts[2])))

    for us in run_queries(args.reader, args.path, queries,
                          use_index=args.use_index):
        print(us)


if __name__ == "__main__":
    main()
