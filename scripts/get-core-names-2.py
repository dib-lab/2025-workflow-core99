#! /usr/bin/env python
"""
Get core names from manysearch results.
"""
import os
import sys
import argparse
import polars as pl

def filter_combined_df(combined_df, frac, threshold_bp, min_mbases):
    sub_df = combined_df.filter((pl.col("mbases") >= min_mbases))

    total = sub_df["metag"].n_unique()
    cutoff = int(frac * total + 1) # don't round the threshold down :sigh:

    group_df = sub_df.group_by('species_with_acc') \
            .agg(pl.len()) \
            .filter(pl.col("len") >= cutoff)

    return (cutoff, total, sub_df, group_df)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('manysearch_csvs', nargs='+')
    p.add_argument('-m', '--metadata-csv', required=True)
    p.add_argument('--fraction-threshold', type=float, default=0.95,
                   help='frequency cutoff for core (>=)')
    p.add_argument('-o', '--output', required=True, help='output CSV')
    p.add_argument('--save-acc-names', help="save acc+names to this file")
    p.add_argument('--save-only-names', help="save only names to this file")
    p.add_argument('--save-core-csv', help="save CSV lines for core to this file")
    p.add_argument('--expect-num-metagenomes', type=int)
    p.add_argument('--min-mbases', default=1000, type=int)
    p.add_argument('--min-intersect-threshold', default=20_000)
    args = p.parse_args()

    total_rows = 0
    dflist = []

    print(f"reading {len(args.manysearch_csvs)} manysearch parquet files total")
    for n, filename in enumerate(args.manysearch_csvs):
        if n % 1000 == 0 and n:
            print(f'...reading {n} of {len(args.manysearch_csvs)} parquets')
        df = (pl.scan_parquet(filename)
              .with_columns(
                  pl.col("match_name").alias("species_with_acc"),
                  pl.col("query_name").alias("metag"),
                  (pl.col("intersect_hashes") * pl.col("scaled")).alias("intersect_bp"),
              )
              .select(["species_with_acc", "metag", "intersect_bp"])
              .filter(pl.col("intersect_bp") >= args.min_intersect_threshold)
              .collect())

        dflist.append(df)
        total_rows += len(df)

    all_df = pl.concat(dflist, how='vertical')

    frac = args.fraction_threshold
    print(f"searching {len(all_df)} rows / {all_df['species_with_acc'].n_unique()} species / {all_df['metag'].n_unique()} metagenomes for core species at {frac}")
    print(f"using only datasets >= {args.min_mbases}MB and unique_intersect_bp >= {args.min_intersect_threshold}")

    meta_df = pl.read_csv(args.metadata_csv)
    meta_df = meta_df.select(['acc', 'mbases'])
    print(f"loaded {len(meta_df)} metadata rows from '{args.metadata_csv}'")

    combo_df = all_df.join(meta_df, left_on='metag', right_on='acc')
    assert len(combo_df) == len(all_df), "should NOT lose rows"

    if args.expect_num_metagenomes:
        assert combo_df['metag'].n_unique() == args.expect_num_metagenomes, combo_df['metag'].n_unique()
    else:
        print(f"WARNING: no number of expected metagenomes set, verify thyself")
    ###

    cutoff, total, filtered_df, final_df = filter_combined_df(combo_df, frac,
                                                 args.min_intersect_threshold,
                                                 args.min_mbases)

    print(f"total number of metags: {total}")
    print(f">= for frequency cutoff {frac}: {cutoff}")
 
    print(f"writing CSV to output '{args.output}'")
    final_df.write_csv(args.output)

    if args.save_acc_names:
        print(f"writing acc names to output '{args.save_acc_names}'")
        names = set(final_df["species_with_acc"].to_list())
        with open(args.save_acc_names, "wt") as fp:
            fp.write("\n".join(names))

    if args.save_only_names:
        print(f"writing names to output '{args.save_only_names}'")
        names = set(final_df["species_with_acc"].to_list())
        names = [ ' '.join(n.split(' ')[1:3]) for n in names ]
        names.sort()
        with open(args.save_only_names, "wt") as fp:
            fp.write("\n".join(names))

    if args.save_core_csv:
        print(f"writing core CSV lines to output '{args.save_core_csv}'")
        names = set(final_df["species_with_acc"].to_list())

        core_df = filtered_df.filter(pl.col("species_with_acc").is_in(names))

        core_df.write_csv(args.save_core_csv)


if __name__ == '__main__':
    sys.exit(main())
