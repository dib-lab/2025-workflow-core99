#! /usr/bin/env python
import os
import sys
import argparse
import polars as pl

def filter_combined_df(combined_df, frac, threshold_bp, min_mbases):
    sub_df = combined_df.filter(((pl.col("mbases") >= min_mbases) & 
                                (pl.col("unique_intersect_bp") >= threshold_bp)))

    total = sub_df["query_name"].n_unique()
    cutoff = int(frac * total + 1) # don't round the threshold down :sigh:

    group_df = sub_df.group_by('match_name') \
            .agg(pl.len()) \
            .filter(pl.col("len") >= cutoff)

    return (cutoff, total, group_df)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('gather_csvs', nargs='+')
    p.add_argument('-m', '--metadata-csv', required=True)
    p.add_argument('--match-name-col', default='name')
    p.add_argument('--fraction-threshold', type=float, default=0.99,
                   help='frequency cutoff for core (>=)')
    p.add_argument('-o', '--output', required=True, help='output CSV')
    p.add_argument('--save-acc-names', help="save acc+names to this file")
    p.add_argument('--save-only-names', help="save only names to this file")
    p.add_argument('--expect-num-metagenomes', type=int)
    p.add_argument('--min-mbases', default=1000, type=int)
    p.add_argument('--min-intersect-threshold', default=10_000)
    args = p.parse_args()

    total_rows = 0
    dflist = []

    print(f"reading {len(args.gather_csvs)} gather parquet files total")
    for n, csvfile in enumerate(args.gather_csvs):
        if n % 1000 == 0 and n:
            print(f'...reading {n} of {len(args.gather_csvs)} parquets')
        df = pl.read_parquet(csvfile)
        # don't pfaff about; support both gather and fast(multi)gather results
        if args.match_name_col != 'match_name':
            df = df.with_columns([pl.col(args.match_name_col).alias("match_name")])
        df = df.select(["match_name", "f_match_orig", "query_name", "unique_intersect_bp", "scaled"])
        
        dflist.append(df)
        total_rows += len(df)

    all_df = pl.concat(dflist, how='vertical')

    frac = args.fraction_threshold
    print(f"searching {len(all_df)} rows / {all_df['match_name'].n_unique()} species / {all_df['query_name'].n_unique()} metagenomes for core species at {frac}")
    print(f"using only datasets >= {args.min_mbases}MB and unique_intersect_bp >= {args.min_intersect_threshold}")

    meta_df = pl.read_csv(args.metadata_csv)
    meta_df = meta_df.select(['acc', 'mbases'])
    print(f"loaded {len(meta_df)} metadata rows from '{args.metadata_csv}'")

    combo_df = all_df.join(meta_df, left_on='query_name', right_on='acc')
    assert len(combo_df) == len(all_df), "should NOT lose rows"

    if args.expect_num_metagenomes:
        assert len(combo_df) == args.expect_num_metagenomes
    else:
        print(f"WARNING: no number of expected metagenomes set, verify thyself")
    ###

    cutoff, total, final_df = filter_combined_df(combo_df, frac,
                                                 args.min_intersect_threshold,
                                                 args.min_mbases)

    print(f"total data set number: {total}")
    print(f">= for frequency cutoff {frac}: {cutoff}")
 
    print(f"writing CSV to output '{args.output}'")
    final_df.write_csv(args.output)

    if args.save_acc_names:
        print(f"writing acc names to output '{args.save_acc_names}'")
        names = set(final_df["match_name"].to_list())
        with open(args.save_acc_names, "wt") as fp:
            fp.write("\n".join(names))

    if args.save_only_names:
        print(f"writing names to output '{args.save_only_names}'")
        names = set(final_df["match_name"].to_list())
        names = [ n.split(' ', 1)[1] for n in names ]
        names.sort()
        with open(args.save_only_names, "wt") as fp:
            fp.write("\n".join(names))


if __name__ == '__main__':
    sys.exit(main())
