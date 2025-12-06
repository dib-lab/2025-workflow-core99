#! /usr/bin/env python
import os
import sys
import argparse
import polars as pl


def main():
    p = argparse.ArgumentParser()
    p.add_argument('gather_csvs', nargs='+')
    p.add_argument('--names', required=True,
                   help='get smallest N for these names')
    p.add_argument('-n', '--num-to-print', type=int, default=5)
    p.add_argument('--match-name-col', default='name')
    p.add_argument('-o', '--output', required=True, help='output CSV')
    p.add_argument('--save-accs-only')
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
        df = df.select(["match_name", "f_match_orig", "query_name"])
        dflist.append(df)
        total_rows += len(df)

    all_df = pl.concat(dflist, how='vertical')

    print(f"searching {len(all_df)} rows / {all_df['match_name'].n_unique()} species / {all_df['query_name'].n_unique()} metagenomes for {args.num_to_print} metagenomes with lowest detection")

    names = set()
    with open(args.names, 'rt') as fp:
        for name in fp:
            names.add(name.strip())
    print(f'read {len(names)} intersect names from {args.names}')

    final_df_list = [] 
    for name in names:
        # pick out just the names we are interested in
        filtered_df = all_df.filter(pl.col("match_name") == name)
        print(f"found {len(filtered_df)} accs for {name} total.")
        assert len(filtered_df)

        filtered_df = filtered_df.sort("f_match_orig").head(args.num_to_print)

        final_df_list.append(filtered_df)

    print(f"writing CSV to output '{args.output}'")
    final_df = pl.concat(final_df_list, how='vertical')
    final_df.write_csv(args.output)

    if args.save_accs_only:
        print(f"writing accessions to output '{args.save_accs_only}'")
        query_names = set(final_df["query_name"].to_list())
        with open(args.save_accs_only, "wt") as fp:
            fp.write("\n".join(query_names))
    


if __name__ == '__main__':
    sys.exit(main())
