#! /usr/bin/env python
import polars as pl
import sys
import argparse

RANKS="superkingdom,phylum,class,order,family,genus,species".split(",")

def main():
    p = argparse.ArgumentParser()
    p.add_argument('names')
    p.add_argument('tax_parquet')
    p.add_argument('-o', '--output-csv', required=True)
    args = p.parse_args()

    names = set()
    with open(args.names, 'rt') as fp:
        names.update([ x.strip() for x in fp ])

    print(f"read {len(names)} species names")

    print(f"reading tax csv {args.tax_parquet}")
    df = (pl.scan_parquet(args.tax_parquet)
          .filter(pl.col("species").is_in(names))
          .select(RANKS)
          .unique()
          ).collect()

    assert len(df) == len(names)

    print(f"writing {len(df)} lineages to {args.output_csv}")
    df.write_csv(args.output_csv)


if __name__ == '__main__':
    sys.exit(main())


    
