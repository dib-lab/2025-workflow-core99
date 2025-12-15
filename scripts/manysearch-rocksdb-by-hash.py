#! /usr/bin/env python
# conda env: sourmash_results
import sys
import os
import argparse
import traceback
import polars as pl
from collections import defaultdict

import sourmash


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--query-sigs', nargs='+')
    p.add_argument('--against-db')
    p.add_argument('-k', '--ksize', default=21, type=int)
    p.add_argument('-s', '--scaled', default=1000, type=int)
    p.add_argument('-o', '--output-directory', required=True)
    args = p.parse_args()

    print(f"opening '{args.against_db}'")
    against_db = sourmash.load_file_as_index(args.against_db)
    against_db = against_db.select(ksize=args.ksize, scaled=args.scaled)
    print(f"len: {len(against_db)}")

    try:
        os.mkdir(args.output_directory)
    except:
        print(f"WARNING: cannot create output directory '{args.output_directory}'. Maybe it already exists? I'm sure it's fine, just fine.")

    hashval_to_queries = defaultdict(list)
    for query_sigfile in args.query_sigs:
        print(f"opening query sigs from: '{query_sigfile}'")
        db = sourmash.load_file_as_index(query_sigfile)
        db = db.select(ksize=args.ksize, scaled=args.scaled)
        print(f"len: {len(db)} sigs")

        did = 0
        for n, ss in enumerate(db.signatures()):
            matches = []
            species_name = ss.name
            outname = os.path.join(args.output_directory,
                                   f"{species_name}.parquet")
            if os.path.exists(outname):
                print(f"{outname} exists; skipping")
                continue

            mh = ss.minhash.downsample(scaled=args.scaled)
            print(species_name, len(mh))

            for hash_i, hashval in enumerate(mh.hashes):
                qq = mh.copy_and_clear()
                qq.add_hash(hashval)
                sq = sourmash.SourmashSignature(qq)

                try:
                    counter = against_db.counter_gather(sq)
                except:
                    print(f"ERROR on species name: {species_name}")
                    traceback.print_exc()
                    continue

                for (size, name) in counter.matches(threshold_hashes=0):
                    dd = dict(hashval=hashval, acc=name, species=species_name)
                    matches.append(dd)
                print(hash_i, len(mh), hashval, len(matches))

            df = pl.DataFrame(matches)
            df.write_parquet(outname)
            print(n, len(db), species_name, len(df))


if __name__ == '__main__':
    sys.exit(main())
