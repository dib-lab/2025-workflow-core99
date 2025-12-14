#! /usr/bin/env python
import sys
import argparse
import sourmash
from collections import Counter
from sourmash.sourmash_args import SaveSignaturesToLocation
import os


def main():
    p = argparse.ArgumentParser()
    p.add_argument('pangenome_db', nargs='+')
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('-o', '--output-directory', required=True)
    args = p.parse_args()

    # load and count.
    seen = set()
    removed = 0
    for filename in args.pangenome_db:
        db = sourmash.load_file_as_index(filename)
        db = db.select(ksize=args.ksize)
        outfile = os.path.join(args.output_directory, os.path.basename(filename))
        with SaveSignaturesToLocation(outfile) as save_sig:
            for n, ss in enumerate(db.signatures()):
                if n % 1000 == 0:
                    print('...', n)

                mh = ss.minhash
                hashes = set(mh.hashes)
                hashes -= seen
                removed += len(mh) - len(hashes)

                mh = mh.copy_and_clear()
                mh.add_many(hashes)

                ss = sourmash.SourmashSignature(mh, name=ss.name)
                save_sig.add(ss)
                
                seen.update(hashes)

    print(removed, len(seen))


if __name__ == '__main__':
    sys.exit(main())
