#! /usr/bin/env python
import sys
import argparse
import screed
import os


def main():
    p = argparse.ArgumentParser()
    p.add_argument('fa')
    p.add_argument('-n', '--num-sequences', default=50)
    args = p.parse_args()

    chunk = 0
    filename = os.path.basename(args.fa) + f".chunk{chunk}"
    print(filename)
    outfp = open(filename, "wt")

    n = 0
    for record in screed.open(args.fa):
        outfp.write(f">{record.name}\n{record.sequence}\n")

        n += 1
        if n == args.num_sequences:
            outfp.close()
            chunk += 1
            filename = os.path.basename(args.fa) + f".chunk{chunk}"
            print(filename)
            outfp = open(filename, "wt")
            n = 0

if __name__ == '__main__':
    sys.exit(main())
