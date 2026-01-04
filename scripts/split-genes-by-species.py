#! /usr/bin/env python
import sys
import argparse
import csv
from collections import defaultdict
import os

import screed


def main():
    p = argparse.ArgumentParser()
    p.add_argument('genes_by_species_csv',
                   help='manually curated genes')
    p.add_argument('fasta_files', nargs='+')
    p.add_argument('-o', '--outdir', help='output directory')
    args = p.parse_args()

    args.outdir = args.outdir.rstrip('/')
    assert args.outdir

    with open(args.genes_by_species_csv, newline='') as fp:
        r = csv.DictReader(fp)
        rows = list(r)

    gene_species = {}
    genes_uniq = set()
    for row in rows:
        is_good = int(row["good"])
        if is_good:
            gene = row["gene_name"]
            species = row["species"]
            gene_species[gene] = species

            assert gene not in genes_uniq, f"{gene} is present at least twice..."
            genes_uniq.add(gene)

    try:
        os.mkdir(args.outdir)
    except:
        pass

    not_found = set(gene_species)
    species_seqs = defaultdict(list)
    for filename in args.fasta_files:
        for record in screed.open(filename):
            name = record.name.split(' ')[0]
            if name in gene_species:
                not_found.remove(name)
                species = gene_species[name]
                species_seqs[species].append(record)

    assert not not_found, not_found
    for species, records in species_seqs.items():
        outfile = f"{args.outdir}/{species}.genes.fa"
        with open(outfile, "wt") as fp:
            for r in records:
                fp.write(f">{r.name}\n{r.sequence}\n")
        print(f"wrote {len(records)} gene sequences to '{outfile}'")


if __name__ == '__main__':
    sys.exit(main())
