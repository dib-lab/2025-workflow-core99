#! /usr/bin/env python
import sys
import argparse
import os
import json
from collections import defaultdict
import pprint

import polars as pl

import taxburst
from taxburst import checks, tree_utils

ranks = ('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('tree_jsons', nargs='+',
                   help='tree files in taxburst JSON format')
    p.add_argument('-t', '--taxonomy', required=True,
                   help='GTDB taxonomy CSV matching species input')
    p.add_argument('-v', '--verbose', action='store_true')
    args = p.parse_args()

    tax_df = pl.read_csv(args.taxonomy)

    n_species_dominant = defaultdict(int)
    n_species_nondom = defaultdict(int)
    n_species_mis = defaultdict(int)
    n_species_not_exist = defaultdict(int)
    n_total_files = defaultdict(int)
    n_genus_match = defaultdict(int)

    for jsonfile in args.tree_jsons:
        this_species = os.path.basename(jsonfile).split('.')[2]
        n_total_files[this_species] += 1

        taxrow = (tax_df.row(named=True,
                             by_predicate=(pl.col("species") == this_species)))
        
        with open(jsonfile, 'rb') as fp:
            tree = json.load(fp)
        nodes = checks.collect_all_nodes(tree)
        if args.verbose:
            print('===')

        rank_vals = {}
        for rank in ranks:
            tracking = defaultdict(int)
            for node in nodes:
                name = node["name"]
                assert name         # don't allow names to be empty
                assert node["rank"]

                if node["rank"] == rank:
                    tracking[node["name"]] += node["count"]

            if not tracking:
                rank_vals[rank] = (None, None, None)
                continue

            total_count = sum(tracking.values())
            tracking = list(sorted(tracking.items(), key=lambda x: -x[1]))

            name, frac = tracking[0]
            dominance = frac / total_count

            rank_vals[rank] = (name, dominance, name == taxrow[rank])

        #pprint.pprint(rank_vals)

        # first, check if mismatch at or above family.
        DONE = False
        for check_rank in ('superkingdom', 'phylum', 'class', 'order', 'family'):
            (name, dom, is_match) = rank_vals[check_rank]
            if name is None:
                #print(f'NO MATCH at {check_rank} for {this_species}')
                n_species_not_exist[this_species] += 1
                DONE = True
                break
            if not is_match:
                print(f'MISMATCH: {taxrow[check_rank]} != {name} for {this_species}')
                n_species_mis[this_species] += 1
                DONE = True
                break

        if DONE:
            continue

        # second, check match at species:
        check_rank = "species"
        (name, dom, is_match) = rank_vals[check_rank]
        if name is None:        # skip
            pass
        elif is_match:
            if dom >= 0.95:
                n_species_dominant[this_species] += 1
            else:
                n_species_nondom[this_species] += 1
            DONE = True
        else:
            # is genus a match?
            (genus_name, genus_dom, genus_is_match) = rank_vals["genus"]
            if genus_is_match: # and genus_dom >= 0.95:
                n_genus_match[this_species] += 1
            else:
                n_species_mis[this_species] += 1
                print(f'MISMATCH_SP: {taxrow[check_rank]} != {name} for {this_species}')
                #print(rank_vals)
            DONE = True

        if DONE:
            continue

        # catchall for now
        n_species_not_exist[this_species] += 1
        

    summary = []
    for species in n_total_files:
        summary.append((n_species_dominant[species],
                        n_species_nondom[species],
                        n_genus_match[species],
                        n_species_mis[species],
                        n_species_not_exist[species],
                        n_total_files[species],
                        species))

    summary.sort(key=lambda x: x[3])

    print('===')
    print('species                                  dom  nondom  genus  mismatch  nocall')
    for dom, nondom, genmatch, mis, not_exist, total, species in summary:
        assert dom + nondom + genmatch + mis + not_exist == total
        print(f"{species:<40} {dom:<7} {nondom:<7} {genmatch:<7} {mis:<7} {not_exist:<7}")


if __name__ == '__main__':
    sys.exit(main())
