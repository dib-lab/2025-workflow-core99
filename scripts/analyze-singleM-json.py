#! /usr/bin/env python
import sys
import argparse
import json
from collections import defaultdict

import taxburst
from taxburst import checks, tree_utils


def main():
    p = argparse.ArgumentParser()
    p.add_argument('tree',
                   help='tree file in taxburst JSON format')
    
    args = p.parse_args()

    treelist = []               # may not need.

    with open(args.tree, 'rb') as fp:
        tree = json.load(fp)
    nodes = checks.collect_all_nodes(tree)
    print(f"loaded {len(nodes)} nodes from '{args.tree}'")

    for rank in 'species', 'genus', 'family':
        print(f"=== {rank} ===")
        tracking = defaultdict(int)
        for node in nodes:
            name = node["name"]
            assert name         # don't allow names to be empty
            assert node["rank"]

            if node["rank"] == rank:
                tracking[node["name"]] += node["count"]

        total_count = sum(tracking.values())
        tracking = list(tracking.items())
        tracking.sort(key=lambda x: -x[1])
        for name, count in tracking:
            frac = count/total_count
            if frac > .01:
                print(f"{count/total_count * 100:>3.1f}% {name}")


if __name__ == '__main__':
    sys.exit(main())
