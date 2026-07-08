#! /usr/bin/env python
import sys
import argparse
import os


import sourmash
from sourmash.save_load import SaveSignaturesToLocation

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--query-sigs', nargs='+', required=True)
    p.add_argument('--metagenomes', required=True,
                   help='list of paths to metagenome signatures')
    p.add_argument('-k', '--ksize', default=21, type=int)
    p.add_argument('--scaled', default=1000, type=float)
    p.add_argument('-o', '--output-dir', required=True)
    args = p.parse_args()

    try:
        os.mkdir(args.output_dir)
    except FileExistsError:
        pass

    scaled = int(args.scaled)
    ksize = args.ksize

    query_mh_list = []
    for q in args.query_sigs:
        db = sourmash.load_file_as_index(q)
        db = db.select(ksize=ksize)
        for ss in db.signatures():
            mh = ss.minhash
            name = ss.name
            query_mh_list.append((name,
                                  mh.downsample(scaled=scaled).flatten()))

    print(f'loaded {len(query_mh_list)} query signatures')

    lines = [ x.strip() for x in open(args.metagenomes) ]
    for metagenome_file in lines:
        sigs = sourmash.load_file_as_signatures(metagenome_file,
                                                ksize=ksize)
        for ss in sigs:
            metag_name = ss.name
            metag_mh = ss.minhash.downsample(scaled=scaled)
            metag_flat_mh = metag_mh.flatten()

            for (query_name, query_mh) in query_mh_list:
                isect_mh = metag_flat_mh.intersection(query_mh)
                isect_mh = isect_mh.inflate(metag_mh)
                print(metag_name, query_name, len(isect_mh))

                isect_name = f'{query_name}.x.{metag_name} isect'
                filename = f'{args.output_dir}/{query_name}.x.{metag_name}.sig.zip'
                isect_ss = sourmash.SourmashSignature(isect_mh, isect_name)

                with SaveSignaturesToLocation(filename) as save_ss:
                    save_ss.add(isect_ss)



if __name__ == '__main__':
    sys.exit(main())
