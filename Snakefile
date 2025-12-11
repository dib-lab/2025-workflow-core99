import csv
import pandas as pd

## Mapping inputs

GRIST_LOWEST='/home/ctbrown/scratch3/2025-grist-annie/outputs.lowest-detection/'
GRIST_RAND100='/home/ctbrown/scratch3/2025-grist-annie/outputs.rand100/'

NAMES = [ x.strip() for x in open('inputs.mapping/names.list') ]
print(f'loaded {len(NAMES)} core species names')
print(NAMES[0])

LOWEST_METAG = [ x.strip() for x in open('inputs.mapping/lowest-detection.metags.txt') ]
print(f'loaded {len(LOWEST_METAG)} metagenome names for lowest-detection mapping')
print(LOWEST_METAG[:3])

RAND_METAG = [ x.strip() for x in open('inputs.mapping/rand_subset.3216.100.txt') ]
print(f'loaded {len(RAND_METAG)} metagenome names for rand mapping')
print(RAND_METAG[:3])

## Branchwater inputs

BW_ROCKSDB='/group/ctbrowngrp5/sra-metagenomes/20241128-k21-s1000'
#BW_METADATA='/group/ctbrowngrp5/sra-metagenomes/20241128-metadata.parquet'

## Databases

GTDB_K51 = '/group/ctbrowngrp5/sourmash-db.new/gtdb-rs226/gtdb-rs226-k51.dna.rocksdb'
GTDB_TAX = '/group/ctbrowngrp5/sourmash-db.new/gtdb-rs226/gtdb-rs226.lineages.sqldb'
EUK_K51='/group/ctbrowngrp5/sourmash-db.new/ncbi-euks-2025.01/ncbi-euks-all-2025.01-k51.dna.rocksdb'
EUK_TAX = '/group/ctbrowngrp5/sourmash-db.new/ncbi-euks-2025.01/ncbi-eukaryotes.2025.01.lineages.sqldb'

SUB_SPECIES=('s__Cryptobacteroides sp900546925',
             's__Phascolarctobacterium_A succinatutens',
             's__Mogibacterium_A kristiansenii',
             's__Prevotella sp002251295',
             's__Prevotella sp000434975',
             's__Holdemanella porci')

###

include: "workflows/process_basic.smk"

include: "workflows/mapping.smk"

include: "workflows/branchwater.smk"
