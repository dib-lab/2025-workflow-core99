import csv
import pandas as pd

## Mapping inputs

GRIST_LOWEST='/home/ctbrown/scratch3/2025-grist-annie/outputs.lowest-detection/'
GRIST_RAND100='/home/ctbrown/scratch3/2025-grist-annie/outputs.rand100/'

NAMES = [ x.strip() for x in open('inputs.mapping/names.list') ]
print(f'loaded {len(NAMES)} core species names')

LOWEST_METAG = [ x.strip() for x in open('inputs.mapping/lowest-detection.metags.txt') ]
print(f'loaded {len(LOWEST_METAG)} metagenome names for lowest-detection mapping')
print(LOWEST_METAG[:3])

RAND_METAG = [ x.strip() for x in open('inputs.mapping/rand_subset.3216.100.txt') ]
print(f'loaded {len(RAND_METAG)} metagenome names for rand mapping')
print(RAND_METAG[:3])

## Branchwater inputs

BW_ROCKSDB='/group/ctbrowngrp5/sra-metagenomes/20241128-k21-s1000'
#BW_METADATA='/group/ctbrowngrp5/sra-metagenomes/20241128-metadata.parquet'

###

include: "workflows/mapping.smk"

include: "workflows/branchwater.smk"
