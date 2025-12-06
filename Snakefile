import csv
import pandas as pd

NAMES = [ x.strip() for x in open('inputs.mapping/names.list') ]
print(f'loaded {len(NAMES)} core species names')

LOWEST_METAG = [ x.strip() for x in open('inputs.mapping/lowest-detection.metags.txt') ]
print(f'loaded {len(LOWEST_METAG)} metagenome names for lowest-detection mapping')
print(LOWEST_METAG[:3])

RAND_METAG = [ x.strip() for x in open('inputs.mapping/rand_subset.3216.100.txt') ]
print(f'loaded {len(RAND_METAG)} metagenome names for rand mapping')
print(RAND_METAG[:3])

include: "workflows/mapping.smk"
