import pandas as pd

DBGTDB = '/home/ctbrown/scratch3/2025-workflow-core99/outputs.cds/cds3/clean-gtdb.species.singleton.k21.rocksdb'
DBMAG= '/home/ctbrown/scratch3/2025-workflow-core99/outputs.cds/cds3/clean-gtdb+bins.species.singleton.k21.rocksdb'

KSIZE =  21

# set list of samples
WORT_METAG = pd.read_csv("/group/ctbrowngrp2/amhorst/2025-pigparadigm/resources/metag-wort-hq.3217.txt", usecols=[0], header=None).squeeze().tolist()
#WORT_METAG=WORT_METAG[:100]

rule do_manysearch:
    input:
        expand('outputs.core2/{m}.x.mags.manysearch.parquet', m=WORT_METAG)

rule manysearch_wc:
    input:
        sig='/group/ctbrowngrp2/amhorst/2025-pigparadigm/resources/wort-pig/{metag}.sig',
        db=DBMAG,
    output:
        'outputs.core2/{metag}.x.mags.manysearch.csv'
    threads: 1
    shell: """
        sourmash scripts manysearch {input.sig} {input.db} --threshold=0 \
           -c {threads} -o {output} -k 21
    """
