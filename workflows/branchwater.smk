# branchwater stuff - figure 3 and associated.
import os.path

rule do_branchwater:
    input:
        expand("outputs.branchwater/manysearch.{species}.parquet", species=NAMES),
        expand("outputs.branchwater/manysearch.cds.{species}.parquet", species=SUB_SPECIES),
        expand("outputs.branchwater/manysearch.cds-100k.{species}.parquet", species=SUB_SPECIES),
        expand("outputs.branchwater/manysearch.cds-500k.{species}.parquet", species=SUB_SPECIES),
        expand("outputs.branchwater/manysearch.cds-10k.{species}.parquet", species=SUB_SPECIES),

rule csv_to_parquet:
    input:
        '{filename}.csv',
    output:
        '{filename}.parquet',
    params:
        dirname=lambda w: os.path.dirname(w.filename)
    shell: """
        scripts/csv-to-parquet.py {input:q} -o {params.dirname:q}
    """

rule search_species:
    input:
        q="inputs.branchwater/{species}.pangenome.sig.zip",
        db=BW_ROCKSDB,
    output:
        "outputs.branchwater/manysearch.{species}.csv",
    conda: "env-sourmash.yml"
    threads: 1
    shell: """
        sourmash scripts manysearch -c {threads} -t 0 -k 21 -s 1000 \
            {input.q:q} {input.db:q} -o {output:q}
    """

rule search_species_cds:
    input:
        q="inputs.branchwater/{species}.cds.sig.zip",
        db=BW_ROCKSDB,
    output:
        "outputs.branchwater/manysearch.cds.{species}.csv",
    conda: "env-sourmash.yml"
    threads: 1
    shell: """
        sourmash scripts manysearch -c {threads} -t 0 -k 21 -s 1000 \
            {input.q:q} {input.db:q} -o {output:q}
    """

rule search_species_cds_100k:
    input:
        q="inputs.branchwater/{species}.cds.sig.zip",
        db=BW_ROCKSDB,
    output:
        "outputs.branchwater/manysearch.cds-100k.{species}.csv",
    conda: "env-sourmash.yml"
    threads: 1
    shell: """
        sourmash scripts manysearch -c {threads} -t 0 -k 21 -s 100_000 \
            {input.q:q} {input.db:q} -o {output:q}
    """

rule search_species_cds_10k:
    input:
        q="inputs.branchwater/{species}.cds.sig.zip",
        db=BW_ROCKSDB,
    output:
        "outputs.branchwater/manysearch.cds-10k.{species}.csv",
    conda: "env-sourmash.yml"
    threads: 1
    shell: """
        sourmash scripts manysearch -c {threads} -t 0 -k 21 -s 10_000 \
            {input.q:q} {input.db:q} -o {output:q}
    """

rule search_species_cds_500k:
    input:
        q="inputs.branchwater/{species}.cds.sig.zip",
        db=BW_ROCKSDB,
    output:
        "outputs.branchwater/manysearch.cds-500k.{species}.csv",
    conda: "env-sourmash.yml"
    threads: 1
    shell: """
        sourmash scripts manysearch -c {threads} -t 0 -k 21 -s 500_000 \
            {input.q:q} {input.db:q} -o {output:q}
    """
