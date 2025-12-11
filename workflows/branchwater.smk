# branchwater stuff - figure 3 and associated.

rule do_branchwater:
    input:
        expand("outputs.branchwater/manysearch.{species}.csv", species=NAMES),
        expand("outputs.branchwater/manysearch.cds.{species}.csv", species=SUB_SPECIES),

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
