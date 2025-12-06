# branchwater stuff - figure 3 and associated.

rule do_branchwater:
    input:
        expand("outputs.branchwater/manysearch.{species}.csv", species=NAMES)

rule search_species:
    input:
        q="inputs.branchwater/{species}.pangenome.sig.zip",
        db=BW_ROCKSDB,
    output:
        "outputs.branchwater/manysearch.{species}.csv",
    conda: "env-sourmash.yml"
    threads: 1
    shell: """
        sourmash scripts manysearch -c {threads} -t 0.1 -k 21 -s 1000 \
            {input.q:q} {input.db:q} -o {output:q}
    """
