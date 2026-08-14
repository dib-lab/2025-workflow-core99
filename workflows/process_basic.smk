# Basic outputs for initial detection of core, and extraction of
# accessions for various purposes.

rule do_core:
    input:
#        'outputs.core/lowest-detection.metags.txt',
        'outputs.core/core-names.list',
#        'outputs.core/core-gather.csv'

rule lowest_detection:
    input:
        names='inputs.mapping/names.with-accs.list',
    output:
        txt='outputs.core/lowest-detection.metags.txt',
        csv='outputs.core/lowest-detection.metags.csv',
    shell: """
        ./scripts/get-lowest-detection.py inputs.pangenome-gather/*.parquet \
            --names {input.names} -o {output.csv} --save-accs {output.txt}
    """

# read manysearch results => output high frequency species
rule get_core:
    input:
        dir="outputs.core2",
        metadata="inputs.metadata/SRA_meta.3216.csv"
    output:
        csv='outputs.core/core-names.csv',
        acc_txt='outputs.core/core-names.acc.list',
        txt='outputs.core/core-names.list',
    shell: """
        scripts/get-core-names-2.py -o {output.csv} \
           --save-only-names {output.txt} --save-acc-names {output.acc_txt} \
           --expect-num=3216 \
           {input.dir}/*.parquet -m {input.metadata}
    """
