rule do_core:
    input:
        'outputs.core/lowest-detection.metags.txt',
        #'outputs.core/names.list',

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

