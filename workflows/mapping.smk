# mapping stuff - figure 1 and associated data.
# NOTES:
#    mapping requires about 30 GB of RAM for the largest pangenome (phaecicola)

wildcard_constraints:
    species="[^/]+"

NAMESPLUS = [ x.strip() for x in open('inputs.cds/names-plus.list') ]
print(f'loaded {len(NAMESPLUS)} coreplus species names')

rule do_map_cds:
    input:
        expand('outputs.mapping/bams.cds.rand/{m}.x.{s}.bam', m=RAND_METAG, s=NAMESPLUS),
        expand('outputs.mapping/bams.cds.highcov/{m}.x.{s}.bam', m=HIGHCOV_METAG, s=NAMESPLUS),
        expand('outputs.mapping/bams.cds.rand/{m}.x.{s}.fastq.gz', m=RAND_METAG, s=NAMESPLUS),
        expand('outputs.mapping/bams.cds.rand/{m}.x.{s}.sig.zip', m=RAND_METAG, s=NAMESPLUS),
        expand('outputs.mapping/bams.cds.rand/{m}.x.{s}.readstats.txt', m=RAND_METAG, s=NAMESPLUS),
        expand('outputs.mapping/bams.cds.rand/{s}.readstats.csv', m=RAND_METAG, s=NAMESPLUS),

rule do_pangenomes:
    input:
        expand('outputs.mapping/genomes/{s}.pangenome.fa.gz', s=NAMES),
        expand('outputs.mapping/genomes/{s}.pangenome.sig.zip', s=NAMES),

rule make_pangenomes:
    input:
        ncbi='outputs.mapping/genomes/{species}.ncbi.d/',
        mags='outputs.mapping/genomes/{species}.ath.d/',
    output:
        'outputs.mapping/genomes/{species}.pangenome.fa.gz',
    shell: """
        gunzip -c {input.ncbi:q}/* | gzip -c > {output:q}
        (cat {input.mags:q}/* | gzip -c >> {output:q}) || true
    """

rule make_pangenome_sigs:
    input:
        'outputs.mapping/genomes/{species}.pangenome.fa.gz',
    output:
        'outputs.mapping/genomes/{species}.pangenome.sig.zip',
    conda: "env-mapping.yml"
    shell: """
        sourmash scripts singlesketch -p dna,k=21,k=31,k=51,scaled=1000 \
            {input:q} -o {output:q} --name {wildcards.species:q}
    """

rule map_index_rand_cds:
    input:
        g="outputs.cds/cds/{s}.cds.fa.gz",
        metag=GRIST_RAND100 + "trim/{m}.trim.fq.gz",
    output:
        bam='outputs.mapping/bams.cds.rand/{m}.x.{s}.bam',
        bai='outputs.mapping/bams.cds.rand/{m}.x.{s}.bam.bai',
    threads: 8
    conda: "env-mapping.yml"
    shell: """
        minimap2 -ax sr -t {threads} {input.g:q} {input.metag:q} | samtools view -b -F 4 - | samtools sort - > {output.bam:q}
        samtools index {output.bam:q}
    """

rule map_index_highcov_cds:
    input:
        g="outputs.cds/cds/{s}.cds.fa.gz",
        metag=GRIST_HIGHCOV + "trim/{m}.trim.fq.gz",
    output:
        bam='outputs.mapping/bams.cds.highcov/{m}.x.{s}.bam',
        bai='outputs.mapping/bams.cds.highcov/{m}.x.{s}.bam.bai',
    threads: 8
    conda: "env-mapping.yml"
    shell: """
        minimap2 -ax sr -t {threads} {input.g:q} {input.metag:q} | samtools view -b -F 4 - | samtools sort - > {output.bam:q}
        samtools index {output.bam:q}
    """

rule make_mapped_read_sigs:
    input:
        '{dir}/{name}.fastq.gz',
    output:
        '{dir}/{name}.sig.zip',
    conda: "env-mapping.yml"
    shell: """
        sourmash scripts singlesketch -p dna,k=21,k=31,k=51,scaled=1000,abund \
            {input:q} -o {output:q} --name {wildcards.name:q}
    """

rule map_readstats:
    input:
        "{filename}.bam",
    output:
        "{filename}.readstats.txt",
    conda: "env-mapping.yml"
    shell: """
       samtools view -c -F 260 {input:q} > {output:q}
    """

rule rand_metags_mf_csv:
    input:
        expand(GRIST_RAND100 + "sigs/{m}.trim.sig.zip", m=RAND_METAG)
    output:
        "outputs.mapping/rand-metags.mf.csv",
    conda: "env-mapping.yml"
    shell: """
        sourmash sig collect --abspath -F csv -o {output} {input}
    """

rule manysearch_rand:
    input:
        pg="outputs.mapping/genomes/{s}.pangenome.sig.zip",
        metags="outputs.mapping/rand-metags.mf.csv",
    output:
        "outputs.mapping/containment/{s}.x.rand-metags.manysearch.csv"
    threads: 32
    conda: "env-mapping.yml"
    shell: """
       sourmash scripts manysearch -k 31 -s 1000 -t 0 -c {threads} \
           {input.pg:q} {input.metags:q} -o {output:q}
    """

rule combine_readstats_rand:
    input:
        expand('outputs.mapping/{{dir}}/{m}.x.{{s}}.readstats.txt', m=RAND_METAG),
    output:
        'outputs.mapping/{dir}/{s}.readstats.csv'
    run:
        readcounts = []
        for txtname in input:
            metag = os.path.basename(txtname).split('.')[0]
            with open(txtname) as fp:
                count = int(fp.readline().strip())
                readcounts.append(dict(metag=metag, count=count))

        df = pd.DataFrame(readcounts)
        df.to_csv(str(output[0]))

rule extract_reads:
    input:
        "{ms}.bam"
    output:
        "{ms}.fastq.gz"
    conda: "env-mapping.yml"
    shell: """
        samtools fastq {input:q} | gzip -9c > {output:q}
    """

rule cat_species_reads:
    input:
        expand('outputs.mapping/bams.rand/{m}.x.{{s}}.fastq.gz', m=RAND_METAG),
    output:
        'outputs.mapping/species-reads.rand/{s}.fastq.gz'
    shell: """
        gunzip -c {input:q} | gzip -9c > {output:q}
    """
    
        
rule seqtk_sample_reads:
    input:
        "outputs.mapping/bams.lowest/{ms}.fastq.gz"
    output:
        "outputs.mapping/bams.lowest/{ms}.fastq.sample{r}"
    conda: "env-mapping.yml"
    shell: """
        seqtk sample -s {wildcards.r} {input:q} 20000 > {output:q} 
    """

rule get_individual_sigs_list:
    input:
        ath='outputs.mapping/lists/{s}.ath-tax.csv',
        gtdb='outputs.mapping/lists/{s}.gtdb-tax.csv',
    output:
        combo='outputs.mapping/sigs/{s}-genomes.picklist.csv',
    shell: """
        head -1 {input[0]:q} > {output:q}
        tail -n +2 {input.ath:q} >> {output:q}
        tail -n +2 {input.gtdb:q} >> {output:q}
    """
        

rule get_individual_sigs:
    input:
        combo='outputs.mapping/sigs/{s}-genomes.picklist.csv',
        dbs=['/home/ctbrown/scratch3/sourmash-midgie-raker/outputs.ath/rename/bin-sketches.renamed.sig.zip', '/group/ctbrowngrp5/sourmash-db.new/gtdb-rs226/gtdb-rs226-k31.dna.zip'],
    output:
        sigs='outputs.mapping/sigs/{s}-genomes.sig.zip',
    conda: "env-mapping.yml"
    shell: """
        sourmash sig cat --picklist {input.combo:q}:ident:identprefix -o {output.sigs:q} {input.dbs:q} -k 31
    """
        

rule get_individual_sigs_k51:
    input:
        combo='outputs.mapping/sigs/{s}-genomes.picklist.csv',
        dbs=['/home/ctbrown/scratch3/sourmash-midgie-raker/outputs.ath/rename/bin-sketches.renamed.sig.zip', '/group/ctbrowngrp5/sourmash-db.new/gtdb-rs226/gtdb-rs226-k51.dna.zip'],
    output:
        sigs='outputs.mapping/sigs/{s}-genomes.k51.sig.zip',
    conda: "env-mapping.yml"
    shell: """
        sourmash sig cat --picklist {input.combo:q}:ident:identprefix -o {output.sigs:q} {input.dbs:q} -k 51
    """
        

rule manysearch_individual_lowest:
    input:
        sigs='outputs.mapping/sigs/{s}-genomes.sig.zip',
        metags="outputs.mapping/lowest-metags.mf.csv",
    output:
        "outputs.mapping/manysearch/individual.{s}.x.lowest-metags.csv",
    threads: 8
    conda: "env-mapping.yml"
    shell: """
        sourmash scripts manysearch -k 31 -s 1000 -t 0 -c {threads} \
            {input.sigs:q} {input.metags:q} -o {output:q}
    """

rule manysearch_individual_rand:
    input:
        sigs='outputs.mapping/sigs/{s}-genomes.sig.zip',
        metags="outputs.mapping/rand-metags.mf.csv",
    output:
        "outputs.mapping/manysearch/individual.{s}.x.rand-metags.csv",
    threads: 8
    conda: "env-mapping.yml"
    shell: """
        sourmash scripts manysearch -k 31 -s 1000 -t 0 -c {threads} \
            {input.sigs:q} {input.metags:q} -o {output:q}
    """

rule manysearch_pg:
    input:
        sigs='outputs.mapping/sigs/{s}-genomes.sig.zip',
        genomes='outputs.mapping/genomes/{s}.pangenome.sig.zip',
    output:
        "outputs.mapping/manysearch/individual.{s}.x.pg.csv",
    threads: 8
    conda: "env-mapping.yml"
    shell: """
        sourmash scripts manysearch -k 31 -s 1000 -t 0 -c {threads} \
            {input.sigs:q} {input.genomes:q} -o {output:q}
    """


rule singlem:
    input:
        '{dir}/{s}.fastq.gz'
    output:
        '{dir}/{s}.profile.tsv'
    threads: 8
    conda: "env-singlem.yml"
    shell: """
         singlem pipe --threads {threads} -1 {input:q} -p {output:q} 
    """

rule taxburst:
    input:
        '{dir}/{s}.profile.tsv'
    output:
        '{dir}/{s}.profile.json'
    threads: 1
    conda: "env-singlem.yml"
    shell: """
         taxburst -F SingleM {input:q} --save-json {output:q}
    """
    
rule gather_k51:
    input:
        q='{dir}/{name}.sig.zip',
        db=[GTDB_K51, EUK_K51],
    output:
        '{dir}/{name}.gather.csv',
    threads: 1
    conda: "env-sourmash.yml"
    shell: """
        sourmash gather -k 51 --scaled 10_000 --threshold-bp=0 \
            {input.q:q} {input.db:q} -o {output:q} 
    """

rule gather_tax:
    input:
        csv='{dir}/{name}.gather.csv',
        db=[GTDB_TAX, EUK_TAX]
    output:
        '{dir}/{name}.gather.with-lineages.csv',
    params:
        outdir='{dir}/'
    threads: 1
    conda: "env-sourmash.yml"
    shell: """
        sourmash tax annotate -t {input.db:q} -g {input.csv:q} \
           -o {params.outdir}
    """
