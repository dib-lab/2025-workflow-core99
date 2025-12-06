# mapping stuff - figure 1 and associated data.
# NOTES:
#    mapping requires about 30 GB of RAM for the largest pangenome (phaecicola)

rule do_mapping:
    input:
        expand('outputs.mapping/lists/{s}.ath-tax.csv', s=NAMES),
        expand('outputs.mapping/lists/{s}.gtdb-tax.csv', s=NAMES),
        expand('outputs.mapping/lists/{s}.gtdb-acc.txt', s=NAMES),
        expand('outputs.mapping/genomes/{s}.ncbi.d/', s=NAMES),
        expand('outputs.mapping/lists/{s}.ath-paths.txt', s=NAMES),
        expand('outputs.mapping/genomes/{s}.ath.d/', s=NAMES),
        expand('outputs.mapping/genomes/{s}.pangenome.fa.gz', s=NAMES),
        expand('outputs.mapping/genomes/{s}.pangenome.sig.zip', s=NAMES),
        expand("outputs.mapping/containment/{s}.x.lowest-metags.manysearch.csv", s=NAMES),
        expand("outputs.mapping/containment/{s}.x.rand-metags.manysearch.csv", s=NAMES),
        expand("outputs.mapping/sigs/{s}-genomes.sig.zip", s=NAMES),
        expand("outputs.mapping/manysearch/individual.{s}.x.rand-metags.csv", s=NAMES),
        expand("outputs.mapping/manysearch/individual.{s}.x.lowest-metags.csv", s=NAMES),
        expand("outputs.mapping/manysearch/individual.{s}.x.pg.csv", s=NAMES),
        expand('outputs.mapping/bams.rand/{m}.x.{s}.bam', m=RAND_METAG, s=NAMES),
        expand('outputs.mapping/bams.rand/{m}.x.{s}.fastq', m=RAND_METAG, s=NAMES),
        expand('outputs.mapping/bams.rand/{m}.x.{s}.sig.zip', m=RAND_METAG, s=NAMES),
        expand('outputs.mapping/bams.rand/{m}.x.{s}.readstats.txt', m=RAND_METAG, s=NAMES),
        expand('outputs.mapping/bams.rand/{s}.readstats.csv', m=RAND_METAG, s=NAMES),
        expand('outputs.mapping/bams.lowest/{m}.x.{s}.bam', m=LOWEST_METAG, s=NAMES),
        expand('outputs.mapping/bams.lowest/{m}.x.{s}.fastq', m=LOWEST_METAG, s=NAMES),
        expand('outputs.mapping/bams.lowest/{m}.x.{s}.sig.zip', m=LOWEST_METAG, s=NAMES),
        expand('outputs.mapping/bams.lowest/{m}.x.{s}.readstats.txt', m=LOWEST_METAG, s=NAMES),
        expand('outputs.mapping/bams.lowest/{s}.readstats.csv', m=LOWEST_METAG, s=NAMES),
        expand('outputs.mapping/species-reads.rand/{s}.fastq.gz', s=NAMES),
        expand('outputs.mapping/singlem.rand/{s}.profile.tsv', s=NAMES),

rule do_map_lowest:
    input:
        expand('outputs.mapping/bams.lowest/{m}.x.{s}.bam', m=LOWEST_METAG, s=NAMES),
        expand('outputs.mapping/bams.lowest/{m}.x.{s}.fastq', m=LOWEST_METAG, s=NAMES),
        expand('outputs.mapping/bams.lowest/{m}.x.{s}.sig.zip', m=LOWEST_METAG, s=NAMES),
        expand('outputs.mapping/bams.lowest/{m}.x.{s}.readstats.txt', m=LOWEST_METAG, s=NAMES),
        expand('outputs.mapping/bams.lowest/{s}.readstats.csv', m=LOWEST_METAG, s=NAMES),

rule do_map_rand:
    input:
        expand('outputs.mapping/bams.rand/{m}.x.{s}.bam', m=RAND_METAG, s=NAMES),
        expand('outputs.mapping/bams.rand/{m}.x.{s}.fastq', m=RAND_METAG, s=NAMES),
        expand('outputs.mapping/bams.rand/{m}.x.{s}.sig.zip', m=RAND_METAG, s=NAMES),
        expand('outputs.mapping/bams.rand/{m}.x.{s}.readstats.txt', m=RAND_METAG, s=NAMES),
        expand('outputs.mapping/bams.rand/{s}.readstats.csv', m=RAND_METAG, s=NAMES),

rule genome_lists:
    input:
        expand('outputs.mapping/lists/{s}.gtdb-acc.txt', s=NAMES),
        expand('outputs.mapping/lists/{s}.ath-paths.txt', s=NAMES),

rule do_pangenomes:
    input:
        expand('outputs.mapping/genomes/{s}.pangenome.fa.gz', s=NAMES),
        expand('outputs.mapping/genomes/{s}.pangenome.sig.zip', s=NAMES),

rule get_ath_mag_names:
    input:
        lineage='/home/ctbrown/scratch3/sourmash-midgie-raker/outputs.ath/rename/bin-sketches.lineages.csv'
    output:
        csv='outputs.mapping/lists/{species}.ath-tax.csv'
    conda: "env-mapping.yml"
    shell: """
        sourmash tax grep {wildcards.species:q} -t {input.lineage} -o {output.csv:q}
    """

rule get_gtdb_tax_names_1:
    input:
        lineage='/group/ctbrowngrp5/sourmash-db.new/gtdb-rs226/gtdb-rs226.lineages.sqldb'
    output:
        csv='outputs.mapping/lists/{species}.gtdb-tax.csv',
    conda: "env-mapping.yml"
    shell: """
        sourmash tax grep {wildcards.species:q} -t {input.lineage} -o {output.csv:q}
    """

rule get_gtdb_tax_names_2:
    input:
        csv='outputs.mapping/lists/{species}.gtdb-tax.csv',
        exclude_csv='/home/ctbrown/scratch3/sourmash-midgie-raker/outputs.ath/host/picklist-exclude.csv'
    output:
        txt='outputs.mapping/lists/{species}.gtdb-acc.txt',
    run:
        exclude = set()
        with open(input.exclude_csv, "rt") as fp:
            for row in csv.DictReader(fp):
                ident = row["ident"].split('.')[0]
                exclude.add(ident)

        outfp = open(output.txt, "wt")
        with open(input.csv, "rt") as fp:
            for row in csv.DictReader(fp):
                ident = row["ident"]
                if ident in exclude:
                    print(f"excluding {ident}")
                else:
                    print(ident, file=outfp)
    

rule get_ncbi_genomes:
    input:
        txt='outputs.mapping/lists/{species}.gtdb-acc.txt',
    output:
        directory('outputs.mapping/genomes/{species}.ncbi.d/')
    conda: "env-mapping.yml"
    shell: """
        get-some-ncbi-genomes --from-file {input:q} -G --output-dir {output:q}
    """

rule get_ath_paths:
    input:
        tax_csv='outputs.mapping/lists/{species}.ath-tax.csv',
        sketch_csv='/home/ctbrown/scratch3/sourmash-midgie-raker/outputs.ath/rename/manysketch-renamed.csv',
        exclude_csv='/home/ctbrown/scratch3/sourmash-midgie-raker/outputs.ath/host/picklist-exclude.csv'
    output:
        csv='outputs.mapping/lists/{species}.ath-paths.txt'
    run:
        exclude = set()
        with open(input.exclude_csv, "rt") as fp:
            for row in csv.DictReader(fp):
                ident = row["ident"].split('.')[0]
                exclude.add(ident)

        paths_by_ident = {}
        with open(input.sketch_csv, 'rt') as fp:
            r = csv.DictReader(fp)
            for row in r:
                ident = row['name'].split(' ')[0]
                path = row['genome_filename']
                paths_by_ident[ident] = path

        with open(output.csv, 'wt') as outfp:
            with open(input.tax_csv, 'rt') as fp:
                r = csv.DictReader(fp)
                for row in r:
                    ident = row['ident']
                    if ident in exclude:
                        print(f"excluding {ident}")
                    else:
                        path = paths_by_ident[ident]
                        outfp.write(f"{path}\n")

rule get_ath_genomes:
    input:
        csv='outputs.mapping/lists/{species}.ath-paths.txt'
    output:
        directory('outputs.mapping/genomes/{species}.ath.d/')
    shell: """
        mkdir -p {output:q}
        cp $(cat {input.csv:q}) {output:q} || true
    """

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

rule map_index_lowest:
    input:
        g="outputs.mapping/genomes/{s}.pangenome.fa.gz",
        metag=GRIST_LOWEST+"trim/{m}.trim.fq.gz",
    output:
        bam='outputs.mapping/bams.lowest/{m}.x.{s}.bam',
        bai='outputs.mapping/bams.lowest/{m}.x.{s}.bam.bai',
    threads: 8
    conda: "env-mapping.yml"
    shell: """
        minimap2 -ax sr -t {threads} {input.g:q} {input.metag:q} | samtools view -b -F 4 - | samtools sort - > {output.bam:q}
        samtools index {output.bam:q}
    """

rule map_index_rand:
    input:
        g="outputs.mapping/genomes/{s}.pangenome.fa.gz",
        metag=GRIST_RAND100 + "trim/{m}.trim.fq.gz",
    output:
        bam='outputs.mapping/bams.rand/{m}.x.{s}.bam',
        bai='outputs.mapping/bams.rand/{m}.x.{s}.bam.bai',
    threads: 8
    conda: "env-mapping.yml"
    shell: """
        minimap2 -ax sr -t {threads} {input.g:q} {input.metag:q} | samtools view -b -F 4 - | samtools sort - > {output.bam:q}
        samtools index {output.bam:q}
    """

rule make_mapped_read_sigs:
    input:
        'outputs.mapping/bams.{dir}/{name}.fastq',
    output:
        'outputs.mapping/bams.{dir}/{name}.sig.zip',
    conda: "env-mapping.yml"
    shell: """
        sourmash scripts singlesketch -p dna,k=21,k=31,k=51,scaled=1000 \
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

rule lowest_metags_mf_csv:
    input:
        expand(GRIST_LOWEST + "sigs/{m}.trim.sig.zip", m=LOWEST_METAG)
    output:
        "outputs.mapping/lowest-metags.mf.csv",
    conda: "env-mapping.yml"
    shell: """
        sourmash sig collect --abspath -F csv -o {output} {input}
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

rule manysearch_lowest:
    input:
        pg="outputs.mapping/genomes/{s}.pangenome.sig.zip",
        metags="outputs.mapping/lowest-metags.mf.csv",
    output:
        "outputs.mapping/containment/{s}.x.lowest-metags.manysearch.csv"
    threads: 32
    conda: "env-mapping.yml"
    shell: """
       sourmash scripts manysearch -k 31 -s 1000 -t 0 -c {threads} \
           {input.pg:q} {input.metags:q} -o {output:q}
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
        expand('outputs.mapping/bams.rand/{m}.x.{{s}}.readstats.txt', m=RAND_METAG),
    output:
        'outputs.mapping/bams.rand/{s}.readstats.csv'
    run:
        readcounts = []
        for txtname in input:
            metag = os.path.basename(txtname).split('.')[0]
            with open(txtname) as fp:
                count = int(fp.readline().strip())
                readcounts.append(dict(metag=metag, count=count))

        df = pd.DataFrame(readcounts)
        df.to_csv(str(output[0]))

rule combine_readstats_lowest:
    input:
        expand('outputs.mapping/bams.lowest/{m}.x.{{s}}.readstats.txt', m=LOWEST_METAG),
    output:
        'outputs.mapping/bams.lowest/{s}.readstats.csv'
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
        "{ms}.fastq"
    conda: "env-mapping.yml"
    shell: """
        samtools fastq {input:q} > {output:q}
    """

rule cat_species_reads:
    input:
        expand('outputs.mapping/bams.rand/{m}.x.{{s}}.fastq', m=RAND_METAG),
    output:
        'outputs.mapping/species-reads.rand/{s}.fastq.gz'
    shell: """
        cat {input:q} | gzip -9c > {output:q}
    """
    
        
rule seqtk_sample_reads:
    input:
        "outputs.mapping/bams.lowest/{ms}.fastq"
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
        'outputs.mapping/species-reads.rand/{s}.fastq.gz'
    output:
        'outputs.mapping/singlem.rand/{s}.profile.tsv'
    threads: 8
    conda: "env-singlem.yml"
    shell: """
         singlem pipe --threads {threads} -1 {input:q} -p {output:q} 
    """
