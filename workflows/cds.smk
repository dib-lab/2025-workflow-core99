# coding sequence stuff - figure 1 and associated data.
# TODO:
#  - automate cd-hit clustering of all + cluster extraction
# NOTES:
#    mapping requires about 30 GB of RAM for the largest pangenome (phaecicola)

wildcard_constraints:
    species="[^/]+"

NAMESPLUS = [ x.strip() for x in open('inputs.mapping/names-plus.list') ]
print(f'loaded {len(NAMESPLUS)} coreplus species names')

SPECIES_WITH_GENES, =glob_wildcards('outputs.cds/cds3-genes/{s}.genes.fa')

ncbi_genomes = []
ncbi_species = []
ath_genomes = []
ath_species = []

for species in NAMESPLUS:
    g, = glob_wildcards(f'outputs.cds/genomes/{species}.ncbi.d/{{g}}.fna.gz')
    for x in g:
        ncbi_species.append(species)
        ncbi_genomes.append(x)

    g, = glob_wildcards(f'outputs.cds/genomes/{species}.ath.d/{{g}}.fasta')

    for x in g:
        ath_species.append(species)
        ath_genomes.append(x)

rule do_cds:
    input:
        #expand('outputs.cds/genomes/{s}.ncbi.d/', s=NAMES),
        #expand('outputs.cds/lists/{s}.ath-paths.txt', s=NAMES),
        #expand('outputs.cds/genomes/{s}.ath.d/', s=NAMES),
        expand('outputs.cds/cds/{species}.cds.fa.gz', species=NAMESPLUS),
        expand('outputs.cds/cds/{species}.cds.sig.zip', species=NAMESPLUS),
        expand('outputs.cds/cds/gtdb-only/{species}.cds.fa.gz', species=NAMESPLUS),
        expand('outputs.cds/cds/gtdb-only/{species}.cds.sig.zip', species=NAMESPLUS),
        expand('outputs.cds/cds/ath-only/{species}.cds.fa.gz', species=NAMESPLUS),
        expand('outputs.cds/cds/ath-only/{species}.cds.sig.zip', species=NAMESPLUS),
        'outputs.cds/cds3/clean-gtdb+bins.species.singleton.k21.sig.zip',
        'outputs.cds/cds3/clean-gtdb.species.singleton.k21.sig.zip',
        expand('outputs.cds/cds3/{s}.cds3.sig.zip', s=NAMESPLUS),
        expand('outputs.cds/cds3/{s}.cds4.sig.zip', s=NAMESPLUS),
        expand('outputs.cds/cds3/singleton_pg/{s}.sig.zip', s=NAMESPLUS),
        expand("outputs.cds/cds3/manysearch.cds{c}.3216.csv", c=[3,4,5,6]),

rule foo:
    input:
        expand('outputs.cds/cds3/outputs.branchwater.cds{c}/manysearch.{s}.parquet', s=NAMESPLUS,
               c=[3,4,5,6]),
        expand('outputs.cds/cds3/outputs.minsig/{s}.cds3.min{m}.sig.zip', s=NAMESPLUS, m=[10, 20, 50, 60, 70, 80, 90]),
        expand('outputs.cds/cds3/outputs.branchwater.cds3.min50/manysearch.{s}.parquet', s=NAMESPLUS),
        "outputs.cds/cds3/manysearch.cds3.min50.3216.csv",
        expand('outputs.cds/cds3/outputs.minsig/{s}.cds3.min50.fa', s=NAMESPLUS),
        expand('outputs.cds/cds3/outputs.minsig/{s}.cds3.min50.dedup.fa', s=NAMESPLUS),
        expand('outputs.cds/cds3/outputs.minsig/{s}.cds3.min50.dedup2.fa', s=MIN50_NAMES),
        expand('outputs.cds/cds3/{s}.cds3.min50.dedup.sig.zip', s=MIN50_NAMES),
        expand('outputs.cds/cds3/{s}.cds3.singleclust.sig.zip', s=MIN50_NAMES),
        expand('outputs.cds/cds3/outputs.branchwater.cds3.min50.dedup/manysearch.{s}.parquet', s=MIN50_NAMES),
        "outputs.cds/cds3/manysearch.cds3.min50.dedup.3216.csv",
        expand('outputs.cds/cds3/outputs.branchwater.cds3.singleclust/manysearch.{s}.parquet', s=MIN50_NAMES),
        "outputs.cds/cds3/manysearch.cds3.singleclust.3216.csv",
        expand("outputs.cds/cds3-genes/{s}.sig.zip", s=SPECIES_WITH_GENES),
        "outputs.cds/cds3-genes/manysearch.3216.csv",
        expand("outputs.cds/cds3-genes/outputs.branchwater/manysearch.{s}.parquet", s=SPECIES_WITH_GENES),

####

rule do_prokka_ncbi:
    input:
        expand('outputs.cds/genomes/{species}.ncbi.prokka/{g}.prokka.d', zip, species=ncbi_species, g=ncbi_genomes)

rule do_prokka_ath:
    input:
        expand('outputs.cds/genomes/{species}.ath.prokka/{g}.prokka.d', zip, species=ath_species, g=ath_genomes)

rule genome_lists:
    input:
        expand('outputs.cds/lists/{s}.gtdb-acc.txt', s=NAMES),
        expand('outputs.cds/lists/{s}.ath-paths.txt', s=NAMES),

rule get_ath_mag_names:
    input:
        lineage='/home/ctbrown/scratch3/sourmash-midgie-raker/outputs.ath/rename/bin-sketches.lineages.csv'
    output:
        csv='outputs.cds/lists/{species}.ath-tax.csv'
    conda: "env-mapping.yml"
    shell: """
        sourmash tax grep {wildcards.species:q} -t {input.lineage} -o {output.csv:q}
    """

rule get_gtdb_tax_names_1:
    input:
        lineage='/group/ctbrowngrp5/sourmash-db.new/gtdb-rs226/gtdb-rs226.lineages.sqldb'
    output:
        csv='outputs.cds/lists/{species}.gtdb-tax.csv',
    conda: "env-mapping.yml"
    shell: """
        sourmash tax grep {wildcards.species:q} -t {input.lineage} -o {output.csv:q}
    """

rule get_gtdb_tax_names_2:
    input:
        csv='outputs.cds/lists/{species}.gtdb-tax.csv',
        exclude_csv='/home/ctbrown/scratch3/sourmash-midgie-raker/outputs.ath/host/picklist-exclude.csv'
    output:
        txt='outputs.cds/lists/{species}.gtdb-acc.txt',
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
        txt='outputs.cds/lists/{species}.gtdb-acc.txt',
    output:
        directory('outputs.cds/genomes/{species}.ncbi.d/')
    conda: "env-mapping.yml"
    shell: """
        get-some-ncbi-genomes --from-file {input:q} -G --output-dir {output:q}
    """

rule get_ath_paths:
    input:
        tax_csv='outputs.cds/lists/{species}.ath-tax.csv',
        sketch_csv='/home/ctbrown/scratch3/sourmash-midgie-raker/outputs.ath/rename/manysketch-renamed.csv',
        exclude_csv='/home/ctbrown/scratch3/sourmash-midgie-raker/outputs.ath/host/picklist-exclude.csv'
    output:
        csv='outputs.cds/lists/{species}.ath-paths.txt'
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
        csv='outputs.cds/lists/{species}.ath-paths.txt'
    output:
        directory('outputs.cds/genomes/{species}.ath.d/')
    shell: """
        mkdir -p {output:q}
        cp $(cat {input.csv:q}) {output:q} || true
    """


rule do_prokka_ncbi_wc:
    input:
        g='outputs.cds/genomes/{species}.ncbi.d/{g}.fna.gz'
    output:
        g=temporary("outputs.cds/genomes/{species}.ncbi.prokka/{g}.fa"),
        dir=directory('outputs.cds/genomes/{species}.ncbi.prokka/{g}.prokka.d')
    threads: 8
    conda: "env-prokka.yml"
    shell: """
        gunzip -c {input.g:q} > {output.g:q}
        prokka --outdir {output.dir:q} {output.g:q} --fast --cpus {threads}
    """

rule do_prokka_ath_wc:
    input:
        g='outputs.cds/genomes/{species}.ath.d/{g}.fasta'
    output:
        dir=directory('outputs.cds/genomes/{species}.ath.prokka/{g}.prokka.d')
    threads: 8
    conda: "env-prokka.yml"
    shell: """
        prokka --outdir {output.dir:q} {input.g:q} --fast --cpus {threads}
    """

def get_ath_prokka_dirs(w):
    ath_dirs = []
    for (species, genome) in zip(ath_species, ath_genomes):
        if species == w.species:
            ath_dir = f'outputs.cds/genomes/{species}.ath.prokka/{genome}.prokka.d'
            ath_dirs.append(ath_dir)

    if len(ath_dirs) == 0:
        print(f"WARNING, no AtH directories for '{w.species}'")

    return ath_dirs

def get_ncbi_prokka_dirs(w):
    ncbi_dirs = []
    for (species, genome) in zip(ncbi_species, ncbi_genomes):
        if species == w.species:
            ncbi_dir = f'outputs.cds/genomes/{species}.ncbi.prokka/{genome}.prokka.d'
            ncbi_dirs.append(ncbi_dir)

    assert ncbi_dirs, f"no Prokka directories for {w.species}"

    return ncbi_dirs

rule prokka_cds_wc:
    input:
        ancient(get_ath_prokka_dirs),
        ancient(get_ncbi_prokka_dirs),
    output:
        'outputs.cds/cds/{species}.cds.fa.gz'
    shell: """
        find {input:q} -name *.ffn -exec cat {{}} \\; | gzip > {output:q}
    """

rule prokka_cds_wc_gtdb_only:
    input:
        ancient(get_ncbi_prokka_dirs),
    output:
        'outputs.cds/cds/gtdb-only/{species}.cds.fa.gz'
    shell: """
        find {input:q} -name *.ffn -exec cat {{}} \\; | gzip > {output:q}
    """

rule prokka_cds_wc_ath_only:
    input:
        ancient(get_ath_prokka_dirs),
    output:
        'outputs.cds/cds/ath-only/{species}.cds.fa.gz'
    shell: """
        find {input:q} -name *.ffn -exec cat {{}} \\; | gzip > {output:q}
    """

rule make_prokka_cds_sig_zip:
    input:
        '{dir}/{name}.cds.fa.gz',
    output:
        '{dir}/{name}.cds.sig.zip',
    conda: "env-mapping.yml"
    shell: """
        sourmash scripts singlesketch -p dna,k=21,k=31,k=51,scaled=1000 \
            {input:q} -o {output:q} --name {wildcards.name:q}
    """

rule singleton_gtdb_bins:
    input:
        db_k21='../sourmash-midgie-raker/outputs.ath/host/clean-gtdb+bins.species.k21.sig.zip',
    output:
        'outputs.cds/cds3/clean-gtdb+bins.species.singleton.k21.sig.zip',
    shell: """
        scripts/remove-multihash-by-sig.py {input} -o . -k 21
        mv clean-gtdb+bins.species.k21.sig.zip outputs.cds/cds3/clean-gtdb+bins.species.singleton.k21.sig.zip
    """

rule singleton_gtdb:
    input:
        db_k21='../sourmash-midgie-raker/outputs.ath/host/clean-gtdb.species.k21.sig.zip',
    output:
        'outputs.cds/cds3/clean-gtdb.species.singleton.k21.sig.zip',
    shell: """
        scripts/remove-multihash-by-sig.py {input} -o . -k 21
        mv clean-gtdb.species.k21.sig.zip outputs.cds/cds3/clean-gtdb.species.singleton.k21.sig.zip
    """

rule extract_singleton_pg:
    input:
        db_k21='outputs.cds/cds3/clean-gtdb+bins.species.singleton.k21.sig.zip',
    output:
        'outputs.cds/cds3/singleton_pg/{s}.sig.zip'
    params:
        exact=lambda w: f"{w.s} " # add space to require exact matching
    shell: """
        sourmash sig grep {params.exact:q} {input:q} -o {output:q}
    """

# cds3: singleton CDS.
rule cds3_isect:
    input:
        pg='outputs.cds/cds3/singleton_pg/{s}.sig.zip',
        cds='outputs.cds/cds/{s}.cds.sig.zip',
    output:
        'outputs.cds/cds3/{s}.cds3.sig.zip',
    shell: """
        sourmash sig intersect {input:q} -o {output:q} --set-name {wildcards.s:q} -k 21
    """

# cds4: GTDB only
rule cds4_isect:
    input:
        pg='outputs.cds/cds3/singleton_pg/{s}.sig.zip',
        cds='outputs.cds/cds/gtdb-only/{s}.cds.sig.zip',
    output:
        'outputs.cds/cds3/{s}.cds4.sig.zip',
    shell: """
        sourmash sig intersect {input:q} -o {output:q} --set-name {wildcards.s:q} -k 21
    """

# cds5: AtH-GTDB
rule cds5_isect:
    input:
        pg='outputs.cds/cds3/singleton_pg/{s}.sig.zip',
        cds_gtdb='outputs.cds/cds/gtdb-only/{s}.cds.sig.zip',
        cds_combo='outputs.cds/cds/{s}.cds.sig.zip',
    output:
        tmp=temporary("temp/{s}.cds5_isect.sig.zip"),
        cds5='outputs.cds/cds3/{s}.cds5.sig.zip',
    shell: """
        sourmash sig subtract {input.cds_combo:q} {input.cds_gtdb:q} -o \
           {output.tmp:q} -k 21
        sourmash sig intersect {input.pg:q} {output.tmp:q} -o {output.cds5:q} --set-name {wildcards.s:q} -k 21 
    """

# cds6: GTDB-AtH
rule cds6_isect:
    input:
        pg='outputs.cds/cds3/singleton_pg/{s}.sig.zip',
        cds_gtdb='outputs.cds/cds/gtdb-only/{s}.cds.sig.zip',
        cds_ath='outputs.cds/cds/ath-only/{s}.cds.sig.zip',
    output:
        tmp=temporary("temp/{s}.cds6_isect.sig.zip"),
        cds5='outputs.cds/cds3/{s}.cds6.sig.zip',
    shell: """
        sourmash sig subtract {input.cds_gtdb:q} {input.cds_ath:q} -o \
           {output.tmp:q} -k 21
        sourmash sig intersect {input.pg:q} {output.tmp:q} -o {output.cds5:q} --set-name {wildcards.s:q} -k 21 
    """

        
rule species_cds_mf:
    input:
        expand("outputs.cds/cds3/{s}.cds{{c}}.sig.zip", s=NAMES)
    output:
        "outputs.cds/cds3/species.cds{c}.mf.csv"
    shell: """
        sourmash sig collect -F csv {input:q} -o {output:q} --abspath
    """
        
rule species_min_mf:
    input:
        expand("outputs.cds/cds3/outputs.minsig/{s}.cds{{c}}.min{{m}}.sig.zip", s=NAMES)
    output:
        "outputs.cds/cds3/species.cds{c}.min{m,^[\\.]+}.mf.csv"
    shell: """
        sourmash sig collect -F csv {input:q} -o {output:q} --abspath
    """
        
rule species_min_mf_dedup:
    input:
        expand("outputs.cds/cds3/{s}.cds3.min50.dedup.sig.zip", s=MIN50_NAMES)
    output:
        "outputs.cds/cds3/species.cds3.min50.dedup.mf.csv"
    shell: """
        sourmash sig collect -F csv {input:q} -o {output:q} --abspath
    """
        
rule species_min_mf_singleclust:
    input:
        expand("outputs.cds/cds3/{s}.cds3.singleclust.sig.zip", s=MIN50_NAMES)
    output:
        "outputs.cds/cds3/species.cds3.singleclust.mf.csv"
    shell: """
        sourmash sig collect -F csv {input:q} -o {output:q} --abspath
    """
        
rule search_species_3216:
    input:
        q="outputs.cds/cds3/species.cds{c}.mf.csv",
        db="3216.manifest.csv",
    output:
        "outputs.cds/cds3/manysearch.cds{c}.3216.csv",
    conda: "env-sourmash.yml"
    threads: 32
    shell: """
        sourmash scripts manysearch -c {threads} -t 0 -k 21 -s 1000 \
            {input.q:q} {input.db:q} -o {output:q}
    """

rule search_species_3216_min50_dedup:
    input:
        q="outputs.cds/cds3/species.cds3.min50.dedup.mf.csv",
        db="3216.manifest.csv",
    output:
        "outputs.cds/cds3/manysearch.cds3.min50.dedup.3216.csv",
    conda: "env-sourmash.yml"
    threads: 32
    shell: """
        sourmash scripts manysearch -c {threads} -t 0 -k 21 -s 1000 \
            {input.q:q} {input.db:q} -o {output:q}
    """

rule search_species_3216_singleclust:
    input:
        q="outputs.cds/cds3/species.cds3.singleclust.mf.csv",
        db="3216.manifest.csv",
    output:
        "outputs.cds/cds3/manysearch.cds3.singleclust.3216.csv",
    conda: "env-sourmash.yml"
    threads: 32
    shell: """
        sourmash scripts manysearch -c {threads} -t 0 -k 21 -s 1000 \
            {input.q:q} {input.db:q} -o {output:q}
    """

rule search_species_bw:
    input:
        q="outputs.cds/cds3/{species}.cds{c}.sig.zip",
        db=BW_ROCKSDB,
    output:
        "outputs.cds/cds3/outputs.branchwater.cds{c,^[\\.]+}/manysearch.{species}.csv",
    conda: "env-sourmash.yml"
    threads: 1
    shell: """
        sourmash scripts manysearch -c {threads} -t 0 -k 21 -s 1000 \
            {input.q:q} {input.db:q} -o {output:q}
    """

rule search_species_bw_singleclust:
    input:
        q="outputs.cds/cds3/{species}.cds3.singleclust.sig.zip",
        db=BW_ROCKSDB,
    output:
        "outputs.cds/cds3/outputs.branchwater.cds3.singleclust/manysearch.{species}.csv",
    conda: "env-sourmash.yml"
    threads: 1
    shell: """
        sourmash scripts manysearch -c {threads} -t 0 -k 21 -s 1000 \
            {input.q:q} {input.db:q} -o {output:q}
    """

rule search_species_bw_min:
    input:
        q="outputs.cds/cds3/outputs.minsig/{species}.cds{c}.min{m}.sig.zip",
        db=BW_ROCKSDB,
    output:
        "outputs.cds/cds3/outputs.branchwater.cds{c}.min{m,^[\\.]+}/manysearch.{species}.csv",
    conda: "env-sourmash.yml"
    threads: 1
    shell: """
        sourmash scripts manysearch -c {threads} -t 0 -k 21 -s 1000 \
            {input.q:q} {input.db:q} -o {output:q}
    """

rule search_species_bw_min_dedup:
    input:
        q="outputs.cds/cds3/{species}.cds{c}.min{m}.dedup.sig.zip",
        db=BW_ROCKSDB,
    output:
        "outputs.cds/cds3/outputs.branchwater.cds{c}.min{m}.dedup/manysearch.{species}.csv",
    conda: "env-sourmash.yml"
    threads: 1
    shell: """
        sourmash scripts manysearch -c {threads} -t 0 -k 21 -s 1000 \
            {input.q:q} {input.db:q} -o {output:q}
    """

rule screen_min90:
    input:
        sig='outputs.cds/cds3/{n}.sig.zip',
        metags='outputs.mapping/rand-metags.mf.csv',
    output:
        'outputs.cds/cds3/outputs.minsig/{n}.min{m}.sig.zip'
    shell: '''
        scripts/screen-sigs-x-metags.py {input.sig:q} {input.metags:q} \
            -m {wildcards.m} -o {output:q} -k 21
    '''
    
rule kmers_wc:
    input:
        sig='outputs.cds/cds3/outputs.minsig/{s}.cds3.min50.sig.zip',
        seqs='outputs./cds/{s}.cds.fa.gz',
    output:
        fa=touch('outputs.cds/cds3/outputs.minsig/{s}.cds3.min50.fa'),
    shell: """
        sourmash sig kmers --sig {input.sig:q} --sequences {input.seqs:q} \
            --save-sequences {output.fa:q} -k 21 || true
    """
        
    
rule cdhit_wc:
    input:
        fa='outputs.cds/cds3/outputs.minsig/{s}.cds3.min50.fa',
    output:
        fa='outputs.cds/cds3/outputs.minsig/{s}.cds3.min50.dedup.fa',
    threads: 16
    shell: """
        cd-hit -c 1.0 -i {input.fa:q} -o {output:q} -T {threads} -M 5000
    """
        
rule sketch_dedup:
    input:
        fa='outputs.cds/cds3/outputs.minsig/{s}.cds3.min50.dedup.fa',
    output:
        fa='outputs.cds/cds3/{s}.cds3.min50.dedup.sig.zip',
    shell: """
        sourmash sketch dna -p k=21,k=31,k=51 {input.fa:q} -o {output.fa:q} \
           --name {wildcards.s:q} -f
    """

rule sketch_singleclust:
    input:
        fa='outputs.cds/singleclust/{s}.cds3.min50.dedup.fa',
    output:
        fa='outputs.cds/cds3/{s}.cds3.singleclust.sig.zip',
    shell: """
        sourmash sketch dna -p k=21,k=31,k=51 {input.fa:q} -o {output.fa:q} \
           --name {wildcards.s:q} -f
    """

rule map_index_rand_cds:
    input:
        g="outputs.cds/singleclust/{s}.cds3.min50.dedup.fa",
        metag="outputs.mapping/bams.cds.rand/{m}.x.{s}.fastq.gz",
    output:
        bam='outputs.cds/singleclust/bam/{m}.x.{s}.bam',
        bai='outputs.cds/singleclust/bam/{m}.x.{s}.bam.bai',
    threads: 8
    conda: "env-mapping.yml"
    shell: """
        minimap2 -ax sr -t {threads} {input.g:q} {input.metag:q} | samtools view -b -F 4 - | samtools sort - > {output.bam:q}
        samtools index {output.bam:q}
    """

rule map_coverage_cds:
    input:
        bam='outputs.cds/singleclust/bam/{m}.x.{s}.bam',
        bai='outputs.cds/singleclust/bam/{m}.x.{s}.bam.bai',
        fa='outputs.cds/singleclust/{s}.cds3.min50.dedup.fa',
    output:
        'outputs.cds/singleclust/bam/{m}.x.{s}.coverage.txt'
    threads: 8
    conda: "env-mapping.yml"
    shell: """
        samtools coverage {input.bam:q} {input.fa:q} > {output:q}
    """

rule map_coverage_cds_depth:
    input:
        bam='outputs.cds/singleclust/bam/{m}.x.{s}.bam',
        bai='outputs.cds/singleclust/bam/{m}.x.{s}.bam.bai',
        fa='outputs.cds/singleclust/{s}.cds3.min50.dedup.fa',
    output:
        'outputs.cds/singleclust/bam/{m}.x.{s}.depth.txt'
    threads: 8
    conda: "env-mapping.yml"
    shell: """
        samtools depth -aa {input.bam:q} {input.fa:q} > {output:q}
    """

rule extract_species_genes:
    input:
        "outputs.cds/singleclust/species-genes.csv",
        # singleclust FA files
    output:
        directory("outputs.cds/cds3-genes/")
    shell: """
        scripts/split-genes-by-species.py outputs.cds/singleclust/species-genes.csv outputs.cds/singleclust/*.fa -o {output}
    """

rule sketch_species_genes:
    input:
        "outputs.cds/cds3-genes/{s}.genes.fa",
    output:
        "outputs.cds/cds3-genes/{s}.sig.zip",
    shell: """
        sourmash sketch dna {input:q} -o {output:q} --name {wildcards.s:q} \
           -p k=21,k=31,k=51
    """

rule genes_mf:
    input:
        expand("outputs.cds/cds3-genes/{s}.sig.zip", s=SPECIES_WITH_GENES),
    output:
        "outputs.cds/cds3-genes/cds3-genes.mf.csv",
    shell: """
        sourmash sig collect -F csv {input:q} -o {output:q}
    """

rule search_species_3216_genes:
    input:
        q="outputs.cds/cds3-genes/cds3-genes.mf.csv",
        db="3216.manifest.csv",
    output:
        "outputs.cds/cds3-genes/manysearch.3216.csv",
    conda: "env-sourmash.yml"
    threads: 32
    shell: """
        sourmash scripts manysearch -c {threads} -t 0 -k 21 -s 1000 \
            {input.q:q} {input.db:q} -o {output:q}
    """

rule search_species_bw_genes:
    input:
        q="outputs.cds/cds3-genes/{s}.sig.zip",
        db=BW_ROCKSDB,
    output:
        "outputs.cds/cds3-genes/outputs.branchwater/manysearch.{s}.csv",
    conda: "env-sourmash.yml"
    threads: 1
    shell: """
        sourmash scripts manysearch -c {threads} -t 0 -k 21 -s 1000 \
            {input.q:q} {input.db:q} -o {output:q}
    """
