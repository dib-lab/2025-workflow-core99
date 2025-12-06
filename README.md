# 2025-workflow-core99

Workflow for analysis of matches to pig core microbiome @ 99%.

## 'genome-grist' inputs for mapping workflow

The mapping workflow takes two batches of metagenomes
(see `inputs.mapping/lowest-detection.metags.txt` and
`inputs.mapping/rand_subset.3216.100.txt`) and maps them against the
18 pangenomes in `inputs.mapping/names.list`. These metagenomes are
downloaded and trimmed and sketched with genome-grist. The config
files for genome grist are in `inputs.grist/conf-lowest-detection.yml`
and `inputs.grist/conf-rand-100-core.yml`, and you can run them with:

```
genome grist run conf-rand-100-core.yml smash_reads -j 24
```

## Mapping workflow

For Figure 1:

```
snakemake --use-conda -j 8 do_mapping
```

All outputs will be placed in `outputs.mapping/`.

Notes:
- you'll need about 40 GB of RAM to map to the largest pangenome,
  'Phocaeicola vulgatus', with 8 threads. So if you run with -j 24,
  you'll need to allocate 120 GB of RAM total.
