# 2025-workflow-core99

Workflow for analysis of matches to pig core microbiome @ 99%.

## Mapping workflow

For Figure 1:

```
snakemake --use-conda -j 8 do_mapping
```

All outputs will be placed in `outputs.mapping/`.
