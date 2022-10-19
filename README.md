# mev-sctk-mast

This repository contains a WDL-format Cromwell-compatible workflow for executing a differential expression analysis on single-cell RNA-seq data using the MAST tool (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5) as provided through the Single-Cell Toolkit (https://github.com/compbiomed/singleCellTK).

To use, simply fill in the the `inputs.json` with the various inputs and submit to a Cromwell runner. 

Alternatively, you can pull the docker image (https://github.com/web-mev/mev-sctk-mast/pkgs/container/mev-sctk-mast), start the container, and run: 

```
Rscript /opt/software/mast_dge.R \
    -f <path to raw counts tab-delimited file> \
    -o <prefix for output file (string)> \
    -a <"experimental" samples as comma-delimited string> \
    -b <(Optional) "base" samples as comma-delimited string> \
    --experimental_group_name <experimental group name (string)> \
    --base_group_name < (Optional) base group name (string)>
```

If you are performing a basic contrast between two groups of samples/cells, then specify `-a` and `-b` as comma delimited strings (no spaces). If you would like to perform a "biomarker-style" analysis of comparing one group of samples/cells versus the rest, omit the `-b` and `--base_group_name` arguments.