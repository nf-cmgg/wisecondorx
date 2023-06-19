# CenterForMedicalGeneticsGhent/nf-cmgg-wisecondorx: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

```
results/
├── multiqc_reports
├── pipeline_info
├── metrics.txt
└── WisecondorX_09052023.npz
```

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- SampleGender (ngs-bits) - Determine the gender from a BAM/CRAM file
- WiseCondorX convert - Convert CRAM/BAM files to NPZ files
- WiseCondorX newref- Create a new reference from all NPZ files
- MultiQC- Aggregate report describing results and QC from the whole pipeline
- Pipeline information - Report metrics generated during the workflow execution
