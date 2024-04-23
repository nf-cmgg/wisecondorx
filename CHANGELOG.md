# nf-cmgg/wisecondorx: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0dev

### New features

1. The metrics are now in the multiqc report instead of in a separate file

### Changes

1. Updated the pipeline template to nf-core v2.13.1
2. Migrate from nf-validation to nf-schema
3. The `--bin_sizes` parameter now only takes kilobases as input (bases are no longer allowed)

## v1.1.0 - Naive Junior - [15 Sep 2023]

### New features

1. Added a new parameter `--bin_sizes` that takes a comma-delimited list of bin sizes to create references for. This will make it possible to create references for multiple bin sizes at once.

## v1.0.1 - Helpful Apprentice - [19 June 2023]

### New features

1. Added a gender metrics output file (can be turned off by using `--no_metrics`)
2. Added a new field to the samplesheet (`gender`). This optional field is used for the generation of the metrics file. The gender will be determined by the pipeline if this field is empty.

### Improvements

1. Added parameter information to the docs

## v1.0.0 - Clumsy Apprentice - [22 May 2023]

Initial release of nf-cmgg/wisecondorx, created with the [nf-core](https://nf-co.re/) template.
