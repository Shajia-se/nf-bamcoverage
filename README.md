# nf-bamcoverage

`nf-bamcoverage` creates normalized BigWig tracks from `nf-chipfilter` BAM files using deepTools `bamCoverage`.

## Inputs

Recommended:

```bash
--samples_master /path/to/samples_master.csv
--bam_input_dir /path/to/chipfilter_output
```

Alternative:

- `--bam_list`: TSV with `sample_id<TAB>bam_path`
- `--bam_pattern`: direct BAM glob

## Outputs

```text
${project_folder}/${bamcoverage_output}/bigwig/
  sample.bw
  sample.bamCoverage.log.txt
```

## Run

```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --bam_input_dir /path/to/chipfilter_output \
  --project_folder /path/to/output_project
```

Actual execution should be tested where Nextflow is installed.
