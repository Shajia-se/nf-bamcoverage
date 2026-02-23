# nf-bamcoverage

Generate bigWig files from BAM using deepTools `bamCoverage`.

## Input Modes (Priority Order)

1. `--bam_list` (TSV: `sample_id<TAB>bam_path`)
2. `--samples_master` auto mode
3. `--bam_pattern` glob mode

## Mode 1: `bam_list`

TSV format (no header):
```text
sample_id<TAB>/path/to/sample.clean.bam
```

## Mode 2: Auto from `samples_master`

Required column:
```text
sample_id
```
Optional columns used:
```text
is_control,enabled
```

Auto behavior:
- resolve BAM from `${bam_input_dir}/${sample_id}*.clean.bam`
- include/exclude controls via `bamcoverage_include_controls` (default `true`)

## Mode 3: `bam_pattern`

Example:
```bash
--bam_pattern '/path/to/*.clean.bam'
```

## Key Parameters

- `normalizeUsing` (default: `RPGC`)
- `effectiveGenomeSize` (default: `2468088461` for mm39)
- `binSize` (default: `20`)
- `minMappingQuality` (default: `10`)
- defaults aligned with previous scripts:
  - `extendReads=150`
  - `centerReads=true`
  - `ignoreDuplicates=true`
  - `smoothLength=60`

## Output

Output directory: `${project_folder}/${bamcoverage_output}/bigwig`

Per sample:
- `<sample>.bw`
- `<sample>.bamCoverage.log.txt`

## Run

Auto from `samples_master`:
```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --bam_input_dir /path/to/nf-chipfilter/chipfilter_output
```

Pattern mode:
```bash
nextflow run main.nf -profile hpc \
  --bam_pattern '/path/*.clean.bam'
```

Resume:
```bash
nextflow run main.nf -profile hpc -resume
```
