# nf-bamcoverage

Generate bigWig files from BAM using deepTools `bamCoverage`.

## Inputs

Provide one of:

- `--bam_pattern` (glob for BAM files)
- `--bam_list` (TSV: `sample_id<TAB>bam_path`)

Optional:

- `--blacklist` (BED)
- `--normalizeUsing` (default: `RPGC`)
- `--effectiveGenomeSize` (default: `2468088461` for mm39)
- `--binSize` (default: `20`)
- `--minMappingQuality` (default: `10`)

## Defaults (aligned with your previous scripts)

- `--extendReads 150`
- `--centerReads`
- `--ignoreDuplicates`
- `--smoothLength 60`

## Output

Output directory: `${project_folder}/${bamcoverage_output}/bigwig`

Per sample:

- `<sample>.bw`
- `<sample>.bamcoverage.log.txt`

## Run

```bash
nextflow run main.nf -profile hpc \
  --bam_pattern '/ictstr01/groups/idc/projects/uhlenhaut/jiang/pipelines/nf-chipfilter/chipfilter_output/*.clean.bam'
```

With explicit output folder:

```bash
nextflow run main.nf -profile hpc \
  --bam_pattern '/path/*.clean.bam' \
  --bamcoverage_output bamcoverage_output
```

Resume:

```bash
nextflow run main.nf -profile hpc -resume
```
