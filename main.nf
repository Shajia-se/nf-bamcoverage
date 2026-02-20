#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process bamcoverage {

  tag { "${sample_id}" }

  publishDir "${params.project_folder}/${params.bamcoverage_output}/bigwig", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(bam)

  output:
    tuple val(sample_id), path("${sample_id}.bw"), emit: bw
    tuple val(sample_id), path("${sample_id}.bamcoverage.log.txt"), emit: log

  script:
    def optArgs = []
    if (params.extendReads)       optArgs << "--extendReads ${params.extendReads}"
    if (params.ignoreDuplicates)  optArgs << "--ignoreDuplicates"
    if (params.centerReads)       optArgs << "--centerReads"
    if (params.smoothLength)      optArgs << "--smoothLength ${params.smoothLength}"
    if (params.minMappingQuality) optArgs << "--minMappingQuality ${params.minMappingQuality}"
    if (params.samFlagExclude)    optArgs << "--samFlagExclude ${params.samFlagExclude}"
    if (params.scaleFactor)       optArgs << "--scaleFactor ${params.scaleFactor}"
    if (params.blacklist)         optArgs << "--blackListFileName ${params.blacklist}"
    if (params.extra_args)        optArgs << params.extra_args.toString().trim()
    def optArgString = optArgs.join(' ')

    """
    set -euo pipefail
    mkdir -p tmp
    export TMPDIR=\$PWD/tmp
    export TEMP=\$PWD/tmp
    export TMP=\$PWD/tmp
    export MPLCONFIGDIR=\$PWD/tmp/matplotlib
    mkdir -p "\$MPLCONFIGDIR"

    bam_local="${bam}"
    bam_abs=\$(readlink -f "\$bam_local" || echo "\$bam_local")

    # Ensure bamCoverage can see a BAM index in work dir
    if [[ -f "\${bam_local}.bai" ]]; then
      :
    elif [[ -f "\${bam_local%.bam}.bai" ]]; then
      ln -sf "\${bam_local%.bam}.bai" "\${bam_local}.bai"
    elif [[ -f "\${bam_abs}.bai" ]]; then
      ln -sf "\${bam_abs}.bai" "\${bam_local}.bai"
    elif [[ -f "\${bam_abs%.bam}.bai" ]]; then
      ln -sf "\${bam_abs%.bam}.bai" "\${bam_local}.bai"
    else
      echo "ERROR: BAM index not found for \$bam_local" >&2
      exit 1
    fi

    bamCoverage \
      --bam ${bam} \
      --outFileName ${sample_id}.bw \
      --outFileFormat bigwig \
      --binSize ${params.binSize} \
      --normalizeUsing ${params.normalizeUsing} \
      --effectiveGenomeSize ${params.effectiveGenomeSize} \
      --numberOfProcessors ${task.cpus} \
      ${optArgString} \
      &> ${sample_id}.bamCoverage.log.txt

    ls -lh ${sample_id}.bw 1>&2
    """
}

workflow {
  Channel
    .fromPath(params.bam_pattern, checkIfExists: true)
    .ifEmpty { error "No BAM files found with --bam_pattern: ${params.bam_pattern}" }
    .map { bam ->
      tuple(bam.baseName, bam)
    }
    .set { ch_bams_from_pattern }

  def ch_bams = ch_bams_from_pattern

  if (params.bam_list) {
    ch_bams = Channel
      .fromPath(params.bam_list, checkIfExists: true)
      .splitCsv(sep: '\t', header: false)
      .map { row ->
        def sample_id = row[0].toString().trim()
        def bam_path  = row[1].toString().trim()
        tuple(sample_id, file(bam_path))
      }
  }

  if (params.require_bai) {
    ch_bams = ch_bams.map { sample_id, bam ->
      def bai1 = file("${bam}.bai")
      def bai2 = file(bam.toString().replaceFirst(/\.bam$/, '.bai'))
      if (!bai1.exists() && !bai2.exists()) {
        error "Missing BAM index for ${bam}. Provide .bai or set --require_bai false."
      }
      tuple(sample_id, bam)
    }
  }

  bamcoverage(ch_bams)
}
