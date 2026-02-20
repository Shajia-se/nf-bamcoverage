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
    def blacklist_arg = params.blacklist ? "--blackListFileName ${params.blacklist}" : ""
    def extend_arg    = params.extendReads ? "--extendReads ${params.extendReads}" : ""
    def ignore_dup    = params.ignoreDuplicates ? "--ignoreDuplicates" : ""
    def center_arg    = params.centerReads ? "--centerReads" : ""
    def smooth_arg    = params.smoothLength ? "--smoothLength ${params.smoothLength}" : ""
    def min_mapq_arg  = params.minMappingQuality ? "--minMappingQuality ${params.minMappingQuality}" : ""
    def sam_flag_arg  = params.samFlagExclude ? "--samFlagExclude ${params.samFlagExclude}" : ""
    def scale_arg     = params.scaleFactor ? "--scaleFactor ${params.scaleFactor}" : ""
    def extra         = params.extra_args ?: ""

    """
    set -euo pipefail

    bamCoverage \\
      --bam ${bam} \\
      --outFileName ${sample_id}.bw \\
      --outFileFormat bigwig \\
      --binSize ${params.binSize} \\
      --normalizeUsing ${params.normalizeUsing} \\
      --effectiveGenomeSize ${params.effectiveGenomeSize} \\
      ${extend_arg} \\
      ${ignore_dup} \\
      ${center_arg} \\
      ${smooth_arg} \\
      ${min_mapq_arg} \\
      ${sam_flag_arg} \\
      ${scale_arg} \\
      ${blacklist_arg} \\
      --numberOfProcessors ${task.cpus} \\
      ${extra} \\
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
