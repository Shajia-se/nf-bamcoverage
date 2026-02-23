#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process bamcoverage {

  tag { "${sample_id}" }

  publishDir "${params.project_folder}/${params.bamcoverage_output}/bigwig", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(bam)

  output:
    tuple val(sample_id), path("${sample_id}.bw"), emit: bw
    tuple val(sample_id), path("${sample_id}.bamCoverage.log.txt"), emit: log

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
  def ch_bams

  if (params.bam_list) {
    ch_bams = Channel
      .fromPath(params.bam_list, checkIfExists: true)
      .splitCsv(sep: '\t', header: false)
      .map { row ->
        def sample_id = row[0].toString().trim()
        def bam_path  = row[1].toString().trim()
        tuple(sample_id, file(bam_path))
      }
  } else if (params.samples_master) {
    def master = file(params.samples_master)
    assert master.exists() : "samples_master not found: ${params.samples_master}"

    def header = null
    def records = []
    master.eachLine { line, n ->
      if (!line?.trim()) return
      def cols = line.split(',', -1)*.trim()
      if (n == 1) {
        header = cols
      } else {
        def rec = [:]
        header.eachWithIndex { h, i -> rec[h] = i < cols.size() ? cols[i] : '' }
        records << rec
      }
    }

    assert header : "samples_master header not found: ${params.samples_master}"
    assert header.contains('sample_id') : "samples_master missing required column: sample_id"

    def isEnabled = { rec ->
      def v = rec.enabled?.toString()?.trim()?.toLowerCase()
      (v == null || v == '' || v == 'true')
    }
    def isControl = { rec ->
      rec.is_control?.toString()?.trim()?.toLowerCase() == 'true'
    }

    def includeControls = (params.bamcoverage_include_controls == null) ? true : params.bamcoverage_include_controls
    def bamDir = file(params.bam_input_dir)
    assert bamDir.exists() : "bam_input_dir not found: ${params.bam_input_dir}"

    def rows = records
      .findAll { rec -> isEnabled(rec) }
      .findAll { rec -> includeControls ? true : !isControl(rec) }
      .collect { rec ->
        def sid = rec.sample_id?.toString()?.trim()
        if (!sid) return null
        def hits = bamDir.listFiles()?.findAll { f ->
          f.isFile() && f.name.endsWith('.clean.bam') && (f.name == "${sid}.clean.bam" || f.name.startsWith("${sid}_"))
        } ?: []
        if (hits.isEmpty()) throw new IllegalArgumentException("No clean BAM found for sample_id '${sid}' under: ${params.bam_input_dir}")
        if (hits.size() > 1) throw new IllegalArgumentException("Multiple clean BAM files matched sample_id '${sid}': ${hits*.name.join(', ')}")
        [sid: sid, bam: file(hits[0].toString())]
      }
      .findAll { it != null }

    ch_bams = Channel
      .fromList(rows)
      .ifEmpty { error "No BAMs generated from samples_master: ${params.samples_master}" }
      .map { r -> tuple(r.sid, r.bam) }
  } else {
    ch_bams = Channel
      .fromPath(params.bam_pattern, checkIfExists: true)
      .ifEmpty { error "No BAM files found with --bam_pattern: ${params.bam_pattern}" }
      .map { bam ->
        tuple(bam.baseName, bam)
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
