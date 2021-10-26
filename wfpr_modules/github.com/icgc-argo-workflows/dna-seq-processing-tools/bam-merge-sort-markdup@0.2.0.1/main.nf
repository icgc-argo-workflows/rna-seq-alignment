#!/usr/bin/env nextflow

/*
  Copyright (C) 2021,  icgc-argo

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Authors:
    Junjun Zhang
    Linda Xiang
*/

/********************************************************************/
/* this block is auto-generated based on info from pkg.json where   */
/* changes can be made if needed, do NOT modify this block manually */
nextflow.enable.dsl = 2
version = '0.2.0.1'

container = [
    'ghcr.io': 'ghcr.io/icgc-argo-workflows/dna-seq-processing-tools.bam-merge-sort-markdup'
]
default_container_registry = 'ghcr.io'
/********************************************************************/


// universal params go here
params.container_registry = ""
params.container_version = ""
params.container = ""

params.cpus = 1
params.mem = 1  // GB
params.publish_dir = ""  // set to empty string will disable publishDir


// tool specific parmas go here, add / change as needed
params.aligned_lane_bams = ""
params.ref_genome_gz = ""
params.aligned_basename = "grch38-aligned.merged"
params.markdup = true
params.output_format = "cram"
params.lossy = false
params.tempdir = "NO_DIR"

include { getSecondaryFiles } from './wfpr_modules/github.com/icgc-argo/data-processing-utility-tools/helper-functions@1.0.1/main'

process bamMergeSortMarkdup {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"
  publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: params.publish_dir

  cpus params.cpus
  memory "${params.mem} GB"

  input:
    path aligned_lane_bams
    path ref_genome_gz
    path ref_genome_gz_secondary_file
    val tempdir

  output:
    path "${params.aligned_basename}.{bam,cram}", emit: merged_seq
    path "${params.aligned_basename}.{bam.bai,cram.crai}", emit: merged_seq_idx
    path "${params.aligned_basename}.duplicates_metrics.tgz", optional: true, emit: duplicates_metrics

  script:
    arg_markdup = params.markdup ? "-d" : ""
    arg_lossy = params.lossy ? "-l" : ""
    arg_tempdir = tempdir != 'NO_DIR' ? "-t ${tempdir}" : ""
    """
    main.py \
      -i ${aligned_lane_bams} \
      -r ${ref_genome_gz} \
      -n ${params.cpus} \
      -b ${params.aligned_basename} ${arg_markdup} \
      -o ${params.output_format} ${arg_lossy} ${arg_tempdir}
    """
}


// this provides an entry point for this main script, so it can be run directly without clone the repo
// using this command: nextflow run <git_acc>/<repo>/<pkg_name>/<main_script>.nf -r <pkg_name>.v<pkg_version> --params-file xxx
workflow {
  bamMergeSortMarkdup(
    Channel.fromPath(params.aligned_lane_bams, checkIfExists: true).collect(),
    file(params.ref_genome_gz),
    Channel.fromPath(getSecondaryFiles(params.ref_genome_gz, ['fai', 'gzi']), checkIfExists: true).collect(),
    params.tempdir
  )
}
