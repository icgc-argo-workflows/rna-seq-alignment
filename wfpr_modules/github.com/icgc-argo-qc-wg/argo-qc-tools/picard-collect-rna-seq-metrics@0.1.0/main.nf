#!/usr/bin/env nextflow

/*
  Copyright (c) 2021, ICGC ARGO

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

  Authors:
    Linda Xiang
*/

/********************************************************************/
/* this block is auto-generated based on info from pkg.json where   */
/* changes can be made if needed, do NOT modify this block manually */
nextflow.enable.dsl = 2
version = '0.1.0'  // package version

container = [
    'ghcr.io': 'ghcr.io/icgc-argo-qc-wg/argo-qc-tools.picard-collect-rna-seq-metrics'
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
params.aligned_seq = "NO_FILE1"
params.ref_flat = "NO_FILE2"
params.strand = ""
params.ignore_seq = "NO_FILE3"
params.ribosomal_interval_list = "NO_FILE4"


process picardCollectRnaSeqMetrics {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"
  publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: params.publish_dir

  cpus params.cpus
  memory "${params.mem} GB"

  input:  // input, make update as needed
    path aligned_seq
    path ref_flat
    path ignore_seq
    path ribosomal_interval_list
    val strand

  output:  // output, make update as needed
    path "${aligned_seq}.collectrnaseqmetrics.tgz", emit: qc_tar

  script:
    // add and initialize variables here as needed
    arg_strand = strand == '' ? "" : " -s ${strand}"
    arg_ignore_seq = ignore_seq.name.startsWith('NO_FILE') ? "" : "-x ${ignore_seq}"
    arg_ribosomal_interval_list = ribosomal_interval_list.name.startsWith('NO_FILE') ? "" : "-b ${ribosomal_interval_list}"

    """
    main.py \
      -m ${(int) (params.mem * 1000)} \
      -i ${aligned_seq} \
      -r ${ref_flat} \
      ${arg_strand} \
      ${arg_ignore_seq} \
      ${arg_ribosomal_interval_list}
    """
}


// this provides an entry point for this main script, so it can be run directly without clone the repo
// using this command: nextflow run <git_acc>/<repo>/<pkg_name>/<main_script>.nf -r <pkg_name>.v<pkg_version> --params-file xxx
workflow {
  picardCollectRnaSeqMetrics(
    file(params.aligned_seq),
    file(params.ref_flat),
    file(params.ignore_seq),
    file(params.ribosomal_interval_list),
    params.strand
  )
}
