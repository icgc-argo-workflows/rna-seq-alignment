#!/usr/bin/env nextflow

/*
  Copyright (c) 2021, icgc-argo-workflows

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
    Andre Kahles
*/

/********************************************************************/
/* this block is auto-generated based on info from pkg.json where   */
/* changes can be made if needed, do NOT modify this block manually */
nextflow.enable.dsl = 2
version = '0.2.1'

container = [
    'ghcr.io': 'ghcr.io/icgc-argo-workflows/rna-seq-alignment.genome-alignment-hisat2'
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
params.index = "NO_FILE_1/NO_FILE_1" 
params.gtf = "NO_FILE_2" //input/test_annotation.gtf"
params.input_files = ["NO_FILE_3"]//input/TEST-PRO.donor1.donor1_sample1_id.sample_01.b22541e45ff72d9a042e877a0531af0b.lane.bam"
params.tempdir = ""

process icgcArgoRnaSeqAlignmentHISAT2 {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"
  publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: params.publish_dir

  cpus params.cpus
  memory "${params.mem} GB"

  input:  
    val index_base
    path index_parent
    path gtf
    path metadata
    val tempdir
    path input_files

  output:  // output, make update as needed
    path("*.hisat2.Aligned.out.bam"), emit: bam
    path("*novel_splicesites.txt"), emit: junctions
    path("*.hisat2.all_logs.supplement.tar.gz"), emit: logs

  script:
    def tempdir_arg = tempdir != "" ? "--tempdir ${tempdir}" : ""
    """
    python /tools/alignHISAT2.py \\
           --index ${index_base} \\
           --annotation ${gtf} \\
           --metadata ${metadata} \\
           --threads ${params.cpus} \\
           --input-files ${input_files} ${tempdir_arg}\\
           --mem ${params.mem * 1000} > align.log 2>&1
    """
}

// this provides an entry point for this main script, so it can be run directly without clone the repo
// using this command: nextflow run <git_acc>/<repo>/<pkg_name>/<main_script>.nf -r <pkg_name>.v<pkg_version> --params-file xxx
workflow {
  icgcArgoRnaSeqAlignmentHISAT2(
    params.index,
    file(params.index).getParent(),
    file(params.gtf),
    file(params.metadata),
    params.tempdir,
    params.input_files.collect({it -> file(it)}),
  )
}
