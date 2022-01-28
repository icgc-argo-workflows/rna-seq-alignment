#!/usr/bin/env nextflow

/*
  Copyright (C) 2021,  Ontario Institute for Cancer Research
  
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
    Linda Xiang
    Junjun Zhang
*/

/********************************************************************/
/* this block is auto-generated based on info from pkg.json where   */
/* changes can be made if needed, do NOT modify this block manually */
nextflow.enable.dsl = 2
version = '0.5.0.1'

container = [
    'ghcr.io': 'ghcr.io/icgc-argo-workflows/data-processing-utility-tools.payload-gen-seq-experiment'
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
params.experiment_info_tsv = "NO_FILE1"
params.read_group_info_tsv = "NO_FILE2"
params.file_info_tsv = "NO_FILE3"
params.extra_info_tsv = "NO_FILE4"


process payloadGenSeqExperiment {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"
  publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: params.publish_dir

  cpus params.cpus
  memory "${params.mem} GB"

  input:
    path experiment_info_tsv
    path read_group_info_tsv
    path file_info_tsv
    path extra_info_tsv

  output:
    path "*.sequencing_experiment.payload.json", emit: payload

  script:
    args_experiment_info_tsv = !experiment_info_tsv.name.startsWith("NO_FILE") ? "-x ${experiment_info_tsv}" : ""
    args_read_group_info_tsv = !read_group_info_tsv.name.startsWith("NO_FILE") ? "-r ${read_group_info_tsv}" : ""
    args_file_info_tsv = !file_info_tsv.name.startsWith("NO_FILE") ? "-f ${file_info_tsv}" : ""
    args_extra_info_tsv = !extra_info_tsv.name.startsWith("NO_FILE") ? "-e ${extra_info_tsv}" : ""

    """
    main.py \
         ${args_experiment_info_tsv} \
         ${args_read_group_info_tsv} \
         ${args_file_info_tsv} \
         ${args_extra_info_tsv}
    """
}


// this provides an entry point for this main script, so it can be run directly without clone the repo
// using this command: nextflow run <git_acc>/<repo>/<pkg_name>/<main_script>.nf -r <pkg_name>.v<pkg_version> --params-file xxx
workflow {
  payloadGenSeqExperiment(
    file(params.experiment_info_tsv),
    file(params.read_group_info_tsv),
    file(params.file_info_tsv),
    file(params.extra_info_tsv)
  )
}
