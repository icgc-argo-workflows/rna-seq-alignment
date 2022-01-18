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
*/

/********************************************************************/
/* this block is auto-generated based on info from pkg.json where   */
/* changes can be made if needed, do NOT modify this block manually */
nextflow.enable.dsl = 2
version = '0.2.0'

container = [
    'ghcr.io': 'ghcr.io/icgc-argo-workflows/data-processing-utility-tools.json-parser'
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
params.metadata_analysis = ""


process jsonParser {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"
  publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: params.publish_dir

  cpus params.cpus
  memory "${params.mem} GB"

  input:
    path metadata_analysis

  output:
    env STUDY_ID, emit: study_id
    env DONOR_ID, emit: donor_id
    env EXP, emit: experimental_strategy
    env PAIRED, emit: paired
    env ANALYSIS_TOOLS, emit: analysis_tools
    env STRAND, emit: library_strandedness

  script:
    """
    set -euxo pipefail
    VARIABLE1=`cat ${metadata_analysis} | jq -r 'if ([.read_groups[]?] | length) >0 then [.read_groups[] | .is_paired_end] | all | tostring else "NULL" end' | tr -d '\\n'`
    PAIRED=\${VARIABLE1:-'NULL'}
    VARIABLE2=`cat ${metadata_analysis} | jq -r '[.files[] | .info? | .analysis_tools[]?] | unique | join(",")' | tr -d '\\n'`
    ANALYSIS_TOOLS=\${VARIABLE2:-'NULL'}
    VARIABLE3=`cat ${metadata_analysis} | jq -r 'if (.experiment | .library_strandedness?) then .experiment|.library_strandedness else "NULL" end' | tr -d '\\n'`
    STRAND=\${VARIABLE3:-'NULL'}
    STUDY_ID=`cat ${metadata_analysis} | jq -er '.studyId' | tr -d '\\n'`
    DONOR_ID=`cat ${metadata_analysis} | jq -er '.samples[0].donor.donorId' | tr -d '\\n'`
    EXP=`cat ${metadata_analysis} | jq -er '.experiment | .experimental_strategy?  // .library_strategy' | tr -d '\\n'`
    """
}


// this provides an entry point for this main script, so it can be run directly without clone the repo
// using this command: nextflow run <git_acc>/<repo>/<pkg_name>/<main_script>.nf -r <pkg_name>.v<pkg_version> --params-file xxx
workflow {
  jsonParser(
    file(params.metadata_analysis)
  )
}
