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
    Junjun Zhang
*/

/********************************************************************/
/* this block is auto-generated based on info from pkg.json where   */
/* changes can be made if needed, do NOT modify this block manually */
nextflow.enable.dsl = 2
version = '1.0.0.1'

container = [
    'ghcr.io': 'ghcr.io/icgc-argo-workflows/data-processing-utility-tools.cleanup-workdir'
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
params.files_to_delete = "NO_FILE"
params.virtual_dep_flag = true  // default to true, ie, dep is always satisfied


process cleanupWorkdir {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"
  publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: params.publish_dir

  cpus params.cpus
  memory "${params.mem} GB"

  input:
    path files_to_delete  // more accurately, other non-hidden files in the same folder will be deleted as well
    val virtual_dep_flag  // for specifying steps do not produce output files but produce values, set those values here

  output:
    stdout

  script:
    """
    set -euxo pipefail

    IFS=" "
    read -a files <<< "${files_to_delete}"
    for f in "\${files[@]}"
    do
        dir_to_rm=\$(dirname \$(readlink -f \$f))

        if [[ \$dir_to_rm != ${workflow.workDir}/* ]]; then  # skip dir not under workdir, like from input file dir
            echo "Not delete: \$dir_to_rm/*\"
            continue
        fi

        rm -fr \$dir_to_rm/*  # delete all files and subdirs but not hidden ones
        echo "Deleted: \$dir_to_rm/*"
    done
    """
}


// this provides an entry point for this main script, so it can be run directly without clone the repo
// using this command: nextflow run <git_acc>/<repo>/<pkg_name>/<main_script>.nf -r <pkg_name>.v<pkg_version> --params-file xxx
workflow {
  cleanupWorkdir(
    Channel.fromPath(params.files_to_delete),
    params.virtual_dep_flag
  )
}
