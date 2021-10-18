#!/usr/bin/env nextflow

/*
  Copyright (c) 2021, Andre Kahles, ETH Zurich

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

nextflow.enable.dsl = 2
version = '0.1.0'  // package version

// universal params go here, change default value as needed
params.container = ""
params.container_registry = ""
params.container_version = ""
params.cpus = 1
params.mem = 1  // GB
params.publish_dir = ""  // set to empty string will disable publishDir

// tool specific parmas go here, add / change as needed
params.input_file = ""
params.cleanup = true

include { demoFastqc } from './wfpr_modules/github.com/icgc-tcga-pancancer/awesome-wfpkgs1/demo-fastqc@0.2.0/demo-fastqc' params([*:params, 'cleanup': false])
include { cleanupWorkdir; getSecondaryFiles; getBwaSecondaryFiles } from './wfpr_modules/github.com/icgc-argo/demo-wfpkgs/demo-utils@1.3.0/main.nf' params([*:params, 'cleanup': false])


// please update workflow code as needed
workflow RnaSeqAlignmentWf {
  take:  // update as needed
    input_file


  main:  // update as needed
    demoFastqc(input_file)
    if (params.cleanup) { cleanupWorkdir(demoFastqc.out, true) }

  emit:  // update as needed
    output_file = demoFastqc.out.output_file
}


// this provides an entry point for this main script, so it can be run directly without clone the repo
// using this command: nextflow run <git_acc>/<repo>/<pkg_name>/<main_script>.nf -r <pkg_name>.v<pkg_version> --params-file xxx
workflow {
  RnaSeqAlignmentWf(
    file(params.input_file)
  )
}