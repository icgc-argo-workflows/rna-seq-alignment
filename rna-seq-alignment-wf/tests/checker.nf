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

/*
 This is an auto-generated checker workflow to test the generated main template workflow, it's
 meant to illustrate how testing works. Please update to suit your own needs.
*/

nextflow.enable.dsl = 2
version = '0.1.0'  // package version

// universal params
params.publish_dir = ""
params.container = ""
params.container_registry = ""
params.container_version = ""

// tool specific parmas go here, add / change as needed
params.study_id = ""
params.analysis_id = ""
params.ref_genome_fa = ""
params.ref_genome_gtf = ""

params.analysis_metadata = "NO_FILE"
params.experiment_info_tsv = "NO_FILE1"
params.read_group_info_tsv = "NO_FILE2"
params.file_info_tsv = "NO_FILE3"
params.extra_info_tsv = "NO_FILE4"
params.sequencing_files = []


include { RnaSeqAlignmentWf } from '../main'

workflow {
  main:
    RnaSeqAlignmentWf(
      params.study_id,
      params.analysis_id,
      params.ref_genome_fa,
      params.ref_genome_gtf,
      params.analysis_metadata,
      params.experiment_info_tsv,
      params.read_group_info_tsv,
      params.file_info_tsv,
      params.extra_info_tsv,
      params.sequencing_files
    )
}
