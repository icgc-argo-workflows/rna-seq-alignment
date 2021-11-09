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
    Andre Kahles
*/

/*
 This is an auto-generated checker workflow to test the generated main template workflow, it's
 meant to illustrate how testing works. Please update to suit your own needs.
*/

/********************************************************************/
/* this block is auto-generated based on info from pkg.json where   */
/* changes can be made if needed, do NOT modify this block manually */
nextflow.enable.dsl = 2
version = '0.2.0'

container = [
    'ghcr.io': 'ghcr.io/icgc-argo-workflows/rna-seq-alignment.genome-alignment-hisat2'
]
default_container_registry = 'ghcr.io'
/********************************************************************/

// universal params
params.container_registry = ""
params.container_version = ""
params.container = ""

// tool specific parmas go here, add / change as needed
params.index = "NO_FILE_1/NO_FILE_1"
params.gtf = "NO_FILE_2"
params.input_files = ["NO_FILE_3"]
params.input_format = "ubam"
params.sample = "sample_01"
params.pair_status = "paired"
params.expected_output = "tests/expected/sample_01_Aligned.out.bam"

include { icgcArgoRnaSeqAlignmentHISAT2 } from '../alignHISAT2' params(['cleanup': false, *:params]) 

process diff_bam {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"

  input:
    path output_file
    path expected_file

  output:
    stdout()

  script:
    """
    diff <(samtools view --no-PG ${output_file} | sort) <(samtools view --no-PG ${expected_file} | sort) \
      && ( echo "Test PASSED" && exit 0 ) || ( echo "Test FAILED, bam files mismatch." && exit 1 )
    """
}



workflow checker {
  take:
    index_base
    index_parent
    gtf
    input_files
    input_format
    sample
    pair_status
    expected_bam

  main:
    icgcArgoRnaSeqAlignmentHISAT2(
        index_base,
        index_parent,
        gtf,
        input_files,
        input_format,
        pair_status,
        sample,
    )

    diff_bam(
      icgcArgoRnaSeqAlignmentHISAT2.out.bam,
      expected_bam
    )

}


workflow {
  checker(
    params.index,
    file(params.index).getParent(),
    file(params.gtf),
    params.input_files.collect({it -> file(it)}),
    params.input_format,
    params.sample,
    params.pair_status,
    file(params.expected_bam),
  )
}
