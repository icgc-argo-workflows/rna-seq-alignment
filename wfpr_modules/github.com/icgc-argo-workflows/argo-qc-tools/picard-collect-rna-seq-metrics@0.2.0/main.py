#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
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
"""

import os
import sys
import argparse
import subprocess
import csv
import json
from glob import glob
import tarfile
import re

def run_cmd(cmd):
    proc = subprocess.Popen(
                cmd,
                shell=True,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
    stdout, stderr = proc.communicate()

    return (
        stdout.decode("utf-8").strip(),
        stderr.decode("utf-8").strip(),
        proc.returncode
    )

def get_tool_version(toolcmd):
    get_tool_version_cmd = f"{toolcmd} --version 2>&1 | cat"
    stdout, stderr, returncode = run_cmd(get_tool_version_cmd)
    if returncode:
        sys.exit(f"Error: unable to get version info for {toolcmd}.\nStdout: {stdout}\nStderr: {stderr}\n")

    return stdout.strip().split(':')[-1]

def prep_qc_metrics(output_dir, tool_ver):
    qc_metrics = {
        'tool': {
            'name': 'Picard:CollectRnaSeqMetrics',
            'version': tool_ver
        },
        'metrics': {},
        'description': 'RNA alignment metrics for a SAM or BAM file.'
    }

    collect_fields = {
      'PF_BASES': ['pf_bases', 'integer'],
      'PF_ALIGNED_BASES': ['pf_aligned_bases', 'integer'],
      'RIBOSOMAL_BASES': ['ribosomal_bases', 'integer'],
      'CODING_BASES': ['coding_bases', 'integer'],
      'UTR_BASES': ['utr_bases', 'integer'],
      'INTRONIC_BASES': ['intronic_bases', 'integer'],
      'INTERGENIC_BASES': ['intergenic_bases', 'integer'],
      'IGNORED_READS': ['ignored_reads', 'integer'],
      'CORRECT_STRAND_READS': ['correct_strand_reads', 'integer'],
      'INCORRECT_STRAND_READS': ['incorrect_strand_reads', 'integer'],
      'NUM_R1_TRANSCRIPT_STRAND_READS': ['num_r1_transcript_strand_reads', 'integer'],
      'NUM_R2_TRANSCRIPT_STRAND_READS': ['num_r2_transcript_strand_reads', 'integer'],
      'NUM_UNEXPLAINED_READS': ['num_unexplained_reads', 'integer'],
      'PCT_R1_TRANSCRIPT_STRAND_READS': ['pct_r1_transcript_strand_reads', 'float'],
      'PCT_R2_TRANSCRIPT_STRAND_READS': ['pct_r2_transcript_strand_reads', 'float'],
      'PCT_RIBOSOMAL_BASES': ['pct_ribosomal_bases', 'float'],
      'PCT_CODING_BASES': ['pct_coding_bases', 'float'],
      'PCT_UTR_BASES': ['pct_utr_bases', 'float'],
      'PCT_INTRONIC_BASES': ['pct_intronic_bases', 'float'],
      'PCT_INTERGENIC_BASES': ['pct_intergenic_bases', 'float'],
      'PCT_MRNA_BASES': ['pct_mrna_bases', 'float'],
      'PCT_USABLE_BASES': ['pct_usable_bases', 'float'],
      'PCT_CORRECT_STRAND_READS': ['pct_correct_strand_reads', 'float'],
      'MEDIAN_CV_COVERAGE': ['median_cv_coverage', 'float'],
      'MEDIAN_5PRIME_BIAS': ['median_5prime_bias', 'float'],
      'MEDIAN_3PRIME_BIAS': ['median_3prime_bias', 'float'],
      'MEDIAN_5PRIME_TO_3PRIME_BIAS': ['median_5prime_to_3prime_bias', 'float']
    }

    with open(output_dir+'/rna_metrics.txt', 'r') as mytext:
      include = False
      for line in mytext:
        if line.startswith('## METRICS CLASS'):
          include = True
          continue
        elif line.startswith('## HISTOGRAM'):
          include = False
          continue
        elif include:
          if not line.rstrip(): continue
          if line.startswith('PF_BASES'):
            header = line.rstrip().split('\t') 
            continue
          cols = line.rstrip().split('\t')
          break  # should only have one data row for ALL READS
    
    for h, c in zip(header, cols):
      if h not in collect_fields: continue
      if collect_fields[h][1] == 'string':
        field_value = str(c) if c else None
      elif collect_fields[h][1] == 'float':
        field_value = float(c) if c else None
      else: 
        field_value = int(c) if c else None
      qc_metrics['metrics'].update({
        collect_fields[h][0]: field_value
      })

    qc_metrics_file = 'qc_metrics.json'
    with open(qc_metrics_file, "w") as j:
        j.write(json.dumps(qc_metrics, indent=2))

    return qc_metrics_file


def prepare_tarball(seq, qc_metrics, output_dir):

    files_to_tar = [qc_metrics]
    for f in sorted(glob(output_dir+'/*')):
      files_to_tar.append(f)

    tarfile_name = f'{os.path.basename(seq)}.collectrnaseqmetrics.tgz'
    with tarfile.open(tarfile_name, "w:gz") as tar:
      for f in files_to_tar:
        tar.add(f, arcname=os.path.basename(f))


def main():
    """
    Python implementation of tool: picard-collect-rna-seq-metrics
    """

    parser = argparse.ArgumentParser(description='Tool: picard-collect-rna-seq-metrics')
    parser.add_argument('-m', '--jvm-mem', dest='jvm_mem', type=int, default=1000,
                        help='JVM max memory in MB')
    parser.add_argument('-i', '--aligned_seq', dest='aligned_seq', type=str,
                        help='Input SAM or BAM file.', required=True)
    parser.add_argument('-r', '--ref_flat', dest='ref_flat', type=str, required=True,
                        help='Gene annotations in refFlat form.')
    parser.add_argument('-s', '--strand', dest='strand', type=str, default='Unstranded', choices=['First_Stranded', 'Second_Stranded', 'Unstranded'],
                        help='For strand-specific library prep.')
    parser.add_argument('-x', '--ignore_seq', dest='ignore_seq', type=str,
                        help='If a read maps to a sequence specified with this option, all the bases in the read are counted as ignored bases. These reads are not counted as.')
    parser.add_argument('-b', '--ribosomal_interval_list', dest='ribosomal_interval_list', type=str,
                        help='Location of rRNA sequences in genome in interval_list format.')
    parser.add_argument('-t', '--tempdir', dest='tempdir', type=str,
                        help='Specify directory for temporary files')

    args = parser.parse_args()

    if not os.path.isfile(args.aligned_seq):
        sys.exit('Error: specified aligned seq %s does not exist or is not accessible!' % args.aligned_seq)

    if not os.path.isfile(args.ref_flat):
        sys.exit('Error: specified refFalt format gene annotation %s does not exist or is not accessible!' % args.ref_flat)

    if args.strand == 'First_Stranded':
      strand = 'SECOND_READ_TRANSCRIPTION_STRAND'
    elif args.strand == 'Second_Stranded':
      strand = 'FIRST_READ_TRANSCRIPTION_STRAND'
    elif args.strand == 'Unstranded':
      strand = 'NONE'


    jvm_Xmx = f'-Xmx{args.jvm_mem}M'
    # get version info
    toolcmd = f'java {jvm_Xmx} -jar /tools/picard.jar CollectRnaSeqMetrics'
    tool_ver = get_tool_version(toolcmd)

    # create output dir if not exist
    output_dir = 'output'
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)

    # run tool commands
    tool_args = [
        '--INPUT', args.aligned_seq,
        '--OUTPUT', output_dir+'/rna_metrics.txt',
        '--CHART_OUTPUT', output_dir+'/rna_metrics.pdf',
        '--REF_FLAT', args.ref_flat,
        '--STRAND_SPECIFICITY', strand
    ]

    if args.tempdir:
      tool_args += ['--TMP_DIR', args.tempdir]

    if args.ribosomal_interval_list:
      tool_args += ['--RIBOSOMAL_INTERVALS', args.ribosomal_interval_list]
    
    if args.ignore_seq:
      tool_args += ['--IGNORE_SEQUENCE', args.ignore_seq]

    cmd = [toolcmd] + tool_args
    stdout, stderr, returncode = run_cmd(" ".join(cmd))
    if returncode:
        sys.exit(f"Error: 'Picard:CollectRnaSeqMetrics' failed.\nStdout: {stdout}\nStderr: {stderr}\n")

    # parse tool output and put it in qc_metrics.json
    qc_metrics_file = prep_qc_metrics(output_dir, tool_ver)

    # prepare tarball to include output files and qc_metrics.json
    prepare_tarball(args.aligned_seq, qc_metrics_file, output_dir)

if __name__ == "__main__":
    main()

