#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
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
"""

import os
import sys
import json
import argparse
import hashlib
import uuid
import subprocess
import copy
from datetime import date
import re
import tarfile

workflow_full_name = {
    'rna-seq-alignment': 'RNA Seq Alignment'
}

analysis_tools = {
  'star': 'STAR',
  'hisat2': 'HiSAT2'
}

data_type_mapping = {
  #file_type: [dataCategory, dataType, [data_subtypes], [star analysis_tools], [hisat2 analysis_tools]]
  'genome_aln': ['Sequencing Reads', 'Aligned Reads', ['Genome Alignment'], ['STAR'], ['HiSAT2']],
  'transcriptome_aln': ['Sequencing Reads', 'Aligned Reads', ['Transcriptome Alignment'], ['STAR'], ['HiSAT2']],
  'chimeric_aln': ['Sequencing Reads', 'Aligned Reads', ['Chimeric Alignment'], ['STAR'], ['HiSAT2']],
  'splice_junctions': ['Transcriptome Profiling', 'Splice Junctions', [None], ['STAR'], ['HiSAT2']],
  'fastqc': ['Quality Control Metrics', 'Sequencing QC', ['Read Group Metrics'], ['FastQC'], ['FastQC']],
  'collectrnaseqmetrics': ['Quality Control Metrics', 'Aligned Reads QC', ['Alignment Metrics'], ['Picard:CollectRnaSeqMetrics'], ['Picard:CollectRnaSeqMetrics']],
  'duplicates_metrics': ['Quality Control Metrics', 'Aligned Reads QC', ['Duplicates Metrics'], ['biobambam2:bammarkduplicates2'], ['biobambam2:bammarkduplicates2']],
  'supplement': ['Supplement', 'Running Logs', [None], ['STAR'], ['HiSAT2']]
}

def calculate_size(file_path):
    return os.stat(file_path).st_size


def calculate_md5(file_path):
    md5 = hashlib.md5()
    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            md5.update(chunk)
    return md5.hexdigest()

def insert_filename_friendly_rg_id(metadata):
    filename_friendly_rg_ids = set()

    # let's loop it two times, first for the rg id actually doesn't need to convert
    for rg in metadata['read_groups']:
        submitter_read_group_id = rg['submitter_read_group_id']
        filename_friendly_rg_id = "".join([ c if re.match(r"[a-zA-Z0-9\.\-_]", c) else "_" for c in submitter_read_group_id ])

        if filename_friendly_rg_id == submitter_read_group_id:  # no change, safe to use
            rg['filename_friendly_rg_id'] = filename_friendly_rg_id
            filename_friendly_rg_ids.add(filename_friendly_rg_id)

    for rg in metadata['read_groups']:
        submitter_read_group_id = rg['submitter_read_group_id']
        filename_friendly_rg_id = "".join([ c if re.match(r"[a-zA-Z0-9\.\-_]", c) else "_" for c in submitter_read_group_id ])

        if filename_friendly_rg_id == submitter_read_group_id:  # no change, already covered
            continue

        if filename_friendly_rg_id in filename_friendly_rg_ids:  # the converted new friendly ID conflicts with existing one
            for i in range(len(metadata['read_groups'])):
                if not '%s_%s' % (filename_friendly_rg_id, i+1) in filename_friendly_rg_ids:
                    filename_friendly_rg_id += '_%s' % str(i+1)
                    break

        rg['filename_friendly_rg_id'] = filename_friendly_rg_id
        filename_friendly_rg_ids.add(filename_friendly_rg_id)


def get_rg_id_from_ubam_qc(tar, metadata):
    tar_basename = os.path.basename(tar)  # TEST-PR.DO250122.SA610149.D0RE2_1_.6cae87bf9f05cdfaa4a26f2da625f3b2.lane.bam.fastqc.tgz
    md5sum_from_filename = tar_basename.split('.')[-5]
    if not re.match(r'^[a-f0-9]{32}$', md5sum_from_filename):
        sys.exit('Error: ubam naming not expected %s' % tar_basename)

    for rg in metadata.get("read_groups"):
        rg_id_in_bam = rg.get("read_group_id_in_bam") if rg.get("read_group_id_in_bam") else rg.get("submitter_read_group_id")
        seq_file_name = rg.get("file_r1")
        bam_name = seq_file_name if seq_file_name.endswith('.bam') else ''
        md5sum_from_metadata = hashlib.md5(("%s %s" % (bam_name, rg_id_in_bam)).encode('utf-8')).hexdigest()
        if md5sum_from_metadata == md5sum_from_filename:
            return rg.get("filename_friendly_rg_id"), rg.get("submitter_read_group_id")

    # up to this point no match found, then something wrong
    sys.exit('Error: unable to match ubam qc metric tar "%s" to read group id' % tar_basename)


def get_dupmetrics(file_to_upload):
    library = []
    with tarfile.open(file_to_upload, 'r') as tar:
        for member in tar.getmembers():
            if member.name.endswith('.duplicates_metrics.txt'):
                f = tar.extractfile(member)
                cols_name = []
                for r in f:
                    row = r.decode('utf-8')                    
                    if row.startswith('LIBRARY'): 
                        cols_name = row.strip().split('\t')
                        continue
                    if cols_name:
                        if not row.strip(): break
                        metric = {}
                        cols = row.strip().split('\t')
                        for n, c in zip(cols_name, cols):
                            if n == "LIBRARY": metric.update({n: c})
                            elif '.' in c or 'e' in c: metric.update({n: float(c)}) 
                            else: metric.update({n: int(c)})
                        library.append(metric)      
    return library

def get_files_info(file_to_upload, date_str, seq_experiment_analysis_dict, aligner=None):
    file_info = {
        'fileSize': calculate_size(file_to_upload),
        'fileMd5sum': calculate_md5(file_to_upload),
        'fileAccess': 'controlled',
        'info': {}
    }

    experimental_strategy = seq_experiment_analysis_dict['experiment']['experimental_strategy'].lower()
    fname_sample_part = seq_experiment_analysis_dict['samples'][0]['sampleId']

    aligner_or_rgid = aligner.lower() if aligner else None
    submitter_rg_id = None
    if re.match(r'^genome.merged.+?(cram|cram\.crai|bam|bam\.bai)$', file_to_upload): 
      file_type = 'genome_aln'
    elif re.match(r'^transcriptome.merged.+?(cram|cram\.crai|bam|bam\.bai)$', file_to_upload): 
      file_type = 'transcriptome_aln'
    elif re.match(r'^chimeric.merged.+?(cram|cram\.crai|bam|bam\.bai)$', file_to_upload): 
      file_type = 'chimeric_aln'
    elif re.match(r'.+?\.fastqc\.tgz$', file_to_upload):
      file_type = 'fastqc'
      aligner_or_rgid, submitter_rg_id = get_rg_id_from_ubam_qc(file_to_upload, seq_experiment_analysis_dict)
    elif re.match(r'.+?\.collectrnaseqmetrics\.tgz$', file_to_upload):
      file_type = 'collectrnaseqmetrics'
    elif re.match(r'.+?\.duplicates_metrics\.tgz$', file_to_upload):
      file_type = 'duplicates_metrics'
    elif re.match(r'.+?_SJ\.out\.tab$', file_to_upload):
      file_type = 'splice_junctions'
    elif re.match(r'.+?splicesites\.txt$', file_to_upload):
      file_type = 'splice_junctions'
    elif re.match(r'.+?supplement\.tgz$', file_to_upload) or re.match(r'.+?supplement\.tar.gz$', file_to_upload):
      file_type = 'supplement'
    else:
      sys.exit('Error: unknown file type "%s"' % file_to_upload)

    if file_type in ['fastqc', 'collectrnaseqmetrics', 'duplicates_metrics', 'aln_metrics', 'supplement']:
        file_ext = 'tgz'
    elif file_type in ['genome_aln', 'transcriptome_aln', 'chimeric_aln']:
      if file_to_upload.endswith('.bam'):
        file_ext = 'bam'
      elif file_to_upload.endswith('.bam.bai'):
        file_ext = 'bam.bai'
      elif file_to_upload.endswith('.cram'):
        file_ext = 'cram'
      elif file_to_upload.endswith('.cram.crai'):
        file_ext = 'cram.crai'
      else:
        sys.exit('Error: unknown aligned seq extention: %s' % file_to_upload)
    elif file_type in ['splice_junctions']:
      file_ext = 'txt'
    else:
      sys.exit('Error: unknown file type "%s"' % file_type)

    # file naming patterns:
    #   pattern:  <argo_study_id>.<argo_donor_id>.<argo_sample_id>.[rna-seq].<date>.<aligner|rg_id>.<file_type>.<file_ext>
    #   example: TEST-PR.DO250183.SA610229.rna-seq.20200319.star.genome_aln.cram
    new_fname = '.'.join([
                            seq_experiment_analysis_dict['studyId'],
                            seq_experiment_analysis_dict['samples'][0]['donor']['donorId'],
                            fname_sample_part,
                            experimental_strategy,
                            date_str,
                            aligner_or_rgid,
                            file_type,
                            file_ext
                        ])    
    
    file_info['fileName'] = new_fname
    file_info['fileType'] = new_fname.split('.')[-1].upper()

    file_info['info'] = {
        'data_category': data_type_mapping[file_type][0],
        'data_subtypes': data_type_mapping[file_type][2]
    }

    if not aligner:
      file_info['info']['analysis_tools'] = "FastQC"
    elif aligner.lower() == 'star':
      file_info['info']['analysis_tools'] = data_type_mapping[file_type][3]
    elif aligner.lower() == 'hisat2':
      file_info['info']['analysis_tools'] = data_type_mapping[file_type][4]
      

    if new_fname.endswith('.bai') or new_fname.endswith('.crai'):
      file_info['dataType'] = 'Aligned Reads Index'
    else:
      file_info['dataType'] = data_type_mapping[file_type][1]  

    # extract info into payload
    extra_info = {}
    if new_fname.endswith('.tgz'):
      tar = tarfile.open(file_to_upload)
      for member in tar.getmembers():
        if member.name.endswith('qc_metrics.json') or member.name.endswith('.extra_info.json'):
          f = tar.extractfile(member)
          extra_info = json.load(f)
        else:
          if not file_info['info'].get('files_in_tgz'): file_info['info']['files_in_tgz'] = []
          file_info['info']['files_in_tgz'].append(os.path.basename(member.name))

    # retrieve duplicates metrics from the file
    if file_info['info']['data_subtypes'][0] == 'Duplicates Metrics':
      extra_info['metrics'] = {
        'libraries': get_dupmetrics(file_to_upload)
      }
    
    if file_info['info']['data_subtypes'][0] == 'Read Group Metrics':
      extra_info['metrics'].update({'read_group_id': submitter_rg_id})

    if extra_info:
      extra_info.pop('tool', None)
      file_info['info'].update(extra_info)
      
    new_dir = 'out'
    try:
      os.mkdir(new_dir)
    except FileExistsError:
      pass

    dst = os.path.join(os.getcwd(), new_dir, new_fname)
    os.symlink(os.path.abspath(file_to_upload), dst)

    return file_info

def get_sample_info(sample_list):
    samples = copy.deepcopy(sample_list)
    for sample in samples:
        for item in ['info', 'sampleId', 'specimenId', 'donorId', 'studyId']:
            sample.pop(item, None)
            sample['specimen'].pop(item, None)
            sample['donor'].pop(item, None)

    return samples


def main():
    """
    Python implementation of tool: payload-gen-rna-alignment
    """

    parser = argparse.ArgumentParser(description='Tool: payload-gen-rna-alignment')
    parser.add_argument("-f", "--files_to_upload", dest="files_to_upload", type=str, required=True,
                        nargs="+", help="Files to upload")
    parser.add_argument("-a", "--seq_experiment_analysis", dest="seq_experiment_analysis", required=True,
                        help="Input analysis for sequencing experiment", type=str)
    parser.add_argument("-t", "--analysis_type", dest="analysis_type", required=True, help="Specify analysis_type")
    parser.add_argument("-l", "--aligner", dest="aligner", default=None, help="Provide RNA-Seq aligner if files_to_upload are generated from alignment results. Default=None")
    parser.add_argument("-g", "--genome_annotation", dest="genome_annotation", default="GENCODE v38", help="RNA-Seq alignment genome annotation")
    parser.add_argument("-b", "--genome_build", dest="genome_build", default="GRCh38_hla_decoy_ebv", help="RNA-Seq alignment genome build")
    parser.add_argument("-w", "--wf_name", dest="wf_name", required=True, help="Workflow name")
    parser.add_argument("-v", "--wf_version", dest="wf_version", required=True, help="Workflow version")
    parser.add_argument("-r", "--wf_run", dest="wf_run", required=True, help="Workflow run ID")
    parser.add_argument("-s", "--wf_session", dest="wf_session", required=True, help="Workflow session ID")
    args = parser.parse_args()

    with open(args.seq_experiment_analysis, 'r') as f:
        seq_experiment_analysis_dict = json.load(f)

    payload = {
        'analysisType': {
            'name': args.analysis_type
        },
        'studyId': seq_experiment_analysis_dict.get('studyId'),
        'workflow': {
            'workflow_name': workflow_full_name.get(args.wf_name, args.wf_name),
            'workflow_version': args.wf_version,
            'genome_build': args.genome_build,
            'genome_annotation': args.genome_annotation,
            'run_id': args.wf_run,
            'session_id': args.wf_session,
            'inputs': [
                {
                    'analysis_type': 'sequencing_experiment',
                    'input_analysis_id': seq_experiment_analysis_dict.get('analysisId')
                }
            ]
        },
        'files': [],
        'samples': get_sample_info(seq_experiment_analysis_dict.get('samples')),
        'experiment': seq_experiment_analysis_dict.get('experiment')
    }

    if "sequencing_alignment" in args.analysis_type:
      payload['read_group_count'] = seq_experiment_analysis_dict.get('read_group_count')
      payload['read_groups'] = copy.deepcopy(seq_experiment_analysis_dict.get('read_groups'))

    # pass `info` dict from seq_experiment payload to new payload
    if 'info' in seq_experiment_analysis_dict and isinstance(seq_experiment_analysis_dict['info'], dict):
      payload['info'] = seq_experiment_analysis_dict['info']

    if 'library_strategy' in payload['experiment']:
      experimental_strategy = payload['experiment'].pop('library_strategy')
      payload['experiment']['experimental_strategy'] = experimental_strategy

    insert_filename_friendly_rg_id(seq_experiment_analysis_dict)

    # get file of the payload
    date_str = date.today().strftime("%Y%m%d")
    for f in args.files_to_upload:
      file_info = get_files_info(f, date_str, seq_experiment_analysis_dict, args.aligner)
      payload['files'].append(file_info)
    

    with open("%s.%s.payload.json" % (str(uuid.uuid4()), args.analysis_type), 'w') as f:
      f.write(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

