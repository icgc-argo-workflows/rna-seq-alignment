#!/usr/bin/env python3

"""
 Copyright (c) 2019-2021, Ontario Institute for Cancer Research (OICR).

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as published
 by the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with this program. If not, see <https://www.gnu.org/licenses/>.

 Authors:
   Linda Xiang <linda.xiang@oicr.on.ca>
   Junjun Zhang <junjun.zhang@oicr.on.ca>
 """


import sys
import uuid
import json
import csv
import textwrap
from argparse import ArgumentParser


TSV_FIELDS = {
    'experiment': [
        'type', 'program_id', 'submitter_sequencing_experiment_id', 'submitter_donor_id', 'gender',
        'submitter_specimen_id', 'tumour_normal_designation', 'specimen_type', 'specimen_tissue_source', 'submitter_sample_id',
        'sample_type', 'submitter_matched_normal_sample_id', 'sequencing_center', 'platform', 'platform_model',
        'experimental_strategy', 'sequencing_date', 'read_group_count'
    ],
    'read_group': [
        'type', 'submitter_read_group_id', 'read_group_id_in_bam', 'submitter_sequencing_experiment_id', 'platform_unit',
        'is_paired_end', 'file_r1', 'file_r2', 'read_length_r1', 'read_length_r2', 'insert_size', 'sample_barcode', 'library_name'
    ],
    'file': [
        'type', 'name', 'size', 'md5sum', 'path', 'format'
    ]
}


def empty_str_to_null(metadata):
    for k in metadata:
        if k in ['read_groups', 'files']:
            for i in range(len(metadata[k])):
                empty_str_to_null(metadata[k][i])
        if isinstance(metadata[k], str) and metadata[k] in ["", "_NULL_"]:
            metadata[k] = None


def tsv_confomity_check(ftype, tsv):
    expected_fields = TSV_FIELDS[ftype]

    header_processed = False
    with open(tsv, 'r') as t:
        uniq_row = {}
        for l in t:
            l = l.rstrip('\n').rstrip('\r')  # remove trailing newline, remove windows `\r` (just in case)
            if not header_processed:  # it's header
                fields = l.split('\t')
                if len(fields) != len(set(fields)):
                    sys.exit("Error found: Field duplicated in input TSV: %s, offending header: %s\n" % (tsv, l))

                missed_fields = set(expected_fields) - set(fields)
                if missed_fields:  # missing fields
                    sys.exit("Error found: Field missing in input TSV: %s, offending header: %s. Missed field(s): %s\n" % \
                        (tsv, l, ', '.join(missed_fields)))

                unexpected_fields = set(fields) - set(expected_fields)
                if unexpected_fields:  # unexpected fields
                    sys.exit("Error found: Unexpected field in input TSV: %s, offending header: %s. Unexpected field(s): %s\n" % \
                        (tsv, l, ', '.join(unexpected_fields)))

                header_processed = True

            else:  # it's data row
                # at this point we only check whether number of values matches number of expected fields and uniqueness check,
                # later steps will perform more sophisticated content check
                values = l.split('\t')
                if len(expected_fields) != len(values):
                    sys.exit("Error found: number of fields: %s does not match expected: %s, offending data row: %s\n" % \
                        (len(values), len(expected_fields), l))

                if l in uniq_row:
                    sys.exit("Error found: data row repeated in file: %s, offending data row: %s\n" % (tsv, l))
                else:
                    uniq_row[l] = True


def load_all_tsvs(exp_tsv, rg_tsv, file_tsv):
    metadata_dict = {}
    with open(exp_tsv, 'r') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
        if len(rows) != 1:
            sys.exit("Error found: experiment TSV expects exactly one data row, offending file: %s has %s row(s)\n" % \
                (exp_tsv, len(rows)))
        rows[0]['read_group_count'] = int(rows[0]['read_group_count'])
        metadata_dict.update(rows[0])

    with open(rg_tsv, 'r') as f:
        metadata_dict['read_groups'] = []
        for rg in csv.DictReader(f, delimiter='\t'):
            if rg['is_paired_end'].lower() == 'true':
                rg['is_paired_end'] = True
            elif rg['is_paired_end'].lower() == 'false':
                rg['is_paired_end'] = False
            else:
                rg['is_paired_end'] = None

            for field in ('read_length_r1', 'read_length_r2', 'insert_size'):
                if rg[field]:
                    rg[field] = int(rg[field])
                else:
                    rg[field] = None

            metadata_dict['read_groups'].append(rg)

        if len(metadata_dict['read_groups']) == 0:
            sys.exit("Error found: read group TSV does not contain any read group information\n")

    with open(file_tsv, 'r') as f:
        metadata_dict['files'] = []
        for f in csv.DictReader(f, delimiter='\t'):
            if f['size']:
                f['size'] = int(f['size'])
            else:
                f['size'] = None

            metadata_dict['files'].append(f)

        if len(metadata_dict['files']) == 0:
            sys.exit("Error found: file TSV does not contain any file information\n")

    return metadata_dict


def validate_args(args):
    if args.metadata_json and \
            not (args.experiment_info_tsv or args.read_group_info_tsv or args.file_info_tsv):
        return True
    elif not args.metadata_json and \
            (args.experiment_info_tsv and args.read_group_info_tsv and args.file_info_tsv):
        return True
    else:
        sys.exit(textwrap.dedent(
            """
            Usage:
                When '-m' is provided, no other arguments can be used
                When '-m' is not provided, please provide all of these arguments: -x, -r and -f
            """
        ))


def main(metadata, extra_info=dict()):
    empty_str_to_null(metadata)

    payload = {
        'analysisType': {
            'name': 'sequencing_experiment'
        },
        'studyId': metadata.get('program_id'),
        'experiment': {
            'submitter_sequencing_experiment_id': metadata.get('submitter_sequencing_experiment_id'),
            'sequencing_center': metadata.get('sequencing_center'),
            'platform': metadata.get('platform'),
            'platform_model': metadata.get('platform_model'),
            'experimental_strategy': metadata.get('experimental_strategy'),
            'sequencing_date': metadata.get('sequencing_date')
        },
        'read_group_count': metadata.get('read_group_count'),
        'read_groups': [],
        'samples': [],
        'files': []
    }

    # get sample of the payload
    sample = {
        'submitterSampleId': metadata.get('submitter_sample_id'),
        'matchedNormalSubmitterSampleId': metadata.get('submitter_matched_normal_sample_id'),
        'sampleType': metadata.get('sample_type'),
        'specimen': {
            'submitterSpecimenId': metadata.get('submitter_specimen_id'),
            'tumourNormalDesignation': metadata.get('tumour_normal_designation'),
            'specimenTissueSource': metadata.get('specimen_tissue_source'),
            'specimenType': metadata.get('specimen_type')
        },
        'donor': {
            'submitterDonorId': metadata.get('submitter_donor_id'),
            'gender': metadata.get('gender')
        }
    }

    if extra_info:
        if extra_info['sample'].get(sample['submitterSampleId']):
            sample['sampleId'] = extra_info['sample'][sample['submitterSampleId']]
        else:
            sys.exit(f"Provided extra_info_tsv misses mapping for submitter sample ID: {sample['submitterSampleId']}")

        if extra_info['specimen'].get(sample['specimen']['submitterSpecimenId']):
            sample['specimenId'] = extra_info['specimen'][sample['specimen']['submitterSpecimenId']]
            sample['specimen']['specimenId'] = sample["specimenId"]
        else:
            sys.exit(f"Provided extra_info_tsv misses mapping for submitter specimen ID: {sample['specimen']['submitterSpecimenId']}")

        if extra_info['donor'].get(sample['donor']['submitterDonorId']):
            sample['donor']['donorId'] = extra_info['donor'][sample['donor']['submitterDonorId']]
            sample['specimen']['donorId'] = sample['donor']['donorId']
        else:
            sys.exit(f"Provided extra_info_tsv misses mapping for submitter donor ID: {sample['donor']['submitterDonorId']}")

    payload['samples'].append(sample)

    # get file of the payload
    for input_file in metadata.get("files"):
        payload['files'].append(
            {
                'fileName': input_file.get('name'),
                'fileSize': input_file.get('size'),
                'fileMd5sum': input_file.get('md5sum'),
                'fileType': input_file.get('format'),
                'fileAccess': 'controlled',
                'dataType': 'Submitted Reads',
                'info': {
                    'data_category': 'Sequencing Reads'
                }
            }
        )

    for rg in metadata.get("read_groups"):
        rg.pop('type')  # remove 'type' field
        rg.pop('submitter_sequencing_experiment_id')  # remove 'submitter_sequencing_experiment_id' field
        payload['read_groups'].append(rg)

    with open("%s.sequencing_experiment.payload.json" % str(uuid.uuid4()), 'w') as f:
        f.write(json.dumps(payload, indent=2))


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-m", "--metadata-json",
                        help="json file containing experiment, read_group and file information submitted from user")
    parser.add_argument("-x", "--experiment-info-tsv",
                        help="tsv file containing experiment information submitted from user")
    parser.add_argument("-r", "--read-group-info-tsv",
                        help="tsv file containing read_group information submitted from user")
    parser.add_argument("-f", "--file-info-tsv",
                        help="tsv file containing file information submitted from user")
    parser.add_argument("-e", "--extra-info-tsv",
                        help="tsv file containing file information submitted from user")
    args = parser.parse_args()

    validate_args(args)

    if args.metadata_json:
        with open(args.metadata_json, 'r') as f:
            metadata = json.load(f)
    else:
        # fistly TSV format conformity check, if not well-formed no point to continue
        tsv_confomity_check('experiment', args.experiment_info_tsv)
        tsv_confomity_check('read_group', args.read_group_info_tsv)
        tsv_confomity_check('file', args.file_info_tsv)

        # all TSV are well-formed, let's load them
        metadata = load_all_tsvs(
                            args.experiment_info_tsv,
                            args.read_group_info_tsv,
                            args.file_info_tsv
                        )

        # all TSV are well-formed, let's load them
        metadata = load_all_tsvs(args.experiment_info_tsv, args.read_group_info_tsv, args.file_info_tsv)

    extra_info = dict()
    if args.extra_info_tsv:
        with open(args.extra_info_tsv, 'r') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                type = row['type']
                submitter_id = row['submitter_id']
                uniform_id = row['uniform_id']
                if type in extra_info:
                    sys.exit(f"Values in 'type' field duplicated. Offending value: {type}, in file: {args.extra_info_tsv}")
                else:
                    extra_info[type] = dict()

                if submitter_id in extra_info[type]:
                    sys.exit(f"Values in 'submitter_id' field duplicated. Offending value: {submitter_id}, for type: {type}, in file: {args.extra_info_tsv}" )
                else:
                    extra_info[type][submitter_id] = uniform_id

        if 'donor' not in extra_info or 'specimen' not in extra_info or 'sample' not in extra_info:
            sys.exit(f"Provided extra_info_tsv file '{args.extra_info_tsv}' is required to have ID mappings for 'donor', 'specimen' and 'sample'")

    main(metadata, extra_info)
