#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess
import glob
import json

def replace_bam_header(bam, header, out, mem=None):

    jvm_Xmx = "-Xmx%iM" % int(mem) if mem else ""
    try:
        cmd = [
                'java', jvm_Xmx, '-Dpicard.useLegacyParser=false', '-jar', '/tools/picard.jar', 'ReplaceSamHeader', '-I', bam,
                '--HEADER', header, '--OUTPUT', out
              ]

        subprocess.run(cmd, check=True)
    except Exception as e:
        sys.exit("Error: %s. ReplaceSamHeader failed: %s\n" % (e, bam))


def generate_fastq_from_ubam(ubam, outdir, mem=None):

    ### convert bam to fastq
    ### In SamToFastq the default setting for --INCLUDE_NON_PRIMARY_ALIGNMENTS is false. So we will exclude secondary hits.
    ### Further, we will output one fastq file per given read group, but we will only process the read group specified in 
    ### the given metadata. 
    jvm_Xmx = "-Xmx%iM" % int(mem) if mem else ""
    try:
        cmd = [
                'java', jvm_Xmx, '-Dpicard.useLegacyParser=false', '-jar', '/tools/picard.jar', 'SamToFastq', '-I', ubam,
                '--OUTPUT_PER_RG', 'true', '--RG_TAG', 'ID', '--COMPRESS_OUTPUTS_PER_RG', 'true', '--OUTPUT_DIR', outdir
              ]
        subprocess.run(cmd, check=True)
    except Exception as e:
        sys.exit("Error: %s. SamToFastq failed: %s\n" % (e, ubam))

def main():
    """
    Python wrapper for calling STAR aligner on various input configurations,
    assuring that read group information is provided appropriately.

    The relevant input can be provided via a metadata file in json format (via --metadata)
    or using the individual command line flags. A combination of both can be used. In this case 
    individual settings made via command line flag overwrite the values provided in the
    metadata file.
    """

    parser = argparse.ArgumentParser(description='Tool: STAR aligner')
    ### sample metadata
    parser.add_argument('--sample', dest='sample', type=str,
                        help='Sample name / ID to be processed.', default=None)
    parser.add_argument('--readgroup', dest='readgroup', type=str, nargs='+',
                        help='Readgroup name / ID to be processed.', default=None)
    parser.add_argument('--readgroup-orig', dest='readgroup_orig', type=str, nargs='+',
                        help='Readgroup name / ID as present in original data.', default=None)
    parser.add_argument('--pair-status', dest='pair_status', type=str, nargs='+',
                        help='Paired-end status of input samples per read_group. Choices are: single, paired.', default=None)
    parser.add_argument('--input-format', dest='input_format', type=str,
                        help='Format of read in put: fastq or bam', default=None)
    parser.add_argument('--metadata', dest='metadata', type=str,
                        help='Metadata file containing descriptive sample information', default=None)
    ### processing metadata
    parser.add_argument('--input-files', dest='input_files', type=str, nargs='+',
                        help='Input read files in fastq or bam format. For paired end fastq, first mate comes first.', required=True)
    parser.add_argument('--index', dest='index', type=str,
                        help='Path to STAR genome index.', required=True)
    parser.add_argument('--annotation', dest='annotation', type=str,
                        help='Gene annotation (preferably in GTF format).', required=True)
    parser.add_argument('--sjdbOverhang', dest='sjdboverhang', type=int,
                        help='Overhang for splice junction DB. [100]', default=100, required=False)
    parser.add_argument('--threads', dest='threads', type=int,
                        help='Number of threads. [1]', default=1, required=False)
    parser.add_argument('--mem', dest='mem', help="Maximal allocated memory in MB", type=float, default=None)

    args = parser.parse_args()

    ### check whether optional read group substitution is correctly used
    if args.readgroup is not None and args.readgroup_orig is not None:
        assert len(args.readgroup) == len(args.readgroup_orig), 'Error: Number and order of RGIDs in --readgroup and --readgroup-orig need to match'

    ### check whether we are given a metadata file
    ### if a command line parameter is given for a specific setting it 
    ### will overwrite the setting from the metadata file
    if args.metadata:
        with open(args.metadata, 'r') as jsonFile:
            metadata = json.load(jsonFile)
        
        ### make sure we get exactly one sample per metadata file
        assert len(metadata['samples']) == 1, 'Assertion failed: More than one sample provided in metadata file'
        
        ### sample ID
        if args.sample is None and 'samples' in metadata:
            args.sample = metadata['samples'][0]['sampleId']
        ### readgroup info
        ### - we make the assumption that the input files are given for exactly one read group. 
        ###   That is, for ubam, the file must have only a single read group
        ### - in case of fastq input, we use the file name without extensions as read group ID
        if 'read_groups' in metadata:
            if args.readgroup is None:
                args.readgroup = []
                args.readgroup_orig = []
                for rg in sorted(metadata['read_groups']):
                    args.readgroup.append(rg['submitter_read_group_id'])
                    args.readgroup_orig.append(rg['read_group_id_in_bam'])
            ### single/paired
            if args.pair_status is None:
                args.pair_status = []
                for rg in sorted(metadata['read_groups']):
                    if rg['is_paired_end']:
                        args.pair_status.append('paired')
                    else:
                        args.pair_status.append('single')

    ### we make the assumption that all input files have the same format
    ### we derive the input format from the input files
    if args.input_format is None:
        for fname in args.input_files:
            if fname.lower().endswith('gz') or fname.lower().endswith('bz2'):
                fname = fname.rsplit('.', 1)[0]
            fmt = fname.lower().rsplit('.', 1)[1]
            assert fmt in ['bam', 'fastq'], 'Error: unknown input file type: %s. Only fastq and bam are allowed\n' % fmt
            if not args.input_format is None:
                assert fmt == args.input_format, 'Error: input files have inconsistent file formats\n'
            args.input_format = fmt

    ### make sure that no input is missing
    assert args.sample is not None, "Error: sample ID has to be provided via command line or metadata file" 
    assert args.readgroup is not None, "Error: readgroup ID has to be provided via command line or metadata file" 
    assert args.pair_status is not None, "Error: pair_status has to be provided via command line or metadata file" 
    ### make sure the reference files are present
    if not os.path.exists(args.index):
        sys.exit('Error: specified index path %s does not exist or is not accessible!' % args.index)
    if not os.path.exists(args.annotation):
        sys.exit('Error: specified annotation file %s does not exist or is not accessible!' % args.annotation)

    ### handle ubam input
    outdir = '.'
    if args.input_format == 'bam':
        if len(args.input_files) != 1:
            sys.exit('Error: number of input files %s needs to be exactly 1!.' % str(args.input_files))
        generate_fastq_from_ubam(args.input_files[0], outdir, mem=args.mem)
        fqr1 = glob.glob(os.path.join(outdir, '*_1.fastq.gz'))
        fqr2 = []
        for fq in fqr1:
            _fqr2 = fq[:-10] + '2.fastq.gz'
            if os.path.exists(_fqr2):
                fqr2.append(_fqr2)
            else:
                fqr2.append(None)
        rg_input = [os.path.basename(_)[:-11] for _ in fqr1]

        ### we can not have more than one read group per given ubam file
        assert len(rg_input) == 1, 'Error: the number of read groups present in given bam file %s is different from 1' % args.input_files[0]
        rg_input = rg_input[0]
        
        ### make sure that the read group present in the ubam also is given in the metadata
        if rg_input in args.readgroup:
            rg_idx = args.readgroup.index(rg_input)
        elif rg_input in args.readgroup_orig:
            rg_idx = args.readgroup_orig.index(rg_input) 
        else:
            sys.stderr.write('Error: readgroup provided in input bam (%s) is not present in the list of readgroups from metadata (%s) or in the list of alternative original readgroups form metadata (%s)\n' % (args.input_files[0], str(args.readgroup), str(args.readgroup_orig)))
        args.readgroup = args.readgroup[rg_idx]
        args.pair_status = args.pair_status[rg_idx]
        fqr1 = fqr1[rg_idx]
        if args.pair_status == 'paired':
            fqr2 = fqr2[rg_idx]
            assert fqr2 is not None, 'Error: Status of read group %s given as paired, but only single reads could be extracted from %s.' % (args.readgroup, args.input_files[0])
                

    ### handle fastq input
    elif args.input_format == 'fastq':
        assert len(args.input_files) in [1, 2], 'Error: We expect 1 fastq file in single and 2 fastq files in paired mode. Currently given: %i' % len(args.input_files)
        if len(args.pair_status) > 1:
            sys.stderr.write('Warning: Currently %i values are given for --pair-status. Only the first one will be considered.\n' % len(args.pair_status))
        args.pair_status = args.pair_status[0]
        if args.readgroup is not None:
            if len(args.readgroup) > 1:
                sys.stderr.write('Warning: Currently %i values are given for --readgroup. Only the first one will be considered.\n' % len(args.readgroup))
            args.readgroup = args.readgroup[0]
        if args.readgroup_orig is not None:
            if len(args.readgroup_orig):
                sys.stderr.write('Warning: Currently %i values are given for --readgroup-orig. Only the first one will be considered.\n' % len(args.readgroup_orig))
            args.readgroup_orig = args.readgroup_orig[0]
        fqr1 = args.input_files[0]
        fqr2 = None
        if args.pair_status == 'paired': 
            fqr2 = args.input_files[1]
            if len(args.input_files) != 2:
                sys.exit('Error: Paired-end status was given as %s. But files provided were: %s' % (args.pair_status, str(args.input_files)))
        elif args.pair_status == 'single' and len(args.input_files != 1):
            sys.exit('Error: Paired-end status was given as %s. But files provided were: %s' % (args.pair_status, str(args.input_files)))

        ### read group information is derived from fastq file name
        ### - we make the assumption that one (two) fastq files are given in single (paired) mode
        ###   and the read group is derived from the fastq name (by dropping any extensions from the basename)
        if fqr1.lower().endswith('gz') or fqr1.lower().endswith('bz2'):
            rg_input = os.path.basename(fqr1).rsplit('.', 2)[0]
        else:
            rg_input = os.path.basename(fqr1).rsplit('.', 1)[0]
        if args.pair_status == 'paired' and rg_input[-2:] in ['_1', ':1']: 
            rg_input = rg_input[:-2]

        ### if an original read group ID is provided, check whether we can use it to map to the new submitter RGID
        if args.readgroup_orig is not None and args.readgroup is not None:
            assert rg_input in args.readgroup_orig, 'Error: read group ID provided via readgroup_orig (%s) is inconsistent with RGID derived from fastq file name (%s).' % (str(args.readgroup_orig), rg_input)
            rg_idx = args.readgroup_orig.index(rg_input)
            args.readgroup = args.readgroup[rg_idx]
        else:
            args.readgroup = rg_input

    ### this should not happen
    else:
        sys.exit('Error: The input type needs to be either bam or fastq. Currently given: %s' % args.input_format)

    ### figure out correct read command
    read_command = 'cat'
    if fqr1.lower().endswith('.gz'):
        read_command = 'zcat'
    elif fqr1.lower().endswith('.bz2'):
        read_command = 'bzcat'

    ### assemble STAR command
    cmd = ['STAR',
           '--genomeDir',  args.index,
           '--sjdbGTFfile', args.annotation,
           '--runThreadN', str(args.threads),
           '--sjdbOverhang', str(args.sjdboverhang),
           '--outFileNamePrefix', args.sample + '_' + args.readgroup + '_',
           '--readFilesIn', fqr1, fqr2,
           '--outSAMattrRGline', 'ID:%s\tSM:%s' % (args.readgroup, args.sample),
           '--readFilesCommand', read_command,
           '--twopassMode Basic',
           '--outFilterMultimapScoreRange', '1',
           '--outFilterMultimapNmax', '20', 
           '--outFilterMismatchNmax', '10',
           '--alignIntronMax', '500000', 
           '--alignMatesGapMax', '1000000',
           '--sjdbScore', '2',
           '--alignSJDBoverhangMin', '1',
           '--genomeLoad NoSharedMemory',
           '--outFilterMatchNminOverLread', '0.33',
           '--outFilterScoreMinOverLread', '0.33',
           '--outSAMstrandField intronMotif',
           '--outSAMmode Full',
           '--outSAMattributes NH HI NM MD AS XS',
           '--outSAMunmapped Within',
           '--limitSjdbInsertNsj', '2000000',
           '--outSAMtype BAM Unsorted',
           '--outSAMheaderHD', '@HD VN:1.4',
           '--outSAMmultNmax', '1',
          ]
    subprocess.run(' '.join(cmd), shell=True, check=True)

    ### sort output by coordinate
    bam = args.sample + '_' + args.readgroup + '_Aligned.out.bam'
    subprocess.run(f'samtools sort -o {bam}.sorted {bam} && mv {bam}.sorted {bam}', shell=True, check=True)

    ### bundle logs into tarball
    subprocess.run('tar -czf %s_%s.all_logs.supplement.tar.gz *.out align.log' % (args.sample, args.readgroup), shell=True, check=True)

if __name__ == "__main__":
    main()
