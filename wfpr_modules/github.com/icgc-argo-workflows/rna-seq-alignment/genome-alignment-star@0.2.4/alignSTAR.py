#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess
import glob
import json
import uuid
from collections import defaultdict

def generate_fastq_from_ubam(ubam, outdir, pair_status, mem=None):

    ubam_base = os.path.basename(ubam)
    ### collect read groups from bam
    rgs = subprocess.run(f'samtools view -H {ubam} | grep -e "^@RG" | cut -f 2 | sed -e "s/^ID://g"', shell=True, check=True, capture_output=True).stdout.decode('utf-8').strip().split('\n')

    ### convert bam to fastq
    for rg in rgs:
        try:
            if pair_status == 'paired':
                subprocess.run(f'samtools collate -u -O {ubam} | samtools view -h -F 256 -r {rg} | samtools fastq -0 /dev/null -1 {outdir}/{rg}_1.fastq.gz -2 {outdir}/{rg}_2.fastq.gz -s /dev/null', shell=True, check=True)
            else:
                subprocess.run(f'samtools collate -u -O {ubam} | samtools view -h -F 256 -r {rg} | samtools fastq -0 /dev/null | gzip > {outdir}/{rg}_1.fastq.gz', shell=True, check=True)
        except Exception as e:
            sys.exit("Error: %s. samtools fastq failed: %s\n" % (e, ubam))


def get_read_group_info(metadata, args):

    read_groups_info = defaultdict(list)
    read_group_ids = []
    read_group_ids_in_bam = []
    pair_status = None
    if 'read_groups' in metadata:
        for rg in metadata['read_groups']:
                
            experiment = metadata['experiment']
            if 'library_strategy' in experiment:
                experimental_strategy = experiment.pop('library_strategy')
                experiment['experimental_strategy'] = experimental_strategy

            ### read group IDs
            read_group_ids.append(rg['submitter_read_group_id'])
            read_group_ids_in_bam.append(rg.get("read_group_id_in_bam") if rg.get("read_group_id_in_bam") else rg.get("submitter_read_group_id"))
            ### pairing status
            _ps = 'paired' if rg.get('is_paired_end') else 'single'
            if pair_status:
                assert pair_status == _ps, 'Error: Read groups in metadata file have inconsistent pairing state'
            pair_status = _ps
            ### mandatory fields
            read_groups_info['SM'].append(metadata['samples'][0]['sampleId'])
            read_groups_info['LB'].append(rg['library_name'])
            read_groups_info['PU'].append(rg['platform_unit'])
            ### optional fields (do not set if 'None')
            if rg.get('insert_size'): read_groups_info['PI'].append(rg.get('insert_size'))
            if rg.get('sample_barcode'): read_groups_info['BC'].append(rg.get('sample_barcode'))
            if experiment.get('sequencing_center'): read_groups_info['CN'].append(experiment.get('sequencing_center'))
            if experiment.get('platform'):
                if experiment.get('platform').upper() in ['CAPILLARY', 'DNBSEQ', 'HELICOS', 'ILLUMINA', 'IONTORRENT', 'LS454', 'ONT', 'PACBIO', 'SOLID']: 
                    read_groups_info['PL'].append(experiment.get('platform'))
                else:
                    sys.stderr.write('Warning: ignored platform: %s - not conform to SAM standard\n' % experiment.get('platform'))
            if experiment.get('platform_model'): read_groups_info['PM'].append(experiment.get('platform_model'))
            if experiment.get('sequencing_date'): read_groups_info['DT'].append(experiment.get('sequencing_date'))
            ### description
            description = '|'.join([
                                        experiment['experimental_strategy'],
                                        metadata['studyId'],
                                        metadata['samples'][0]['specimenId'],
                                        metadata['samples'][0]['donor']['donorId'],
                                        metadata['samples'][0]['specimen']['specimenType'],
                                        metadata['samples'][0]['specimen']['tumourNormalDesignation']
                                    ])
            read_groups_info['DS'].append(description)
    else:
        for rg in args.readgroups:
            read_group_ids.append(rg)
            read_groups_info['SM'] = args.sample
        pair_status = args.pair_status

    return (read_groups_info, read_group_ids, read_group_ids_in_bam, pair_status)


def main():
    """
    Python wrapper for calling STAR aligner on various input configurations,
    assuring that read group information is provided appropriately.

    The relevant input can be provided via a metadata file in json format (via --metadata)
    or as a very basic version using the individual command line flags. The metadata file always 
    takes precedence. Then provided via command line, only readgroup, sample ID and pair-status 
    can be set.
    """

    parser = argparse.ArgumentParser(description='Tool: STAR aligner')
    ### sample metadata 
    parser.add_argument('--sample', dest='sample', type=str,
                        help='Sample name / ID to be processed.', default=None)
    parser.add_argument('--readgroups', dest='readgroups', type=str, nargs='+',
                        help='Readgroup name / ID to be processed.', default=None)
    parser.add_argument('--pair-status', dest='pair_status', type=str, nargs=1,
                        help='Paired-end status of input samples. Same status for all read groups. Choices are: single, paired.', default=None)
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
    parser.add_argument('-t', '--tempdir', dest='tempdir', type=str, default=None,
                        help='Directory for temporary files [./{sample}__STARtmp]')

    args = parser.parse_args()

    ### check whether we are given a metadata file
    ### if a command line parameter is given for a specific setting it 
    ### will overwrite the setting from the metadata file
    if args.metadata:
        with open(args.metadata, 'r') as jsonFile:
            metadata = json.load(jsonFile)
        
        ### make sure we get exactly one sample per metadata file
        assert len(metadata['samples']) == 1, 'Assertion failed: More than one sample provided in metadata file'
        
        #### sample ID
        if 'samples' in metadata:
            args.sample = metadata['samples'][0]['sampleId']
    else:
        metadata = {}

    #### readgroup info
    (read_groups_info, read_group_ids, read_group_ids_in_bam, pair_status) = get_read_group_info(metadata, args)

    ### we make the assumption that all input files have the same format
    ### we derive the input format from the input files
    input_format = None
    for fname in args.input_files:
        if fname.lower().endswith('gz') or fname.lower().endswith('bz2'):
            fname = fname.rsplit('.', 1)[0]
        fmt = fname.lower().rsplit('.', 1)[1]
        if fmt == 'fq':
            fmt = 'fastq'
        assert fmt in ['bam', 'fastq'], 'Error: unknown input file type: %s. Only fastq and bam are allowed\n' % fmt
        if not input_format is None:
            assert fmt == input_format, 'Error: input files have inconsistent file formats\n'
        input_format = fmt

    ### make sure that no input is missing
    assert args.sample is not None, "Error: sample ID has to be provided via command line or metadata file" 
    assert len(read_group_ids) > 0, "Error: at least one readgroup ID has to be provided via command line or metadata file" 
    assert pair_status is not None, "Error: pair_status has to be provided via command line or metadata file" 
    ### make sure the reference files are present
    if not os.path.exists(args.index):
        sys.exit('Error: specified index path %s does not exist or is not accessible!' % args.index)
    if not os.path.exists(args.annotation):
        sys.exit('Error: specified annotation file %s does not exist or is not accessible!' % args.annotation)

    if args.tempdir:
        try:
            if not os.path.exists(args.tempdir):
                os.makedirs(args.tempdir)
        except Exception as e:
            sys.stderr.write('Error: tempdir %s does not exists and could not be created.\n%s' % (args.tempdir, str(e)))
            sys.exit(1)
        outdir = args.tempdir
    else:
        outdir = '.'
    ### handle ubam input
    if input_format == 'bam':
        ### we iterate over all input files of type bam. we make the assumption that the read group ids 
        ### used between bam files do not overlap
        for fname in args.input_files:
            generate_fastq_from_ubam(fname, outdir, pair_status, mem=args.mem)
        fqr1 = glob.glob(os.path.join(outdir, '*_1.fastq.gz'))
        fqr2 = []
        for fq in fqr1:
            _fqr2 = fq[:-10] + '2.fastq.gz'
            if os.path.exists(_fqr2):
                fqr2.append(_fqr2)
            else:
                fqr2.append(None)
        rgs_in_input = [os.path.basename(_)[:-11] for _ in fqr1]

		### map input files to appropriate read group metadata
        rg_idxs = []
        for rg in rgs_in_input:
            if rg in read_group_ids:
                rg_idxs.append(read_group_ids.index(rg))
            elif rg in read_group_ids_in_bam:
                rg_idxs.append(read_group_ids_in_bam.index(rg))
            else:
                sys.stderr.write('Error: readgroup provided in input bam (%s) is not present in the list of readgroups from metadata (%s) or in the list of alternative original readgroups form metadata (%s)\n' % (args.input_files[0], str(read_group_ids), str(read_group_ids_in_bam)))

        ### bring input files in read_group_order
        fqr1_ = []
        fqr2_ = []
        for rg_idx in rg_idxs:
            fqr1_.append(fqr1[rg_idx])
            if pair_status == 'paired':
                assert fqr2[rg_idx] is not None, 'Error: Status of read group %s given as paired, but only single reads could be extracted from %s.' % (read_group_ids[rg_idx], str(args.input_files))
                fqr2_.append(fqr2[rg_idx])
        fqr1 = fqr1_
        fqr2 = fqr2_
                
    ### handle fastq input
    elif input_format == 'fastq':
        ### if paired end, we expect the first half of files to be mate 1 and the second half to be mate 2 in consistent order
        fqr1 = args.input_files
        if pair_status == 'single':
            fqr2 = None
        else:
            assert len(args.input_files) % 2 == 0, 'Error: We expect an even number of fastq files to be provided in paired mode. Currently given: %i' % len(args.input_files)
            fqr1 = args.input_files[:len(args.input_files) // 2]
            fqr2 = args.input_files[len(args.input_files) // 2:]

        ### read group information is derived from fastq file name
        ### - we make the assumption that one (two) fastq files are given in single (paired) mode
        ###   and the read group is derived from the fastq name (by dropping any extensions from the basename)
        rgs_in_input = []
        rgs_in_input2 = []
        for i in range(len(fqr1)):
            if fqr1[i].lower().endswith('gz') or fqr1[i].lower().endswith('bz2'):
                rgs_in_input.append(os.path.basename(fqr1[i]).rsplit('.', 2)[0])
                if pair_status == 'paired':
                    rgs_in_input2.append(os.path.basename(fqr2[i]).rsplit('.', 2)[0])
            else:
                rgs_in_input.append(os.path.basename(fqr1[i]).rsplit('.', 1)[0])
                if pair_status == 'paired':
                    rgs_in_input2.append(os.path.basename(fqr2[i]).rsplit('.', 1)[0])
            if pair_status == 'paired':
                if rgs_in_input[-1][-2:] in ['_1', ':1']: 
                    rgs_in_input[-1] = rgs_in_input[-1][:-2]
                    rgs_in_input2[-1] = rgs_in_input2[-1][:-2]
                assert rgs_in_input[-1] == rgs_in_input2[-1], 'Error: Compression format or read group inconsistent between mates: %s and %s' % (fqr1[i], fqr2[i])

    ### this should not happen
    else:
        sys.exit('Error: The input type needs to be either bam or fastq. Currently given: %s' % input_format)

    ### figure out correct read command
    read_command = 'cat'
    if fqr1[0].lower().endswith('.gz'):
        read_command = 'zcat'
    elif fqr1[0].lower().endswith('.bz2'):
        read_command = 'bzcat'

    ### assembly string of read group IDs
    rgid_string = []
    for i in range(len(read_group_ids)):
        attr = '\t'.join(['%s:%s' % (k, v[i]) for k, v in read_groups_info.items()])
        if len(attr) > 0:
            attr = '"' + attr + '"'
        rgid_string.append('ID:%s %s' % (read_group_ids[i], attr))

    ### assemble STAR command
    cmd = ['STAR',
           '--genomeDir',  args.index,
           '--sjdbGTFfile', args.annotation,
           '--runThreadN', str(args.threads),
           '--sjdbOverhang', str(args.sjdboverhang),
           '--outFileNamePrefix', args.sample + '_',
           '--readFilesIn', ','.join(fqr1), ','.join(fqr2),
           '--outSAMattrRGline', ' , '.join(rgid_string), 
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
    ### optional tempdir
    if args.tempdir:
        star_tmp = os.path.join(args.tempdir, str(uuid.uuid1()))
        while os.path.exists(star_tmp):
            star_tmp = os.path.join(args.tempdir, str(uuid.uuid1()))
        cmd.extend(['--outTmpDir', star_tmp])
        
    ### run STAR
    subprocess.run(' '.join(cmd), shell=True, check=True)

    ### sort output by coordinate
    bam = args.sample + '_Aligned.out.bam'
    subprocess.run(f'samtools sort -o {bam}.sorted {bam} && mv {bam}.sorted {bam}', shell=True, check=True)

    ### reheader bam file to be compliant with ValidateSamFormat
    subprocess.run(f'samtools reheader -P -c \'sed -e "s/\(^@PG.*CL:STAR\).*/\\1 (see @CO for more information)/g"\' {bam} > {bam}.reheadered.bam && mv {bam}.reheadered.bam {bam}', shell=True, check=True)

    ### bundle logs into tarball
    subprocess.run('tar -czf %s.all_logs.supplement.tar.gz *.out align.log' % (args.sample), shell=True, check=True)

if __name__ == "__main__":
    main()
