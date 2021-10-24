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
    parser.add_argument('--readgroup', dest='readgroup', type=str,
                        help='Readgroup name / ID to be processed.', default=None)
    parser.add_argument('--pair-status', dest='pair_status', type=str,
                        help='Paired-end status of input samples. Choices are: single, paired.', default=None)
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

    ### check whether we are given a metadata file
    ### if a command line parameter is given for a specific setting it 
    ### will overwrite the setting from the metadata file
    if args.metadata:
        with open(args.metadata, 'r') as jsonFile:
            metadata = json.load(jsonFile)
        
        ### make sure we get exactly one read group from one sample
        assert len(metadata['read_groups']) == 1, 'Assertion failed: More than one read group provided in metadata file'
        assert len(metadata['samples']) == 1, 'Assertion failed: More than one sample provided in metadata file'
        
        ### sample ID
        if args.sample is None and 'samples' in metadata:
            args.sample = metadata['samples'][0]['sampleId']
        ### readgroup info
        if 'read_groups' in metadata:
            if args.readgroup is None:
                args.readgroup = metadata['read_groups'][0]['submitter_read_group_id']
            ### single/paired
            if args.pair_status is None:
                args.pair_status = 'paired' if metadata['read_groups'][0]['is_paired_end'] else 'single'
        ### we make the assumption that all input files have the same format
        if args.input_format is None and 'files' in metadata:
            args.input_format = metadata['files'][0]['fileType'].lower()

    ### make sure that no input is missing
    assert args.sample is not None, "Error: sample ID has to be provided via command line or metadata file" 
    assert args.readgroup is not None, "Error: readgroup ID has to be provided via command line or metadata file" 
    assert args.pair_status is not None, "Error: pair_status has to be provided via command line or metadata file" 
    assert args.input_format is not None, "Error: input_format has to be provided via command line or metadata file" 
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
        if args.pair_status == 'paired':
            for fq in fqr1:
                fqr2.append(fq[:-10] + '2.fastq.gz')
                assert os.path.exists(fqr2[-1])
        rgs_bam = [os.path.basename(_)[:-11] for _ in fqr1]
        ### make sure that the read group requested in the parameters is present in the bam
        assert args.readgroup in rgs_bam, 'Error: requested read group (%s) not present in given bam (%s)' % (args.readgroup, args.input_files[0])

    ### handle fastq input
    elif args.input_format == 'fastq':
        fqr1 = [args.input_files[0]]
        fqr2 = []
        if args.pair_status == 'paired': 
            fqr2 = [args.input_files[1]]
            if len(args.input_files) != 2:
                sys.exit('Error: Paired-end status was given as %s. But files provided were: %s' % (args.pair_status, str(args.input_files)))
        elif args.pair_status == 'single' and len(args.input_files != 1):
            sys.exit('Error: Paired-end status was given as %s. But files provided were: %s' % (args.pair_status, str(args.input_files)))
    ### this should not happen
    else:
        sys.exit('Error: The input type needs to be either bam or fastq. Currently given: %s' % args.input_format)

    ### figure out correct read command
    read_command = 'cat'
    for fq in fqr1:
        if fq.lower().endswith('.gz'):
            read_command = 'zcat'
            break
        if fq.lower().endswith('.bz2'):
            read_command = 'bzcat'

    ### assemble STAR command
    cmd = ['STAR',
           '--genomeDir',  args.index,
           '--sjdbGTFfile', args.annotation,
           '--runThreadN', str(args.threads),
           '--sjdbOverhang', str(args.sjdboverhang),
           '--outFileNamePrefix', args.sample + '_' + args.readgroup + '_',
           '--readFilesIn', ','.join(fqr1), ','.join(fqr2),
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

    ### replace original read group line from ubam
    #if args.input_format == 'bam':
    #    bam = args.readgroup + '_Aligned.out.bam'
    #    ### get current header and drop RG info
    #    subprocess.run(f'samtools view -H {bam} --no-PG | grep -v "@RG" > new_header.sam', shell=True, check=True)
    #    ### append old read group info 
    #    for fq in fqr1:
    #        subprocess.run('samtools view -H %s --no-PG | grep -e @RG | grep %s >> new_header.sam' % (args.input_files[0], os.path.basename(fq)[:-11]), shell=True, check=True)
    #    ### replace header
    #    bam_rg = args.sample + '_Aligned.out.rg.bam'
    #    replace_bam_header(bam, 'new_header.sam', bam_rg, mem=args.mem) 
    #    ### clean up
    #    subprocess.run('mv %s %s && rm new_header.sam' % (bam_rg, bam), shell=True, check=True)

if __name__ == "__main__":
    main()
