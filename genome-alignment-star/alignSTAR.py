#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess
import glob

def generate_fastq_from_ubam(ubam, outdir, mem=None):

    jvm_Xmx = "-Xmx%iM" % int(mem) if mem else ""
    try:
        cmd = [
                'java', jvm_Xmx, '-Dpicard.useLegacyParser=false', '-jar', '/tools/picard.jar', 'SamToFastq', '-I', ubam,
                '--OUTPUT_PER_RG', 'true', '--COMPRESS_OUTPUTS_PER_RG', 'true', '--OUTPUT_DIR', outdir
              ]

        subprocess.run(cmd, check=True)
    except Exception as e:
        sys.exit("Error: %s. SamToFastq failed: %s\n" % (e, ubam))

def main():
    """
    Python wrapper for calling STAR aligner on various input configurations,
    assuring that read group information is provided appropriately.
    """

    parser = argparse.ArgumentParser(description='Tool: STAR aligner')
    parser.add_argument('--sample', dest='sample', type=str,
                        help='Input sample name / ID.', required=True)
    parser.add_argument('--index', dest='index', type=str,
                        help='Path to STAR genome index.', required=True)
    parser.add_argument('--annotation', dest='annotation', type=str,
                        help='Gene annotation (preferably in GTF format).', required=True)
    parser.add_argument('--sjdbOverhang', dest='sjdboverhang', type=int,
                        help='Overhang for splice junction DB. [100]', default=100, required=False)
    parser.add_argument('--threads', dest='threads', type=int,
                        help='Number of threads. [1]', default=1, required=False)
    parser.add_argument('--pair-status', dest='pair_status', type=str,
                        help='Paired-end status of input samples. Choices are: single, paired.', required=True)
    parser.add_argument('--input-fastq-r1', dest='fqr1', type=str, nargs='*',
                        help='Input read files in fastq format.First mate(s) or single-end read sample(s).', default=[])
    parser.add_argument('--input-fastq-r2', dest='fqr2', type=str, nargs='*',
                        help='Input read files in fastq format.First mate(s) or single-end read sample(s).', default=[])
    parser.add_argument('--input-ubam', dest='ubam', type=str,
                        help='Input read file in unaligned BAM format.', default='', required=False)
    parser.add_argument('--mem', dest='mem', help="Maximal allocated memory in MB", type=float, default=None)

    args = parser.parse_args()

    #if not os.path.isfile(args.input_file):
    #    sys.exit('Error: specified input file %s does not exist or is not accessible!' % args.input_file)

    #if not os.path.isdir(args.output_dir):
    #    sys.exit('Error: specified output dir %s does not exist or is not accessible!' % args.output_dir)

    ### handle ubam input input
    outdir = '.'
    if len(args.ubam) > 0:
        generate_fastq_from_ubam(args.ubam, outdir, mem=args.mem)
        fqr1 = glob.glob(os.path.join(outdir, '*_1.fastq.gz'))
        fqr2 = []
        if args.pair_status == 'paired':
            for fq in fqr1:
                fqr2.append(fq[:-10] + '2.fastq.gz')
                assert os.path.exists(fqr2[-1])
        rgs = [_[:-10] for _ in fqr1]
    ### handle fastq input
    elif len(args.fqr1) > 0:
        fqr1 = args.fqr1
        fqr2 = args.fqr2
        if args.pair_status == 'paired':
            assert len(args.fqr1) == len(args.fqr2)
        else:
            if len(args.fqr2) > 0:
                sys.stderr.write('Warning: Paired-end status was given as %s. Ignoring all inputs for --input-fastq-r2!\n' % args.pair_status)
        ### get read group names
        rgs = []
        for fq in fqr1:
            rg = fq
            if rg.lower().endswith('.gz'):
                rg = rg[:-3]
            if rg.lower().endswith('.bz2'):
                rg = rg[:-4]
            if rg.lower().endswith('.fastq'):
                rg = rg[:-6]
            rgs.append(rg)
    ### this should not happen
    else:
        sys.exit('Error: need to specify either uBAM input via --input-ubam or fastq input via --input-fastq-r1/--input-fastq-r2!\n')

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
           '--outFileNamePrefix', args.sample + '_',
           '--readFilesIn', ','.join(fqr1), ','.join(fqr2),
           '--outSAMattrRGline', ' , '.join(['ID:%s' % _ for _ in rgs]),
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
           '--limitBAMsortRAM', '7000000000',
           '--outSAMattributes NH HI NM MD AS XS',
           '--outSAMunmapped Within',
           '--limitSjdbInsertNsj', '2000000',
           '--outSAMtype BAM Unsorted',
           '--outSAMheaderHD', '@HD VN:1.4',
           '--outSAMmultNmax', '1',
          ]
    subprocess.run(' '.join(cmd), shell=True, check=True)


if __name__ == "__main__":
    main()
