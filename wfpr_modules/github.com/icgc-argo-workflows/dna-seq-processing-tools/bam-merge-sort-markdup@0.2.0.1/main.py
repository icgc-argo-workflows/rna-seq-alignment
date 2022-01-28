#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
  Copyright (C) 2021,  icgc-argo

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
    Junjun Zhang
    Linda Xiang
"""

import sys
import subprocess
import argparse
from multiprocessing import cpu_count
import json
import os

def run_cmd(cmd):
  stdout, stderr, p, success = '', '', None, True
  try:
    p = subprocess.Popen([cmd],
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          shell=True)
    stdout, stderr = p.communicate()
  except Exception as e:
    print('Execution failed: %s' % e, file=sys.stderr)
    success = False

  if p and p.returncode != 0:
    print('Execution failed, none zero code returned.', file=sys.stderr)
    success = False

  print(stdout.decode("utf-8"))
  print(stderr.decode("utf-8"), file=sys.stderr)

  if not success:
    sys.exit(p.returncode if p.returncode else 1)

  return stdout, stderr

def main():
    """ Main program """
    parser = argparse.ArgumentParser(description='Merge and markdup')
    parser.add_argument('-i','--input-bams', dest='input_bams',
                        type=str, help='Input bam file', nargs='+', required=True)
    parser.add_argument('-b','--output-base', dest='output_base',
                        type=str, help='Output merged file basename', required=True)
    parser.add_argument('-r', '--reference', dest='reference',
                        type=str, help='reference fasta', required=True)
    parser.add_argument('-t', '--tempdir', dest='tempdir', type=str, default=".",
                        help='Specify directory for temporary files')
    parser.add_argument("-n", "--cpus", dest='cpus', type=int, default=cpu_count())
    parser.add_argument("-d", "--mdup", dest='mdup', action='store_true')
    parser.add_argument("-l", "--lossy", dest='lossy', action='store_true')
    parser.add_argument("-o", "--output-format", dest='output_format', default='cram', choices=['bam', 'cram'])

    args = parser.parse_args()

    cmd = []

    if not os.path.isdir(args.tempdir):
        sys.exit('Error: specified tempdir %s does not exist!' % args.tempdir)

    if args.mdup:
        merge = 'bammarkduplicates2 markthreads=%s tmpfile=%s/tmp level=0 O=/dev/stdout M=%s I=%s ' % \
                (str(args.cpus), args.tempdir, args.output_base + ".duplicates_metrics.txt", ' I='.join(args.input_bams))
    else:
        merge = 'samtools merge --no-PG -uf -@ %s /dev/stdout %s ' % (str(args.cpus), ' '.join(args.input_bams))

    if args.lossy:
        cram = 'java -jar /tools/cramtools.jar cram -R %s --capture-all-tags --lossy-quality-score-spec \*8 --preserve-read-names -O %s' % (args.reference, args.output_base + ".cram")
    else:
        cram = 'samtools view -C -T %s -@ %s --write-index /dev/stdin -o %s ' % (args.reference, args.cpus, args.output_base + ".cram")

    bam = 'samtools view -b -h -@ %s --write-index /dev/stdin -o %s##idx##%s ' % (args.cpus, args.output_base + ".bam", args.output_base + ".bam.bai")
    crai1 = 'samtools index -@ %s %s %s ' % (args.cpus, args.output_base + ".cram", args.output_base + ".cram.crai")

    # build command
    if args.output_format == 'bam':
        cmd.append('|'.join([merge, bam]))

    elif args.output_format == 'cram':
        cmd.append('|'.join([merge, cram]))
        if args.lossy: cmd.append(crai1)
    else:
        sys.exit("Unsupported sequence format!")

    for c in cmd:
        run_cmd(c)

    if os.path.isfile(os.path.join(os.getcwd(), args.output_base + ".duplicates_metrics.txt")):
        stdout, _ = run_cmd('bammarkduplicates2  -v 2>&1  | grep "biobambam2 version"')
        version = stdout.decode("utf-8").split(' ')[-1].strip().rstrip('.')
        with open("%s.duplicates_metrics.extra_info.json" % args.output_base, "w") as j:
          j.write(json.dumps({  "tool": "biobambam2:bammarkduplicates2@%s" % version }, indent=2))

        tgz = 'tar czf %s.duplicates_metrics.tgz %s.duplicates_metrics.*' % (args.output_base, args.output_base)
        run_cmd(tgz)


if __name__ == "__main__":
    main()

