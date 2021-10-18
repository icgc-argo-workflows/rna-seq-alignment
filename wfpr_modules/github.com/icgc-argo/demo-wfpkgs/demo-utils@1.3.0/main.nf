#!/usr/bin/env nextflow

/*
 * Copyright (c) 2019-2021, Ontario Institute for Cancer Research (OICR).
 *                                                                                                               
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/*
 * Contributors:
 *   Junjun Zhang <junjun.zhang@oicr.on.ca>
 */

/********************************************************************/
/* this block is auto-generated based on info from pkg.json where   */
/* changes can be made if needed, do NOT modify this block manually */
nextflow.enable.dsl = 2
name = 'demo-utils'
version = '1.3.0'
/********************************************************************/


// load local process (module)
include { cleanupWorkdir } from './local_modules/cleanup-workdir.nf'


// this is kind of like CWL's secondary files
def getSecondaryFiles(main_file, exts){
  if (!(main_file instanceof String)) {
    exit 1, "[getSecondaryFiles] param: main_file must be a string"
  }

  if (!(exts instanceof List)) {
    exit 1, "[getSecondaryFiles] param: exts must be a list of strings"
  }

  def secondaryFiles = []
  for (ext in exts) {
    if (ext.startsWith("^")) {
      ext = ext.replace("^", "")
      parts = main_file.split("\\.").toList()
      parts.removeLast()
      secondaryFiles.add((parts + [ext]).join("."))
    } else {
      secondaryFiles.add(main_file + '.' + ext)
    }
  }
  return secondaryFiles
}


// get specific secondary files for BWA alignment, ensure none is missing
def getBwaSecondaryFiles(main_file){
  return getSecondaryFiles(main_file, ['fai', 'sa', 'bwt', 'ann', 'amb', 'pac', 'alt'])
}
