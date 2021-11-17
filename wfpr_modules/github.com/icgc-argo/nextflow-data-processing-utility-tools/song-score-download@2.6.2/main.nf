#!/usr/bin/env nextflow

/*
  Copyright (c) 2020-2021, Ontario Institute for Cancer Research

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  Authors:
    Alex Lepsa
    Junjun Zhang
*/

nextflow.enable.dsl = 2
version = '2.6.2'

// universal params go here, change default value as needed
params.publish_dir = ""  // set to empty string will disable publishDir

params.max_retries = 5  // set to 0 will disable retry
params.first_retry_wait_time = 1  // in seconds

// tool specific parmas go here, add / change as needed
params.study_id = "TEST-PR"
params.analysis_id = "9940db0f-c100-496a-80db-0fc100d96ac1"

params.api_token = ""

params.song_cpus = 1
params.song_mem = 1  // GB
params.song_url = "https://song.rdpc-qa.cancercollaboratory.org"
params.song_api_token = ""
params.song_container_version = "4.2.1"

params.score_cpus = 1
params.score_mem = 1  // GB
params.score_transport_mem = 1  // GB
params.score_url = "https://score.rdpc-qa.cancercollaboratory.org"
params.score_api_token = ""
params.score_container_version = "5.0.0"


song_params = [
    *:params,
    'cpus': params.song_cpus,
    'mem': params.song_mem,
    'song_url': params.song_url,
    'song_container_version': params.song_container_version,
    'api_token': params.song_api_token ?: params.api_token
]

score_params = [
    *:params,
    'cpus': params.score_cpus,
    'mem': params.score_mem,
    'transport_mem': params.score_transport_mem,
    'song_url': params.song_url,
    'score_url': params.score_url,
    'score_container_version': params.score_container_version,
    'api_token': params.score_api_token ?: params.api_token
]


include { songGetAnalysis as songGet } from './local_modules/song-get-analysis' params(song_params)
include { scoreDownload as scoreDn } from './local_modules/score-download' params(score_params)


// please update workflow code as needed
workflow SongScoreDownload {
  take:  // update as needed
    study_id
    analysis_id

  main:
    songGet(study_id, analysis_id)
    scoreDn(songGet.out.json, study_id, analysis_id)

  emit:
    analysis_json = songGet.out.json
    files = scoreDn.out.files
}


// this provides an entry point for this main script, so it can be run directly without clone the repo
// using this command: nextflow run <git_acc>/<repo>/<pkg_name>/<main_script>.nf -r <pkg_name>.v<pkg_version> --params-file xxx
workflow {
  SongScoreDownload(
    params.study_id,
    params.analysis_id
  )
}
