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
version = '2.7.0'

// universal params go here, change default value as needed
params.publish_dir = ""  // set to empty string will disable publishDir

params.max_retries = 5  // set to 0 will disable retry
params.first_retry_wait_time = 1  // in seconds

// tool specific parmas go here, add / change as needed
params.study_id = "TEST-PR"
params.payload = "NO_FILE"
params.upload = []
params.analysis_id = ""  // optional, analysis must already exist and in UNPUBLISHED state if analysis_id provided

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

include { songSubmit as songSub } from './local_modules/song-submit' params(song_params)
include { songManifest as songMan } from './local_modules/song-manifest' params(song_params)
include { scoreUpload as scoreUp } from './local_modules/score-upload' params(score_params)
include { songPublish as songPub } from './local_modules/song-publish' params(song_params)


workflow SongScoreUpload {
    take:
        study_id
        payload
        upload
        analysis_id

    main:
        if (!analysis_id) {
          // Create new analysis
          songSub(study_id, payload)
          analysis_id = songSub.out
        }

        // Generate file manifest for upload
        songMan(study_id, analysis_id, upload.collect())

        // Upload to SCORE
        scoreUp(analysis_id, songMan.out, upload.collect())

        // Publish the analysis
        songPub(study_id, scoreUp.out.ready_to_publish)

    emit:
        analysis_id = songPub.out.analysis_id
}


// this provides an entry point for this main script, so it can be run directly without clone the repo
// using this command: nextflow run <git_acc>/<repo>/<pkg_name>/<main_script>.nf -r <pkg_name>.v<pkg_version> --params-file xxx

workflow {
  SongScoreUpload(
    params.study_id,
    file(params.payload),
    Channel.fromPath(params.upload),
    params.analysis_id
  )
}
