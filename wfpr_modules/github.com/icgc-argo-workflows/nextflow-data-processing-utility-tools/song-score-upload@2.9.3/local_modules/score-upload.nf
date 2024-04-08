#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// processes resources
params.cpus = 8
params.mem = 20

params.max_retries = 5  // set to 0 will disable retry
params.first_retry_wait_time = 1  // in seconds

// required params w/ default
params.container = "ghcr.io/overture-stack/score"
params.container_version = "latest"
params.transport_mem = 2 // Transport memory is in number of GBs

// optional if secret mounted from pod else required
params.api_token = "" // song/score API token for download process
params.score_force = false
// required params, no default
// --song_url         song url for download process
// --score_url        score url for download process

process scoreUpload {
    maxRetries params.max_retries
    errorStrategy {
        sleep(Math.pow(2, task.attempt) * params.first_retry_wait_time * 1000 as long);  // backoff time doubles before each retry
        return params.max_retries ? 'retry' : 'finish'
    }

    pod = [secret: workflow.runName + "-secret", mountPath: "/tmp/rdpc_secret"]
    
    cpus params.cpus
    memory "${params.mem} GB"
  
    container "${ params.score_container ?: params.container}:${params.score_container_version ?: params.container_version}"

    if (workflow.containerEngine == "singularity") {
        containerOptions "--bind \$(pwd):/score-client/logs"
    } else if (workflow.containerEngine == "docker") {
        containerOptions "-v \$(pwd):/score-client/logs"
    }

    tag "${analysis_id}"

    input:
        val analysis_id
        path manifest
        path upload

    output:
        val analysis_id, emit: ready_to_publish

    script:
        accessToken = params.api_token ? params.api_token : "`cat /tmp/rdpc_secret/secret`"
        forceFlag = params.score_force ? "--force" : ""
        """
        export METADATA_URL=${params.song_url}
        export STORAGE_URL=${params.score_url}
        export TRANSPORT_PARALLEL=${params.cpus}
        export TRANSPORT_MEMORY=${params.transport_mem}
        export ACCESSTOKEN=${accessToken}
        
        score-client upload --manifest ${manifest} ${forceFlag}
        """
}
