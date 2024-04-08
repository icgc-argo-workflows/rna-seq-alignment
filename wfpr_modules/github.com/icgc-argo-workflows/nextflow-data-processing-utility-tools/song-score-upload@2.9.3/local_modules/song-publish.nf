#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// processes resources
params.cpus = 1
params.mem = 1

params.max_retries = 5  // set to 0 will disable retry
params.first_retry_wait_time = 1  // in seconds

// required params w/ default
params.container = "ghcr.io/overture-stack/song-client"
params.container_version = "latest"

// optional if secret mounted from pod else required
params.api_token = "" // song/score API token for download process

// required params, no default
// --song_url         song url for download process
// --score_url        score url for download process

process songPublish {
    maxRetries params.max_retries
    errorStrategy {
        sleep(Math.pow(2, task.attempt) * params.first_retry_wait_time * 1000 as long);  // backoff time doubles before each retry
        return params.max_retries ? 'retry' : 'finish'
    }

    pod = [secret: workflow.runName + "-secret", mountPath: "/tmp/rdpc_secret"]
    
    cpus params.cpus
    memory "${params.mem} GB"
 
    container "${ params.song_container ?: params.container}:${params.song_container_version ?: params.container_version}"

    if (workflow.containerEngine == "singularity") {
        containerOptions "--bind \$(pwd):/song-client/logs"
    } else if (workflow.containerEngine == "docker") {
        containerOptions "-v \$(pwd):/song-client/logs"
    }

    tag "${analysis_id}"
    
    input:
        val study_id
        val analysis_id

    output:
        val analysis_id, emit: analysis_id

    script:
        accessToken = params.api_token ? params.api_token : "`cat /tmp/rdpc_secret/secret`"
        """
        export CLIENT_SERVER_URL=${params.song_url}
        export CLIENT_STUDY_ID=${study_id}
        export CLIENT_ACCESS_TOKEN=${accessToken}

        sing publish -a  ${analysis_id}
        """
}
