#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// processes resources
params.cpus = 8
params.mem = 20

params.publish_dir = ""

params.max_retries = 5  // set to 0 will disable retry
params.first_retry_wait_time = 1  // in seconds

// required params w/ default
params.container_version = "5.0.0"
params.transport_mem = 2 // Transport memory is in number of GBs

// optional if secret mounted from pod else required
params.api_token = "" // song/score API token for download process

// required params, no default
// --song_url         song url for download process
// --score_url        score url for download process

// TODO: Replace with score container once it can download files via analysis_id
process scoreDownload {
    maxRetries params.max_retries
    errorStrategy {
        sleep(Math.pow(2, task.attempt) * params.first_retry_wait_time * 1000 as long);  // backoff time increases exponentially before each retry
        return params.max_retries ? 'retry' : 'finish'
    }

    pod = [secret: workflow.runName + "-secret", mountPath: "/tmp/rdpc_secret"]
    
    cpus params.cpus
    memory "${params.mem} GB"
 
    container "overture/score:${params.container_version}"
    publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: params.publish_dir

    label "scoreDownload"
    tag "${analysis_id}"

    input:
        path analysis
        val study_id
        val analysis_id

    output:
        path analysis, emit: analysis_json
        path 'out/*', emit: files


    script:
        accessToken = params.api_token ? params.api_token : "`cat /tmp/rdpc_secret/secret`"
        """
        export METADATA_URL=${params.song_url}
        export STORAGE_URL=${params.score_url}
        export TRANSPORT_PARALLEL=${params.cpus}
        export TRANSPORT_MEMORY=${params.transport_mem}
        export ACCESSTOKEN=${accessToken}
        
        score-client download --analysis-id ${analysis_id} --study-id ${study_id} --output-dir ./out 
        """
}
