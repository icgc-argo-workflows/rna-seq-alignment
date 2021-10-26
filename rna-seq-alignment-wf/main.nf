#!/usr/bin/env nextflow

/*
  Copyright (c) 2021, icgc-argo-rna-wg

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

  Authors:
    Linda Xiang
*/

nextflow.enable.dsl = 2
name = 'rna-seq-alignment'
version = '0.1.0'  // package version

// universal params go here, change default value as needed
params.container = ""
params.container_registry = ""
params.container_version = ""
params.cpus = 1
params.mem = 1  // GB
params.publish_dir = ""  // set to empty string will disable publishDir

// tool specific parmas go here
params.study_id = ""
params.analysis_id = ""
params.ref_genome_fa = ""
params.ref_genome_gtf = ""
params.genome_annotation = "GENCODE v38"
params.genome_build = "GRCh38_hla_decoy_ebv"
params.ref_flat = ""
params.ignore_seq = "NO_FILE_ignore"
params.ribosomal_interval_list = "NO_FILE_interval"
params.strand = ""


// if provided local mode will be used
params.analysis_metadata = "NO_FILE"
params.experiment_info_tsv = "NO_FILE1"
params.read_group_info_tsv = "NO_FILE2"
params.file_info_tsv = "NO_FILE3"
params.extra_info_tsv = "NO_FILE4"
params.sequencing_files = []

params.cleanup = true
params.max_retries = 5  // set to 0 will disable retry
params.first_retry_wait_time = 1  // in seconds

params.tempdir = "NO_DIR"

params.song_url = ""
params.score_url = ""
params.api_token = ""
params.download = [:]
params.seqDataToLaneBam = [:]
params.starAligner = [:]
params.hisat2Aligner = [:]
params.bamMergeSortMarkdup = [:]
params.uploadAlignment = [:]
params.readGroupUBamQC = [:]
params.payloadGenRnaAlignment = [:]
params.alignedSeqQC = [:]
params.payloadGenRnaSeqQc = [:]
params.uploadQc = [:]
//params.gatkCollectOxogMetrics = [:]

download_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'max_retries': params.max_retries,
    'first_retry_wait_time': params.first_retry_wait_time,
    'song_url': params.song_url,
    'score_url': params.score_url,
    'api_token': params.api_token,
    'publish_dir': params.publish_dir,
    *:(params.download ?: [:])
]

seqDataToLaneBam_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'reads_max_discard_fraction': -1,
    'publish_dir': params.publish_dir,
    *:(params.seqDataToLaneBam ?: [:])
]

starAligner_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'tempdir': params.tempdir ?: 'NO_DIR',
    'publish_dir': params.publish_dir,
    *:(params.starAligner ?: [:])
]

hisat2Aligner_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'tempdir': params.tempdir ?: 'NO_DIR',
    'publish_dir': params.publish_dir,
    *:(params.hisat2Aligner ?: [:])
]

bamMergeSortMarkdup_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'output_format': 'cram',
    'markdup': true,
    'lossy': false,
    'tempdir': params.tempdir ?: 'NO_DIR',
    'aligned_basename': 'genome.merged',
    'publish_dir': params.publish_dir,
    *:(params.bamMergeSortMarkdup ?: [:])
]

readGroupUBamQC_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'publish_dir': params.publish_dir,
    *:(params.readGroupUBamQC ?: [:])
]

payloadGenRnaAlignment_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'publish_dir': params.publish_dir,
    *:(params.payloadGenRnaAlignment ?: [:])
]

alignedSeqQC_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'publish_dir': params.publish_dir,
    *:(params.alignedSeqQC ?: [:])
]

payloadGenRnaSeqQc_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'publish_dir': params.publish_dir,
    *:(params.payloadGenRnaSeqQc ?: [:])
]

uploadAlignment_params = [
    'max_retries': params.max_retries,
    'first_retry_wait_time': params.first_retry_wait_time,
    'cpus': params.cpus,
    'mem': params.mem,
    'song_url': params.song_url,
    'score_url': params.score_url,
    'api_token': params.api_token,
    *:(params.uploadAlignment ?: [:])
]

uploadQc_params = [
    'max_retries': params.max_retries,
    'first_retry_wait_time': params.first_retry_wait_time,
    'cpus': params.cpus,
    'mem': params.mem,
    'song_url': params.song_url,
    'score_url': params.score_url,
    'api_token': params.api_token,
    *:(params.uploadQc ?: [:])
]

include { seqDataToLaneBam as toLaneBam } from "./modules/raw.githubusercontent.com/icgc-argo-workflows/dna-seq-processing-tools/seq-data-to-lane-bam.0.3.3.0/tools/seq-data-to-lane-bam/seq-data-to-lane-bam.nf" params(seqDataToLaneBam_params)
include { icgcArgoRnaSeqAlignmentSTAR as star } from "./wfpr_modules/github.com/icgc-argo-rna-wg/rna-seq-alignment/genome-alignment-star@0.2.0/alignSTAR.nf" params(starAligner_params)
include { icgcArgoRnaSeqAlignmentHISAT2 as hisat2 } from "./wfpr_modules/github.com/icgc-argo-rna-wg/rna-seq-alignment/genome-alignment-hisat2@0.2.0/alignHISAT2.nf" params(hisat2Aligner_params)
include { bamMergeSortMarkdup as merMkdupStar; bamMergeSortMarkdup as merMkdupHisat2} from "./wfpr_modules/github.com/icgc-argo-workflows/dna-seq-processing-tools/bam-merge-sort-markdup@0.2.0.1/main.nf" params(bamMergeSortMarkdup_params)
include { fastqc } from "./wfpr_modules/github.com/icgc-argo-qc-wg/argo-qc-tools/fastqc@0.1.0/main.nf" params(readGroupUBamQC_params)
include { picardCollectRnaSeqMetrics as alignedSeqQcStar;  picardCollectRnaSeqMetrics as alignedSeqQcHisat2} from "./wfpr_modules/github.com/icgc-argo-qc-wg/argo-qc-tools/picard-collect-rna-seq-metrics@0.1.0/main.nf" params(alignedSeqQC_params)
include { getSecondaryFiles as getSec } from './wfpr_modules/github.com/icgc-argo-workflows/data-processing-utility-tools/helper-functions@1.0.1.1/main.nf'
include { cleanupWorkdir as cleanup } from './wfpr_modules/github.com/icgc-argo-workflows/data-processing-utility-tools/cleanup-workdir@1.0.0.1/main.nf'
include { payloadGenSeqExperiment as pGenExp } from './wfpr_modules/github.com/icgc-argo-workflows/data-processing-utility-tools/payload-gen-seq-experiment@0.5.0.1/main.nf'
include { cram2bam as cram2bamStar; cram2bam as cram2bamHisat2 } from './wfpr_modules/github.com/icgc-argo-workflows/dna-seq-processing-tools/cram2bam@0.1.0/main.nf'
include { payloadGenRnaAlignment as pGenAlnStar;  payloadGenRnaAlignment as pGenAlnHisat2; payloadGenRnaAlignment as pGenAlnStarSj;  payloadGenRnaAlignment as pGenAlnHisat2Sj; payloadGenRnaAlignment as pGenQcStar; payloadGenRnaAlignment as pGenQcHisat2 } from './wfpr_modules/github.com/icgc-argo-workflows/data-processing-utility-tools/payload-gen-rna-alignment@0.1.0/main.nf'


// please update workflow code as needed
workflow RnaSeqAlignmentWf {
  take:  
    study_id
    analysis_id
    ref_genome_fa
    ref_genome_gtf
    analysis_metadata
    experiment_info_tsv
    read_group_info_tsv
    file_info_tsv
    extra_info_tsv
    sequencing_files


  main:
    // detect local mode or not
    local_mode = false
    if ((!analysis_metadata.startsWith("NO_FILE") || !experiment_info_tsv.startsWith("NO_FILE")) && sequencing_files.size() > 0){
        local_mode = true
        if (!params.publish_dir) {
            exit 1, "You specified local sequencing data as input, please also set `params.publish_dir` to keep the output."
        }
        log.info "Run the workflow using local input sequencing data, alignment results will be in: ${params.publish_dir}"

        if (!analysis_metadata.startsWith("NO_FILE")) {
            if (!experiment_info_tsv.startsWith("NO_FILE") ||
                !read_group_info_tsv.startsWith("NO_FILE") ||
                !file_info_tsv.startsWith("NO_FILE") ||
                !extra_info_tsv.startsWith("NO_FILE")
            )  {
                log.info "Use analysis metadata JSON as input, will ignore input: 'experiment_info_tsv', 'read_group_info_tsv', 'file_info_tsv', 'extra_info_tsv'"
            }
            analysis_metadata = file(analysis_metadata)
        } else if (!experiment_info_tsv.startsWith("NO_FILE") &&
                  !read_group_info_tsv.startsWith("NO_FILE") &&
                  !file_info_tsv.startsWith("NO_FILE") &&
                  !extra_info_tsv.startsWith("NO_FILE")
            ) {
            pGenExp(
                file(experiment_info_tsv),
                file(read_group_info_tsv),
                file(file_info_tsv),
                file(extra_info_tsv)
            )
            analysis_metadata = pGenExp.out.payload
        } else {
            exit 1, "To run the workflow using local inputs, please specify metadata in JSON using params.analysis_metadata or metadata in TSVs using params.experiment_info_tsv, params.read_group_info_tsv, params.file_info_tsv and params.extra_info_tsv"
        }

        sequencing_files = Channel.fromPath(sequencing_files)
    } else if (study_id && analysis_id) {
        // download files and metadata from song/score (analysis type: sequencing_experiment)
        log.info "Run the workflow using input sequencing data from SONG/SCORE, alignment results will be uploaded to SONG/SCORE as well"
        dnld(study_id, analysis_id)
        analysis_metadata = dnld.out.analysis_json
        sequencing_files = dnld.out.files
    } else {
        exit 1, "To use sequencing data from SONG/SCORE as input, please provide `params.study_id`, `params.analysis_id` and other SONG/SCORE params.\n" +
            "Or please provide `params.analysis_metadata` (or `params.experiment_info_tsv`, `params.read_group_info_tsv`, `params.file_info_tsv` and `params.extra_info_tsv`) and `params.sequencing_files` from local files as input."
    }

    // preprocessing input data (BAM or FASTQ) into read group level unmapped BAM (uBAM)
    toLaneBam(analysis_metadata, sequencing_files.collect())

    // perform ubam QC
    fastqc(toLaneBam.out.lane_bams.flatten())

    // use scatter to run STAR alignment for each ubam in parallel
    star(file(params.ref_genome_index_star), file(ref_genome_gtf), toLaneBam.out.lane_bams.flatten(), params.input_format, 
        params.pair_status, params.sample, params.sjdboverhang)  // just to run after rgQC

    // use scatter to run HISAT2 alignment for each ubam in parallel
    hisat2(params.ref_genome_index_hisat2, file(params.ref_genome_index_hisat2).getParent(), file(params.ref_genome_gtf), 
           toLaneBam.out.lane_bams.flatten(), params.input_format, params.pair_status, params.sample)  // just to run after rgQC

    // collect aligned lane bams for merge and markdup
    merMkdupStar(star.out.bam.collect(), file(ref_genome_fa + '.gz'),
        Channel.fromPath(getSec(ref_genome_fa + '.gz', ['fai', 'gzi']), checkIfExists: true).collect(),
        bamMergeSortMarkdup_params.tempdir)
    
    // collect aligned lane bams for merge and markdup
    merMkdupHisat2(hisat2.out.bam.collect(), file(ref_genome_fa + '.gz'),
        Channel.fromPath(getSec(ref_genome_fa + '.gz', ['fai', 'gzi']), checkIfExists: true).collect(),
        bamMergeSortMarkdup_params.tempdir)

    // generate payload for aligned seq (analysis type: sequencing_alignment)
    pGenAlnStar(merMkdupStar.out.merged_seq.concat(merMkdupStar.out.merged_seq_idx).collect(),
        analysis_metadata, 'STAR', 'sequencing_alignment', params.genome_annotation, params.genome_build, name, version)

    // generate payload for splice junctions (analysis type: splice_junctions)
    pGenAlnStarSj(star.out.junctions.collect(),
        analysis_metadata, 'STAR', 'splice_junctions', params.genome_annotation, params.genome_build, name, version)

    // generate payload for aligned seq (analysis type: sequencing_alignment)
    pGenAlnHisat2(merMkdupHisat2.out.merged_seq.concat(merMkdupHisat2.out.merged_seq_idx).collect(),
        analysis_metadata, 'HiSAT2', 'sequencing_alignment', params.genome_annotation, params.genome_build, name, version)

    // generate payload for splice junctions (analysis type: splice_junctions)
    pGenAlnHisat2Sj(hisat2.out.junctions.collect(),
        analysis_metadata, 'HiSAT2', 'splice_junctions', params.genome_annotation, params.genome_build, name, version)

    // // upload aligned file and metadata to song/score
    // def alnAnalysisId
    // if (!local_mode) {
    //     upAln(study_id, pGenRnaAln.out.payload, pGenRnaAln.out.alignment_files.collect())
    //     alnAnalysisId = upAln.out.analysis_id
    // } else {
    //     alnAnalysisId = 'Unknown'
    // }

    // convert cram to bam
    cram2bamStar(pGenAlnStar.out.cram.flatten().first(), file(ref_genome_fa + '.gz'), 
             Channel.fromPath(getSec(ref_genome_fa + '.gz', ['fai', 'gzi']), checkIfExists: true).collect())

    // convert cram to bam
    cram2bamHisat2(pGenAlnHisat2.out.cram.flatten().first(), file(ref_genome_fa + '.gz'), 
             Channel.fromPath(getSec(ref_genome_fa + '.gz', ['fai', 'gzi']), checkIfExists: true).collect())

    // perform aligned seq QC
    alignedSeqQcStar(cram2bamStar.out.output_bam, file(params.ref_flat),
        file(params.ignore_seq), file(params.ribosomal_interval_list), params.strand)  // run after upAln

    // perform aligned seq QC
    alignedSeqQcHisat2(cram2bamHisat2.out.output_bam, file(params.ref_flat),
        file(params.ignore_seq), file(params.ribosomal_interval_list), params.strand)  // run after upAln

    // prepare song payload for qc metrics
    pGenQcStar(alignedSeqQcStar.out.qc_tar.concat(
            merMkdupStar.out.duplicates_metrics,
            fastqc.out.qc_tar).collect(),
            analysis_metadata, 'STAR', 'qc_metrics', 
            params.genome_annotation, params.genome_build, name, version)

    // prepare song payload for qc metrics
    pGenQcHisat2(alignedSeqQcHisat2.out.qc_tar.concat(
            merMkdupHisat2.out.duplicates_metrics).collect(),
            analysis_metadata, 'HiSAT2', 'qc_metrics', 
            params.genome_annotation, params.genome_build, name, version)

    // // upload aligned file and metadata to song/score
    // if (!local_mode) {
    //     upQc(study_id, pGenRnaSeqQc.out.payload, pGenRnaSeqQc.out.qc_files.collect())
    // }

    // if (params.cleanup && !local_mode) {
    //     cleanup(
    //         sequencing_files.concat(toLaneBam.out, star.out,  merMkdupStar.out,
    //             alignedSeqQC.out, oxog.out, rgQC.out).collect(),
    //         upAln.out.analysis_id.concat(upQc.out.analysis_id).collect())  // wait until upAln and upQc is done
    // } else if (params.cleanup && local_mode) {
    //     cleanup(
    //         sequencing_files.concat(toLaneBam.out, star.out, merMkdupStar.out,
    //             alignedSeqQC.out, oxog.out, rgQC.out).collect(), pGenRnaSeqQc.out.payload)
    // }

  emit:  // update as needed
    star_aln_payload = pGenAlnStar.out.payload
    star_aln_files = pGenAlnStar.out.files_to_upload
    star_qc_metrics_payload = pGenQcStar.out.payload
    star_qc_metrics_files = pGenQcStar.out.files_to_upload
    hisat2_aln_payload = pGenAlnHisat2.out.payload
    hisat2_aln_files = pGenAlnHisat2.out.files_to_upload
    hisat2_qc_metrics_payload = pGenQcHisat2.out.payload
    hisat2_qc_metrics_files = pGenQcHisat2.out.files_to_upload
}


// this provides an entry point for this main script, so it can be run directly without clone the repo
// using this command: nextflow run <git_acc>/<repo>/<pkg_name>/<main_script>.nf -r <pkg_name>.v<pkg_version> --params-file xxx
workflow {
  RnaSeqAlignmentWf(
    params.study_id,
    params.analysis_id,
    params.ref_genome_fa,
    params.ref_genome_gtf,
    params.analysis_metadata,
    params.experiment_info_tsv,
    params.read_group_info_tsv,
    params.file_info_tsv,
    params.extra_info_tsv,
    params.sequencing_files
  )
}