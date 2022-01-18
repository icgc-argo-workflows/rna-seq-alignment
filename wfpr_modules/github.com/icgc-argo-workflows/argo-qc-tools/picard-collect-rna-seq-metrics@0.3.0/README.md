

# Usage
```
Use command to invoke the tool directly

nextflow run main.nf --aligned_seq *** --ref_flat *** --ignore_seq *** --strand *** --ribosomal_interval_list ***

Arguments:
    --aligned_seq              Input SAM or BAM file. Required
    --ref_flat                 Gene annotations in refFlat form. Required 
    --ignore_seq               Seqence specified are ignored. Default: []
    --strand                   For strand-specific library prep. Default: "NONE". Allowed values are: ["NONE",
                               "FIRST_READ_TRANSCRIPTION_STRAND", "SECOND_READ_TRANSCRIPTION_STRAND"]
    --ribosomal_interval_list  Location of rRNA sequences in genome, in interval_list format. Default: null

```

# How to get the output files
You can specify the param `publish_dir` to a folder like `myout`, 
```
nextflow run main.nf --aligned_seq *** --ref_flat *** --ignore_seq *** --strand *** --ribosomal_interval_list *** --publish_dir myout
```
You will find the output files under the given folder. E.g.,
```
myout
└── picardCollectRnaSeqMetrics
    └── <aligned_seq>.collectrnaseqmetrics.tgz
```

# Output explanation
```
tar -tzf <aligned_seq>.collectrnaseqmetrics.tgz
        --rna_metrics.txt         Original tool output
        --qc_metrics.json         Information retrieved from the original output
```

# How to get the params.ref_flat.
Param `ref_flat` provide the gene annotations in refFlat form. You can download the file from [URL](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz)


# How to generate the params.ribosomal_interval_list
Param `ribosomal_interval_list` provide the locations of rRNA sequences in genome in interval_list format. If not specified no bases will be identified as being ribosomal. You can find more details about the format [here](https://gatk.broadinstitute.org/hc/en-us/articles/360037057492-CollectRnaSeqMetrics-Picard-#--RIBOSOMAL_INTERVALS)


You can use the provided script [`make_rRNA.sh`](./make_rRNA.sh) to generate the `ribosomal_interval_list` file which is suitable for `Picard:CollectRnaSeqMetrics`

# How to run the testing jobs
In order to run the testing jobs successfully, you will need to 
- Have `root` or `sudo` privileges of the machine
- Have `docker` installed and been up running

You can following the steps to do the testing
```
cd tests
nextflow run checker.nf -params-file test-job-1.json
nextflow run checker.nf -params-file test-job-2.json
```

