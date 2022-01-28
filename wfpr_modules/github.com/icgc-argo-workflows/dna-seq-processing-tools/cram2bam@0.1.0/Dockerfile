FROM ubuntu:20.04

LABEL org.opencontainers.image.source https://github.com/icgc-argo-workflows/dna-seq-processing-tools

RUN apt-get update -y && \
    apt-get install -y samtools

ENTRYPOINT ["/usr/bin/env"]

CMD ["/bin/bash"]
