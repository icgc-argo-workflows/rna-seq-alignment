FROM quay.io/icgc-argo/dna-seq-processing-tools:base-docker.0.2.1

LABEL org.opencontainers.image.source https://github.com/icgc-argo-workflows/dna-seq-processing-tools

ENV PATH="/tools:${PATH}"

COPY *.py /tools/

ENTRYPOINT ["/usr/bin/env"]

CMD ["/bin/bash"]
