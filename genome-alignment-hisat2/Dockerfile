FROM ratschlab/wf-base-miniconda4.8.3:v1.0.1

LABEL org.opencontainers.image.source https://github.com/icgc-argo-rna-wg/rna-seq-alignment

ENV PATH="/tools:${PATH}"

COPY *.py /tools/

ENTRYPOINT ["/usr/bin/env"]

CMD ["/bin/bash"]
