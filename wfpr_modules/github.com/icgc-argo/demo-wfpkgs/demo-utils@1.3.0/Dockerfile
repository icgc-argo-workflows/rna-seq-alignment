FROM ubuntu:20.04

LABEL org.opencontainers.image.source https://github.com/icgc-argo/demo-wfpkgs

RUN groupadd -g 1000 ubuntu && \
    useradd -l -u 1000 -g ubuntu ubuntu && \
    install -d -m 0755 -o ubuntu -g ubuntu /home/ubuntu

USER ubuntu

CMD ["/bin/bash"]
