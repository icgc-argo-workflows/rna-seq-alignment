FROM ubuntu:20.04

LABEL org.opencontainers.image.source https://github.com/icgc-argo-workflows/data-processing-utility-tools

RUN groupadd -g 1000 ubuntu && \
    useradd -l -u 1000 -g ubuntu ubuntu && \
    install -d -m 0755 -o ubuntu -g ubuntu /home/ubuntu

CMD ["/bin/bash"]
