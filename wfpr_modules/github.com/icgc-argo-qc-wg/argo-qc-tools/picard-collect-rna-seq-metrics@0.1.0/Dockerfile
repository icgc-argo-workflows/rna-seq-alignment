FROM ubuntu:20.04

LABEL org.opencontainers.image.source https://github.com/icgc-argo-qc-wg/argo-qc-tools

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y && \
    apt-get install -y software-properties-common python3-pip python3-dev curl r-base

RUN mkdir /tools

# Install OpenJDK-8
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk && \
    apt-get install -y ant && \
    apt-get clean

# Picard installation
ARG picard_version=2.25.6
RUN curl -sSL -o /tools/picard.jar --retry 10 https://github.com/broadinstitute/picard/releases/download/${picard_version}/picard.jar

ENV PATH="/tools:${PATH}"

COPY *.py /tools/

RUN groupadd -g 1000 ubuntu && \
    useradd -l -u 1000 -g ubuntu ubuntu && \
    install -d -m 0755 -o ubuntu -g ubuntu /home/ubuntu

USER ubuntu

ENTRYPOINT ["/usr/bin/env"]

CMD ["/bin/bash"]
