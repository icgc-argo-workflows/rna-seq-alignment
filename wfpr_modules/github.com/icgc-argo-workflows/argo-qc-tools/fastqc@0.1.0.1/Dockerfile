FROM ubuntu:20.04

LABEL org.opencontainers.image.source https://github.com/icgc-argo-workflows/argo-qc-tools

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update \
  && apt-get install -y software-properties-common

RUN apt-get update \
	&& apt-get install -y openjdk-8-jre \
	&& rm -rf /var/lib/apt/lists/*

ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/

RUN apt-get -qq update \
  && apt-get -y upgrade \
  && apt install -y wget libfindbin-libs-perl software-properties-common unzip

# install fastqc v0.11.9
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip -O /tmp/fastqc.zip \
  && unzip /tmp/fastqc.zip -d /opt/ \
  && rm /tmp/fastqc.zip \
  && chmod 775 /opt/FastQC/fastqc

ENV PATH="/opt/FastQC/:${PATH}"

ENV PATH="/tools:${PATH}"

COPY *.py /tools/

RUN groupadd -g 1000 ubuntu && \
    useradd -l -u 1000 -g ubuntu ubuntu && \
    install -d -m 0755 -o ubuntu -g ubuntu /home/ubuntu

USER ubuntu

ENTRYPOINT ["/usr/bin/env"]

CMD ["/bin/bash"]
