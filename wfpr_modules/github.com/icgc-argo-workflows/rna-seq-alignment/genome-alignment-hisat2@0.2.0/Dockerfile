FROM ratschlab/wf-base-miniconda4.8.3:v1.0.1

LABEL org.opencontainers.image.source https://github.com/icgc-argo-workflows/icgc-argo-rna-seq-alignment
LABEL org.opencontainers.image.authors Andre Kahles (andre.kahles@inf.ethz.ch)
LABEL org.opencontainers.image.title ICGC ARGO RNA-Seq alignment image

ENV PATH="/opt/local:${PATH}"

RUN apt-get update
RUN apt-get -y install wget zip make gcc g++ libcurl4-openssl-dev libgomp1 bzip2 curl libncurses5-dev libbz2-dev liblzma-dev autotools-dev automake libz-dev pkg-config libtool m4 

# install HISAT2
RUN mkdir -p /opt/local/src/HISAT2 && \
    cd /opt/local/src/HISAT2 && \
    wget -O hisat2-2.2.1-source.zip https://cloud.biohpc.swmed.edu/index.php/s/fE9QCsX3NH4QwBi/download && \
    unzip hisat2-2.2.1-source.zip && \
    cd hisat2-2.2.1 && \
    make && \
    cp hisat2 hisat2-align-s hisat2-align-l hisat2-build hisat2-build-s hisat2-build-l hisat2-inspect hisat2-inspect-s hisat2-inspect-l /usr/local/bin/ && \
    cp *.py /opt/local/

# install samtools 1.10
RUN cd /tmp \
    && curl -sSL -o samtools-1.10.tar.bz2 --retry 10 https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 \
    && bunzip2 -c samtools-1.10.tar.bz2 |tar xf - \
    && cd samtools-1.10 \
    && ./configure --prefix=/usr/local \
    && make \
    && make install

# Install OpenJDK-8
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk && \
    apt-get install -y ant && \
    apt-get clean;

RUN mkdir -p /tools

# Picard installation
RUN curl -sSL -o /tools/picard.jar --retry 10 https://github.com/broadinstitute/picard/releases/download/2.21.8/picard.jar

### get run script
COPY alignHISAT2.py /tools

### set up environment
ENV PATH="$PATH:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:"
ENV LC_ALL=C

RUN groupadd -g 1000 ubuntu &&\
    useradd -l -u 1000 -g ubuntu ubuntu &&\
    install -d -m 0755 -o ubuntu -g ubuntu /home/ubuntu

USER ubuntu

ENTRYPOINT ["/usr/bin/env"]

CMD ["/bin/bash"]
