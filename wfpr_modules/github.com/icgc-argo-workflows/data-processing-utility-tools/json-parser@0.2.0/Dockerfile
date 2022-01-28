FROM python:3.7.5-slim-buster

LABEL org.opencontainers.image.source https://github.com/icgc-argo-workflows/data-processing-utility-tools

RUN apt-get update && apt-get install -y procps jq

ENV PATH="/tools:${PATH}"

WORKDIR /tools

COPY . .
