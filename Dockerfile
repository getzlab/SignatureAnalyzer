FROM gcr.io/broad-getzlab-workflows/base_image:v0.0.6
MAINTAINER Zachary Everton

WORKDIR build
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install --upgrade setuptools

WORKDIR /app
COPY . ./signatureanalyzer
#RUN pip install "signatureanalyzer==0.0.9"
