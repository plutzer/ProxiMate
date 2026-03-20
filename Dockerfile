# syntax=docker/dockerfile:1

FROM ubuntu:20.04

LABEL maintainer="plutzer@wustl.edu"

# Set environment variable to prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install necessary build tools and dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    perl \
    wget \
    unzip \
    libboost-all-dev \
    libnlopt-dev

# RUN apt-get update && apt-get install -y perl

ADD ./SAINTexpress-custom/bin/SAINTexpress-int /bin/SAINTexpress-int_oldimp

COPY SAINTexpress-custom /SAINTexpress-custom

COPY SAINTexpress_v3.6.3__2018-03-09 /SAINTexpress_v3.6.3__2018-03-09
RUN find /SAINTexpress_v3.6.3__2018-03-09 \( -name "*.sh" -o -name "configure" -o -name "bootstrap" -o -name "b2" -o -name "bjam" \) -exec chmod +x {} +
RUN mkdir -p /SAINTexpress_v3.6.3__2018-03-09/bin

RUN make -C /SAINTexpress_v3.6.3__2018-03-09
RUN mv /SAINTexpress_v3.6.3__2018-03-09/bin/SAINTexpress-int /bin/SAINTexpress-int_default
RUN mv /SAINTexpress_v3.6.3__2018-03-09/bin/SAINTexpress-spc /bin/SAINTexpress-spc

# Overwrite the default SAINTexpress with the custom version
COPY /SAINTexpress-custom/SAINT-MRF-int/*.cpp /SAINTexpress_v3.6.3__2018-03-09/SAINT-MRF-int/
COPY /SAINTexpress-custom/SAINT-MRF-int/*.hpp /SAINTexpress_v3.6.3__2018-03-09/SAINT-MRF-int/
# COPY /SAINTexpress-custom/SAINT-MRF-int/Makefile /SAINTexpress_v3.6.3__2018-03-09/SAINT-MRF-int/
RUN make -C /SAINTexpress_v3.6.3__2018-03-09/SAINT-MRF-int clean
RUN make -C /SAINTexpress_v3.6.3__2018-03-09/SAINT-MRF-spc clean

# Build the project
RUN make -C /SAINTexpress_v3.6.3__2018-03-09/SAINT-MRF-int
RUN cp /SAINTexpress_v3.6.3__2018-03-09/bin/SAINTexpress-int /bin/SAINTexpress-int

RUN apt-get update && apt-get install -y python3-pip

COPY requirements.txt ../requirements.txt
RUN pip3 install -r ../requirements.txt

# Copy CORUM dataset (can't be downloaded programmatically) and download the rest
COPY Datasets /Datasets
COPY Scripts/preprocess_biogrid.py /Scripts/preprocess_biogrid.py
COPY Scripts/setup_datasets.py /Scripts/setup_datasets.py
RUN python3 /Scripts/setup_datasets.py --output-dir /Datasets --skip corum

# Copy over the folders with code and tests
COPY Scripts /Scripts
COPY GUI /GUI
COPY run_pipeline.sh /run_pipeline.sh

ENV PYTHONPATH=/Scripts:/Scripts/GOGO
ENV PATH="/Scripts/GOGO:${PATH}"

# Create a directory for output files
RUN mkdir -p /Outputs

# When the container starts, start the GUI
EXPOSE 3838
CMD ["python3", "-u", "/GUI/app.py"]

