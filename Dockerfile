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
    libnlopt-dev \
    dos2unix

# RUN apt-get update && apt-get install -y perl

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

# Python 3.12 is built from source: focal's apt tops out at 3.9 and the deadsnakes
# PPA publishes nothing for this release, while the pinned requirements need >= 3.11.
# Building here rather than on a newer base image keeps gcc 9, which is what the
# vendored Boost 1.57 and nlopt 2.3 compiled above still build under.
# altinstall installs the interpreter as python3.12; the /usr/local/bin symlinks are
# what make `python3` and `pip3` resolve to it for every later stage and for the
# subprocess calls in GUI/app.py and run_pipeline.sh.
ARG PYTHON_VERSION=3.12.14
RUN apt-get update && apt-get install -y \
        zlib1g-dev libssl-dev libffi-dev libbz2-dev \
        libreadline-dev libsqlite3-dev liblzma-dev ca-certificates \
 && wget -q https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz \
 && tar -xf Python-${PYTHON_VERSION}.tgz \
 && cd Python-${PYTHON_VERSION} \
 && ./configure --quiet \
 && make -j"$(nproc)" \
 && make altinstall \
 && cd / && rm -rf /Python-${PYTHON_VERSION} /Python-${PYTHON_VERSION}.tgz \
 && python3.12 -m ensurepip --upgrade \
 && python3.12 -m pip install --upgrade pip \
 && ln -sf /usr/local/bin/python3.12 /usr/local/bin/python3 \
 && ln -sf /usr/local/bin/pip3.12 /usr/local/bin/pip3

COPY requirements.txt ../requirements.txt
RUN pip3 install -r ../requirements.txt

# CORUM is tracked in git; the other databases are downloaded here unless the
# build context already holds them.  setup_datasets.py keeps existing files, which
# is how the Release workflow gives the amd64 and arm64 builds one snapshot.
COPY Datasets /Datasets
COPY Scripts/preprocess_biogrid.py /Scripts/preprocess_biogrid.py
COPY Scripts/setup_datasets.py /Scripts/setup_datasets.py
RUN python3 /Scripts/setup_datasets.py --output-dir /Datasets --skip corum

# Copy over the folders with code and tests
COPY Scripts /Scripts
COPY GUI /GUI
COPY run_pipeline.sh /run_pipeline.sh
RUN dos2unix /run_pipeline.sh

ENV PYTHONPATH=/Scripts:/Scripts/GOGO
ENV PATH="/Scripts/GOGO:${PATH}"

# Identifies the code in every run.json.  The image carries no .git, so this is
# the only place the version can come from:
#   docker build --build-arg PROXIMATE_VERSION=$(git rev-parse --short HEAD) ...
ARG PROXIMATE_VERSION=unknown
ENV PROXIMATE_VERSION=${PROXIMATE_VERSION}
LABEL org.opencontainers.image.source="https://github.com/plutzer/ProxiMate"
LABEL org.opencontainers.image.version="${PROXIMATE_VERSION}"

# Logging defaults, overridable with `docker run -e`.
ENV LOG_LEVEL=INFO
ENV PROXIMATE_LOG_DIR=/Outputs
# Subprocess output is piped; without this it would sit in a block buffer
# instead of streaming to the terminal while a stage runs.
ENV PYTHONUNBUFFERED=1

# Create a directory for output files
RUN mkdir -p /Outputs

# When the container starts, serve the GUI (3838) and the MCP tools (3839) from one
# process.  The MCP port carries no authentication: publish it on the loopback
# interface only, e.g. -p 127.0.0.1:3839:3839.
EXPOSE 3838 3839
CMD ["python3", "-u", "/GUI/server.py"]

