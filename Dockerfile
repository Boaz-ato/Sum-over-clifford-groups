FROM --platform=linux/amd64 python:3.11

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    gcc \
    g++ \
    make \
    wget \
    build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir /opt/project

COPY requirements.txt requirements.txt

RUN pip3 install --upgrade pip


RUN pip install --no-cache-dir -r requirements.txt

# Install the Python QuaEC package directly from GitHub
RUN pip install --no-cache-dir git+https://github.com/cgranade/python-quaec.git

RUN sed -i 's/from collections import Sequence/from collections.abc import Sequence/' /usr/local/lib/python3.11/site-packages/qecc/paulicollections.py

# Install GLPK from source
RUN wget http://ftp.gnu.org/gnu/glpk/glpk-5.0.tar.gz && \
    tar -xzf glpk-5.0.tar.gz && \
    cd glpk-5.0 && \
    ./configure && \
    make && \
    make install && \
    ldconfig && \
    cd .. && \
    rm -rf glpk-5.0 glpk-5.0.tar.gz

WORKDIR /opt/project
COPY . /opt/project
ENV PYTHONPATH=/opt/project
