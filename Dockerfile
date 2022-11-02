FROM intel/oneapi-hpckit
RUN apt-get update && apt-get install -y hdf5-tools hdf5-helpers libhdf5-dev
RUN mkdir /opt/download /work
WORKDIR /opt/download
RUN git clone --depth 1 https://github.com/avr-aics-riken/TextParser \
    && mkdir TextParser/build \
    && cd TextParser/build \
    && cmake -DINSTALL_DIR=/opt/TextParser .. \
    && make -j4 install \
    && rm -r /opt/download/TextParser

WORKDIR /work