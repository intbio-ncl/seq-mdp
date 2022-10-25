FROM ubuntu

RUN apt-get update && \
        DEBIAN_FRONTEND="noninteractive" apt-get install -y \
        nano \
        python3.7 \
        python3-pip \
        graphviz \
        git && \
    pip3 install \
        networkx \
        tqdm==4.31.1 \
        pandas \
        numpy \
        biopython \
        matplotlib \
        pydot

RUN git clone https://github.com/chrisAta/seq-mdp.git && \
    cd /seq-mdp && \
    git submodule update --init --recursive

RUN git clone https://github.com/chrisAta/seq-sublight.git
