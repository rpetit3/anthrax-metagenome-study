FROM rpetit3/nextconda-base
MAINTAINER robbie.petit@gmail.com

# Install dependencies and via Bioconda
# Update Nextflow and Conda
RUN apt-get -qq update \
    && apt-get -qq -y --no-install-recommends install \
        gcc \
        g++ \
        less \
        libicu-dev \
        libsqlite3-dev \
        libxml2-dev \
        make \
        python3-tk \
        sqlite3 \
        zlib1g-dev \
    && nextflow self-update \
    && conda upgrade conda \
    && conda install -y art=2016.06.05 \
    && conda install -y biopython=1.70 \
    && conda install -y blast=2.7.1 \
    && conda install -y jellyfish=2.2.6 \
    && conda install -y mash=2.0 \
    && conda install -y sra-tools=2.8.2 \
    && conda install -y perl-bioperl=1.6.924 \
    && conda install -y  quicktree=2.2 \
    && conda clean --all --yes \
    && pip install --upgrade pip setuptools \
    && pip install numpy \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log /tmp/* /var/tmp/* \
    && mkdir -p /data /tmp/install

# Final Programs
# Mashtree
RUN cd /tmp/install \
    && curl -sSL https://github.com/lskatz/mashtree/archive/v0.32.tar.gz -o mashtree.tar.gz \
    && tar -xzf mashtree.tar.gz \
    && cd mashtree-0.32 \
    && perl Makefile.PL \
    && make install \
    && ln -s /usr/local/bin/mashtree /usr/local/bin/mashtree.pl \
    && rm -rf /tmp/install

# Final touches
COPY data /opt/data
COPY results /opt/results
COPY scripts /tmp/scripts
RUN chmod 755 /tmp/scripts/* \
    && mv /tmp/scripts/* /usr/local/bin \
    && rm -rf /tmp/*

WORKDIR /data
