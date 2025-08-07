FROM condaforge/mambaforge:latest

# For opencontainers label definitions, see:
#    https://github.com/opencontainers/image-spec/blob/master/annotations.md
LABEL org.opencontainers.image.title="MultiRTC"
LABEL org.opencontainers.image.description="ISCE3-based RTC generation for multiple sensors"
LABEL org.opencontainers.image.vendor="Forrest Williams"
LABEL org.opencontainers.image.licenses="BSD-3-Clause"

ARG DEBIAN_FRONTEND=noninteractive
ENV PYTHONDONTWRITEBYTECODE=true

RUN apt-get update && apt-get install -y --no-install-recommends unzip vim && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

ARG CONDA_UID=1000
ARG CONDA_GID=1000

RUN groupadd -g "${CONDA_GID}" --system conda && \
    useradd -l -u "${CONDA_UID}" -g "${CONDA_GID}" --system -d /home/conda -m  -s /bin/bash conda && \
    chown -R conda:conda /opt && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> /home/conda/.profile && \
    echo "conda activate base" >> /home/conda/.profile

USER ${CONDA_UID}
SHELL ["/bin/bash", "-l", "-c"]
WORKDIR /home/conda/

RUN mkdir -p ./isce3/isce3_build && \
    mkdir -p ./isce3/isce3_install && \
    git clone --branch docker https://github.com/forrestfwilliams/MultiRTC.git ./multirtc && \
    git clone --branch pfa https://github.com/forrestfwilliams/isce3.git ./isce3/isce3_src


RUN mamba env create -f ./multirtc/environment.pfa.yml && \
    conda clean -afy && \
    conda activate multirtc && \
    sed -i 's/conda activate base/conda activate multirtc/g' /home/conda/.profile

# RUN CC=clang CXX=clang++ cmake -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DCMAKE_INSTALL_PREFIX=/home/conda/isce3/isce3_install /home/conda/isce3/isce3_src make > make.tx

# ENTRYPOINT ["/hyp3-isce2/src/hyp3_isce2/etc/entrypoint.sh"]
CMD ["/bin/bash"]
