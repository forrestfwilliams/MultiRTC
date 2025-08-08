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
    git clone --branch pfa https://github.com/forrestfwilliams/isce3.git ./isce3/isce3_src

COPY --chown=${CONDA_UID}:${CONDA_GID} . ./multirtc/

RUN mamba env create -f ./multirtc/environment.pfa.yml && \
    conda clean -afy && \
    conda activate multirtc && \
    sed -i 's/conda activate base/conda activate multirtc/g' /home/conda/.profile && \
    python -m pip install --no-cache-dir ./multirtc && \
    conda remove --force -y isce3

RUN cd ./isce3/isce3_build && \
    cmake -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DCMAKE_INSTALL_PREFIX=/home/conda/isce3/isce3_install /home/conda/isce3/isce3_src && \
    make -j"$(nproc)" && \
    make install && \
    cd ../..

RUN rm -rf ./isce3/isce3_build && rm -rf ./isce3/isce3_src

ENV ISCE_INSTALL=/home/conda/isce3/isce3_install \
    PATH=$ISCE_INSTALL/bin:$ISCE_INSTALL/packages/nisar/workflows:$PATH \
    PYTHONPATH=$ISCE_INSTALL/packages:$ISCE_INSTALL/lib:$PYTHONPATH \
    LD_LIBRARY_PATH=$ISCE_INSTALL/lib64:$LD_LIBRARY_PATH \
    DYLD_LIBRARY_PATH=$ISCE_INSTALL/lib:$DYLD_LIBRARY_PATH \
    GDAL_VRT_ENABLE_PYTHON=YES

ENTRYPOINT ["/home/conda/multirtc/src/multirtc/etc/entrypoint.sh"]
CMD ["-h"]
