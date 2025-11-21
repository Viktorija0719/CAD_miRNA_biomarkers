FROM mambaorg/micromamba:1.5.10

WORKDIR /work

COPY environment.yml /tmp/environment.yml
RUN micromamba create -y -f /tmp/environment.yml -n snakemake_mirdeep2 && \
    micromamba clean --all --yes

SHELL ["bash", "-lc"]
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/envs/snakemake_mirdeep2/bin:$PATH

# miRDeep2 extra file
COPY extra/Rfam_for_miRDeep.fa /opt/conda/envs/snakemake_mirdeep2/Rfam_for_miRDeep.fa

# OPTIONAL: Python + R stack for downstream analysis
COPY requirements.txt /tmp/requirements.txt
COPY install.R /tmp/install.R

# Use env’s python / R (they’re on PATH)
RUN python -m pip install --no-cache-dir -r /tmp/requirements.txt || true
RUN Rscript /tmp/install.R || true

# Patch miRDeep2 sanity_check scripts
RUN perl -pi -e 's|^#!/usr/bin/env perl -W|#!/usr/bin/env perl|' \
      /opt/conda/envs/snakemake_mirdeep2/bin/sanity_check_* || true

CMD ["bash"]
