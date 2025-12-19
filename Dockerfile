FROM mambaorg/micromamba:latest
# https://github.com/mamba-org/mamba

# Create environment
COPY env.yml /tmp/env.yml

RUN micromamba create -f /tmp/env.yml -y -v && \
    micromamba clean --all --yes

SHELL ["micromamba", "run", "-n", "spatialflow", "/bin/bash", "-c"]
