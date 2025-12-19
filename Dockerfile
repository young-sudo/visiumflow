FROM mambaorg/micromamba:latest
# https://github.com/mamba-org/mamba

# Create environment
COPY env.yml /tmp/env.yml

RUN micromamba create -f /tmp/env.yml -y && \
    micromamba clean --all -y
RUN micromamba run -n spatialflow pip install --no-deps "spatialdata[extra]==0.4.0"

SHELL ["micromamba", "run", "-n", "spatialflow", "/bin/bash", "-c"]
