FROM hsqpqbrc/sctk-base:latest

SHELL ["/bin/bash", "-c"]

# Add the sctk seurat clustering script
ADD mast_dge.R /opt/software

# add to the PATH. Otherwise, conda can swallow return codes
ENV PATH="/opt/software/miniconda/envs/sctk_env/bin:${PATH}"

ENTRYPOINT ["/bin/bash"]