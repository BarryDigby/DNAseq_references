FROM nfcore/base:1.9
LABEL authors="Barry Digby" \
      description="Docker image containing all software requirements for prepping Germline VC"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml python=2.7.15 && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/DNAseq_references/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name DNAseq_references > DNAseq_references.yml
