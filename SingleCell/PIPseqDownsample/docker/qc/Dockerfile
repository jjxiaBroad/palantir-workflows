# Dockerfile for the PIPseq batch downsample-and-combine pipeline
FROM continuumio/miniconda3:24.1.2-0

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Create conda environment with Python packages
RUN conda install -y -c conda-forge -c bioconda \
    python=3.10 \
    pandas \
    numpy \
    scipy \
    h5py \
    anndata \
    scanpy \
    && conda clean -a -y

# Copy requirements.txt and install Python packages
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Set PATH for scripts
ENV PATH="/app/bin:${PATH}"

# Give all users full access to /app
RUN chmod -R 777 /app

CMD ["/bin/bash"]
