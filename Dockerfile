FROM mambaorg/micromamba:1.5.10-noble

USER root

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    git-lfs \
    sudo \
    curl \
    wget \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
    
USER $MAMBA_USER

RUN micromamba create \
    -y -p /opt/conda/envs/bioenv \
    -c conda-forge \
    -c bioconda \
    nextflow \
    ucsc-genepredtobed==482 \
    ucsc-gtftogenepred==482 \
    samtools=1.22.1 \
    duckdb-cli=1.3.2 \
    bedtools=2.31.1 \
    seqkit=2.10.1 \
    python=3.12.12 \
    pip && \
    micromamba clean --all -y

    
# Put that env first on PATH so binaries are available without activation
ENV PATH=/opt/conda/envs/bioenv/bin:$PATH

RUN nextflow -version
    
# Keeps Python from generating .pyc files in the container
ENV PYTHONDONTWRITEBYTECODE=1

# Turns off buffering for easier container logging
ENV PYTHONUNBUFFERED=1

# Set the working directory
WORKDIR /app

# Copy the requirements file and install Python dependencies
COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

# Copy the application code
COPY . .

# Set the entrypoint to the main script
ENTRYPOINT ["python"]
CMD ["-c", "import time; time.sleep(float('inf'))"]