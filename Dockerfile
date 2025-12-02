FROM python:3.12.11-bullseye

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
    git \
    git-lfs \
    sudo \
    openjdk-17-jre \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /tmp

ARG MAMBA_VERSION=2.3.3

# Install micromamba (arch-specific)
RUN ARCH=$(uname -m) && \
    if [ "$ARCH" = "x86_64" ]; then \
        MAMBA_URL="https://micro.mamba.pm/api/micromamba/linux-64/${MAMBA_VERSION}"; \
    elif [ "$ARCH" = "aarch64" ]; then \
        MAMBA_URL="https://micro.mamba.pm/api/micromamba/linux-aarch64/${MAMBA_VERSION}"; \
    else \
        echo "Unsupported architecture: $ARCH" && exit 1; \
    fi && \
    curl -L "$MAMBA_URL" | tar -xvj bin/micromamba && \
    mv bin/micromamba /usr/local/bin/micromamba && \
    rm -rf bin

# 🔸 Force micromamba to use plain repodata.json instead of .zst
RUN mkdir -p /etc/mamba && \
  printf '%s\n' \
    "channels:" \
    "  - conda-forge" \
    "  - bioconda" \
    "repodata_fns:" \
    "  - repodata.json" \
    > /etc/mamba/.mambarc

# Create an env in a fixed path and install samtools + friends
RUN /usr/local/bin/micromamba clean -a && \
    /usr/local/bin/micromamba create -y \
      -p /opt/conda/envs/bioenv \
      -c conda-forge -c bioconda \
      ucsc-genepredtobed==482 \
      ucsc-gtftogenepred==482 \
      samtools=1.22.1 \
      duckdb-cli=1.3.2 \
      bedtools=2.31.1 \
      seqkit=2.10.1 && \
    /usr/local/bin/micromamba clean --all -y

# Put that env first on PATH so binaries are available without activation
ENV PATH=/opt/conda/envs/bioenv/bin:$PATH

# Keeps Python from generating .pyc files in the container
ENV PYTHONDONTWRITEBYTECODE=1

# Turns off buffering for easier container logging
ENV PYTHONUNBUFFERED=1

# Create a non-root user for development
ARG USERNAME=devuser
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# Set the working directory
WORKDIR /app

# Copy the requirements file and install Python dependencies
COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

# Copy the application code
COPY . .

RUN chown -R $USERNAME:$USERNAME /app

# Install nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod go+rx /usr/local/bin/nextflow && \
    nextflow -version && \
    echo "Nextflow installed successfully"

# Set the entrypoint to the main script
ENTRYPOINT ["python"]
CMD ["-c", "import time; time.sleep(float('inf'))"]