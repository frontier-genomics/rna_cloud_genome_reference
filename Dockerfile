FROM python:3.12.11-bullseye

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
    samtools \
    bedtools \
    tabix \
    git \
    git-lfs \
    sudo \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install DuckDB CLI - detect architecture and download appropriate binary
RUN ARCH=$(uname -m) && \
    if [ "$ARCH" = "x86_64" ]; then \
        DUCKDB_ARCH="amd64"; \
    elif [ "$ARCH" = "aarch64" ] || [ "$ARCH" = "arm64" ]; then \
        DUCKDB_ARCH="arm64"; \
    else \
        echo "Unsupported architecture: $ARCH" && exit 1; \
    fi && \
    wget -q https://github.com/duckdb/duckdb/releases/latest/download/duckdb_cli-linux-${DUCKDB_ARCH}.zip \
    && unzip duckdb_cli-linux-${DUCKDB_ARCH}.zip \
    && mv duckdb /usr/local/bin/ \
    && chmod +x /usr/local/bin/duckdb \
    && rm duckdb_cli-linux-${DUCKDB_ARCH}.zip && \
    /usr/local/bin/duckdb --version && \
    echo "DuckDB CLI installed successfully"

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

# Set the entrypoint to the main script
ENTRYPOINT ["python"]
CMD ["-c", "import time; time.sleep(float('inf'))"]