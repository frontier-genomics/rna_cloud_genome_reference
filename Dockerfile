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