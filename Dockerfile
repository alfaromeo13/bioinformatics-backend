FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt update -qq && apt install -y \
    python3 python3-pip python3-dev \
    build-essential gfortran libgfortran5 \
    libstdc++6 libc6 libx11-6 libxext6 libxrender1 libxt6 libxi6 \
    openmpi-bin libopenmpi-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir --break-system-packages \
    numpy biopython matplotlib scipy pandas seaborn

COPY charmm_program /usr/local/charmm_program

ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

WORKDIR /IIME
CMD ["python3", "IIME.py", "--help"]
