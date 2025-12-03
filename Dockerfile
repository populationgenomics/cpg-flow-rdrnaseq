FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.134.cpg2-1

ENV PYTHONDONTWRITEBYTECODE=1

WORKDIR /rdrnaseq

COPY src src/
COPY LICENSE pyproject.toml README.md ./

# pip install but don't retain the cache files
RUN pip install --no-cache-dir .
RUN micromamba install -y -c bioconda -c conda-forge \
    # 2. List all required bioinformatics tools
    star samtools sambamba fastp cutadapt \
    outrider fraser \
    picard \
    # 3. List general utilities (optional but recommended)
    less wget curl \
    && micromamba clean --all --force-pkgs