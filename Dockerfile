FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.137.cpg1-2

ENV PYTHONDONTWRITEBYTECODE=1

ENV VERSION=2.0.0

WORKDIR /rdrnaseq

COPY src src/
COPY LICENSE pyproject.toml README.md ./

# pip install but don't retain the cache files
RUN pip install --no-cache-dir .
