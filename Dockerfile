FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.138.cpg1-1

ENV PYTHONDONTWRITEBYTECODE=1

ENV VERSION=2.2.6

WORKDIR /rdrnaseq

COPY src src/
COPY LICENSE pyproject.toml README.md ./

# pip install but don't retain the cache files
RUN pip install --no-cache-dir .
