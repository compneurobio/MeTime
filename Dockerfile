FROM rocker/r-ver:4.3.3

LABEL org.opencontainers.image.title="MeTime"
LABEL org.opencontainers.image.description="Container image for the MeTime R package"
LABEL org.opencontainers.image.source="https://github.com/compneurobio/MeTime"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    libxt-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libglpk40 \
    libgit2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libjpeg-dev \
    libtiff5-dev \
    libwebp-dev \
    git \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

RUN R -q -e 'install.packages(c("BiocManager", "remotes"), repos = "https://cloud.r-project.org")'

WORKDIR /opt/MeTime
COPY . /opt/MeTime

RUN R -q -e 'options(repos = BiocManager::repositories()); remotes::install_deps(dependencies = TRUE)' \
 && R CMD INSTALL --no-manual --no-build-vignettes .

CMD ["R", "-q", "-e", "library(MeTime); data('humet_object'); print(packageVersion('MeTime'))"]
