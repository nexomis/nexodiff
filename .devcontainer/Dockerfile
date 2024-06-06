# See here for image contents: https://github.com/microsoft/vscode-dev-containers/tree/v0.231.6/containers/debian/.devcontainer/base.Dockerfile

# [Choice] Debian version (use bullseye on local arm64/Apple Silicon): bullseye, buster
ARG VARIANT="bullseye"
FROM mcr.microsoft.com/vscode/devcontainers/base:0-${VARIANT} 

RUN apt-get -y update \
  && apt-get -y install --no-install-recommends \
  software-properties-common \
  gnupg \
  tree \
  git \
  ed \ 
  less \
  locales \
  vim-tiny \
  wget \
  ca-certificates \
  fonts-texgyre \
  libc6 \
  libgcc1 \
  libstdc++6 \
  libgmp10 \
  libmpfr6 \
  libmpfr-dev \
  cmake \
  python-is-python3 \
  python3 \
  python3-pip \
  python3-yaml \
  pandoc \
  texlive-latex-base \
  texlive-fonts-extra

RUN pip3 install \
  notebook \
  h5py \
  numpy \
  scipy

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-key  \
  '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'

RUN echo "deb https://cloud.r-project.org/bin/linux/debian bullseye-cran40/" \
  > /etc/apt/sources.list.d/bullseye-cran40.list

RUN apt-get -y update \
  && apt-get -y install --no-install-recommends \
  libopenblas0-pthread \
  littler \
  r-cran-docopt \
  r-cran-littler \
  r-base \
  r-base-dev \
  r-base-core \
  r-recommended \
  r-cran-littler \
  r-cran-devtools \
  r-cran-tidyverse \
  r-cran-biocmanager \
  r-bioc-deseq2 \
  r-cran-cowplot \
  r-cran-digest \
  r-cran-dplyr \
  r-cran-factoextra \
  r-cran-furrr \
  r-cran-future \
  r-cran-ggally \
  r-cran-ggplot2 \
  r-cran-gridextra \
  r-cran-logging \
  r-cran-matrixstats \
  r-cran-openxlsx \
  r-cran-parallelly \
  r-cran-purrr \
  r-cran-r6 \
  r-cran-readr \
  r-bioc-rhdf5 \
  r-cran-rstan \
  r-cran-scales \
  r-cran-stringr \
  r-cran-tidyr \
  r-cran-tidyselect \
  r-cran-yaml \
  r-cran-pvclust \
  nodejs \
  npm

RUN Rscript -e 'install.packages("cli")' \
  && Rscript -e 'install.packages("ggh4x")' \
  && Rscript -e 'install.packages("ggmulti")' \
  && Rscript -e 'install.packages("ggrepel")' \
  && Rscript -e 'BiocManager::install("apeglm")'

# I do not know why those are required ... ( maybe in R_old ?)

RUN Rscript -e 'install.packages("tzdb")' \
 && Rscript -e 'install.packages("multcompView")' \
 && Rscript -e 'install.packages("emmeans")' \
 && Rscript -e 'install.packages("png")' \
 && Rscript -e 'install.packages("vroom")' \
 && Rscript -e 'BiocManager::install("KEGGREST")' \
 && Rscript -e 'BiocManager::install("Biostrings")' \
 && npm install -g standalone-html

RUN Rscript -e 'install.packages("languageserver")' \
 && Rscript -e 'install.packages("IRkernel")' \
 && Rscript -e 'IRkernel::installspec(user = FALSE)' \
 && Rscript -e 'install.packages("pkgload")' \
 && Rscript -e 'install.packages("tibble")'

RUN  rm -rf /tmp/downloaded_packages/ /tmp/*.rds
