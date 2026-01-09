
FROM mcr.microsoft.com/devcontainers/base:noble

ENV R_VERSION=4.4.2
ENV BIOC_VERSION=3.20

# Install a specific version of using posit release of gdebi packages.
RUN export DEBIAN_FRONTEND=noninteractive \ 
  && apt-get update \
  && apt-get -y install --no-install-recommends \
    gdebi-core curl ca-certificates apt-utils apt-file \
  && . /etc/os-release \
  && curl -O https://cdn.rstudio.com/r/${ID}-$(echo $VERSION_ID | sed 's/\.//g' )/pkgs/r-${R_VERSION}_1_amd64.deb \
  && gdebi --n r-${R_VERSION}_1_amd64.deb \
  && apt-get -y clean \
  && apt-get -y purge gdebi \
  && apt-get -y autoremove \
  && rm -rf /var/lib/apt/lists/* \
  && rm r-${R_VERSION}_1_amd64.deb \
  && ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R \
  && ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript \
  && R_HOME=$(Rscript -e 'R.home()' | sed -e 's/\[1\] //g' | sed -e 's/"//g') \
  && R_PROFILE=$R_HOME/etc/Rprofile.site \
  && echo "options(BioC_mirror = \"https://packagemanager.posit.co/bioconductor\")" >> $R_PROFILE \
  && echo "options(BIOCONDUCTOR_CONFIG_FILE = \"https://packagemanager.posit.co/bioconductor/config.yaml\")" >> $R_PROFILE \
  && echo "options(repos = c(CRAN = \"https://packagemanager.posit.co/cran/__linux__/${VERSION_CODENAME}/latest\"))" >> $R_PROFILE \
  && Rscript -e 'install.packages(c("BiocManager", "rspm"), repos = "https://cloud.r-project.org", clean = TRUE)' \
  && mkdir -p $R_HOME/site-library \
  && echo ".library <- \"${R_HOME}/site-library\"" >> $R_PROFILE

RUN export DEBIAN_FRONTEND=noninteractive \
  && apt-get update && apt-get install -y \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  && wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.7.30/quarto-1.7.30-linux-amd64.deb \
  && dpkg -i quarto-1.7.30-linux-amd64.deb  \
  && rm -rf /var/lib/apt/lists/*

# Install nvm, Node.js v24, and allure for all users
RUN export DEBIAN_FRONTEND=noninteractive \
  && curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.3/install.sh | bash \
  && export NVM_DIR="$HOME/.nvm" \
  && [ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh" \
  && nvm install 24 \
  && nvm alias default 24 \
  && npm install -g allure \
  && echo 'export NVM_DIR="$HOME/.nvm"' >> /etc/bash.bashrc \
  && echo '[ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh"' >> /etc/bash.bashrc \
  && echo '[ -s "$NVM_DIR/bash_completion" ] && \. "$NVM_DIR/bash_completion"' >> /etc/bash.bashrc

COPY install.r requirements.txt requirements.dev.txt /root/

RUN Rscript /root/install.r /root/requirements.txt \
  && Rscript /root/install.r /root/requirements.dev.txt