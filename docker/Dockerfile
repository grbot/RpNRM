# Install R version 3.6.3
FROM r-base:4.2.2

# Install Ubuntu packages
ARG DEBIAN_FRONTEND=noninteractive
RUN rm -rf /var/lib/apt/lists/* && apt-get update --fix-missing && apt-get upgrade -y && apt-get install -y  \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libmagick++-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    ssh \
    curl \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=3 \
    && conda update conda \
    && conda clean --all --yes

# There is an issue with compatibility with new OpenSSL library, need to access the older version in buildinf devtools
RUN wget http://nz2.archive.ubuntu.com/ubuntu/pool/main/o/openssl/libssl1.1_1.1.1f-1ubuntu2.16_amd64.deb
RUN dpkg -i libssl1.1_1.1.1f-1ubuntu2.16_amd64.deb

RUN wget --no-verbose https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb

RUN apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

# Install R packages that are required
RUN R -e "install.packages(c('shiny', 'shinydashboard'), repos='http://cran.rstudio.com/')"

RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('rmarkdown')"
RUN R -e "BiocManager::install('DECIPHER')"
RUN R -e "BiocManager::install('phyloseq')"
RUN R -e "BiocManager::install('ggtree')"
RUN R -e "BiocManager::install('phytools')"
RUN R -e "BiocManager::install('phangorn')"
RUN R -e "BiocManager::install('roxygen2')"
RUN R -e "BiocManager::install('rversions')"
RUN R -e "BiocManager::install('devtools')"
RUN R -e "BiocManager::install('magick')"
RUN R -e "BiocManager::install('shinydashboard')"
RUN R -e "BiocManager::install('ape')"
RUN R -e "BiocManager::install('magrittr')"
RUN R -e "BiocManager::install('ggplot2')"
RUN R -e "BiocManager::install('Biostrings')"
RUN R -e "BiocManager::install('dplyr')"
RUN R -e "BiocManager::install('readr')"
RUN R -e "BiocManager::install('shinyjs')"
RUN R -e "BiocManager::install('caper')"
RUN R -e "BiocManager::install('geiger')"
RUN R -e "BiocManager::install('DT')"
RUN R -e "BiocManager::install('dashboardthemes')"

# Setup hyphy
ENV PATH /opt/conda/bin:$PATH
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

ENV PATH /usr/local/envs/hyphy/bin/:$PATH

# Copy configuration files into the Docker image
COPY R/shiny-server.conf  /etc/shiny-server/
RUN mkdir -p /srv/shiny-server/rpnrm
#RUN mkdir  -p /srv/shiny-server/rpnrm/
RUN mkdir -p /srv/shiny-server/rpnrm/data
RUN mkdir -p /srv/shiny-server/rpnrm/scripts
COPY R/rroot.R /srv/shiny-server/rpnrm/
COPY R/app.R /srv/shiny-server/rpnrm/
COPY scripts/aic.sh /srv/shiny-server/rpnrm/
COPY scripts/modeltest.bf /srv/shiny-server/rpnrm/
COPY scripts/REROOT.sh /srv/shiny-server/rpnrm/
COPY scripts/likeli.bf /srv/shiny-server/rpnrm/
#COPY data/seq.fasta /srv/shiny-server/rpnrm/data
#COPY data/seq.nwk /srv/shiny-server/rpnrm/data

# Make the ShinyApp available at port 80
EXPOSE 3838

RUN chown -R shiny:shiny /srv/shiny-server

# Copy further configuration files into the Docker image
COPY R/shiny-server.sh /usr/bin/shiny-server.sh
ENV SHINY_LOG_STDERR=1
RUN ["chmod", "+x", "/usr/bin/shiny-server.sh"]

RUN ["chmod", "+w", "/srv/shiny-server/rpnrm"]
RUN ["chown","shiny:shiny","-R","/srv/shiny-server/rpnrm"]
RUN ["chown","shiny:shiny","-R","/var/lib/shiny-server"]

CMD ["/usr/bin/shiny-server.sh"]
