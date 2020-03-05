# Install R version 3.5.1
FROM r-base:3.5.1

# Install Ubuntu packages
RUN apt-get -qq update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev/unstable \
    libxt-dev \
    libssl-dev \
    libmagick++-dev \
    ssh \
    curl \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=3 \
    && conda update conda \
    && conda clean --all --yes

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
RUN R -e "BiocManager::install('DECIPHER')"
RUN R -e "BiocManager::install('phyloseq')"
RUN R -e "BiocManager::install('ggtree')"
RUN R -e "BiocManager::install('phytools')"
RUN R -e "BiocManager::install('roxygen2')"
RUN R -e "BiocManager::install('rversions')"
RUN R -e "BiocManager::install('devtools')"
RUN R -e "BiocManager::install('magick')"
RUN R -e "BiocManager::install('shinydashboard')"
RUN R -e "BiocManager::install('rmarkdown')"

# Setup hyphy
ENV PATH /opt/conda/bin:$PATH
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

ENV PATH /usr/local/envs/hyphy/bin/:$PATH

# Copy configuration files into the Docker image
COPY R/shiny-server.conf  /etc/shiny-server/shiny-server.conf
RUN mkdir /srv/shiny-server/app
COPY R/app.R /srv/shiny-server/app
COPY scripts/REROOT.sh /srv/shiny-server/app
COPY scripts/likeli.bf /srv/shiny-server/app
COPY data/bestTree.nwk /srv/shiny-server/app
COPY data/HIVM.fasta /srv/shiny-server/app
COPY data/HIVM.fasta.nwk /srv/shiny-server/app
RUN ["chmod", "+w", "/srv/shiny-server/app"]
RUN ["chown","shiny:shiny","-R","/srv/shiny-server/app"]

# Make the ShinyApp available at port 80
EXPOSE 80

# Copy further configuration files into the Docker image
COPY R/shiny-server.sh /usr/bin/shiny-server.sh

RUN ["chmod", "+x", "/usr/bin/shiny-server.sh"]

CMD ["/usr/bin/shiny-server.sh"]
