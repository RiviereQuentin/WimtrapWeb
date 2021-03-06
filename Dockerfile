FROM bioconductor/bioconductor_docker:devel
RUN apt-get update && apt-get install -y  git-core libcurl4-openssl-dev libgit2-dev libicu-dev libssl-dev libxml2-dev make pandoc pandoc-citeproc zlib1g-dev && rm -rf /var/lib/apt/lists/*
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" >> /usr/local/lib/R/etc/Rprofile.site
RUN R -e 'install.packages("remotes")'
RUN R -e 'install.packages("utils")'
RUN R -e 'install.packages("curl")'
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install(version = "3.15")'
RUN Rscript -e 'options(repos = BiocManager::repositories());BiocManager::install("RiviereQuentin/Wimtrap", dependencies = TRUE, build_vignettes = FALSE, force = TRUE)'
RUN Rscript -e 'remotes::install_version("glue",upgrade="never", version = "1.6.1")'
RUN Rscript -e 'remotes::install_version("processx",upgrade="never", version = "3.5.2")'
RUN Rscript -e 'remotes::install_version("htmltools",upgrade="never", version = "0.5.2")'
RUN Rscript -e 'remotes::install_version("data.table",upgrade="never", version = "1.14.2")'
RUN Rscript -e 'remotes::install_version("bslib",upgrade="never", version = "0.3.1")'
RUN Rscript -e 'remotes::install_version("shiny",upgrade="never", version = "1.7.1")'
RUN Rscript -e 'remotes::install_version("config",upgrade="never", version = "0.3.1")'
RUN Rscript -e 'remotes::install_version("attempt",upgrade="never", version = "0.3.1")'
RUN Rscript -e 'remotes::install_version("testthat",upgrade="never", version = "3.1.2")'
RUN Rscript -e 'remotes::install_version("spelling",upgrade="never", version = "2.2")'
RUN Rscript -e 'remotes::install_version("thematic",upgrade="never", version = "0.1.2.1")'
RUN Rscript -e 'remotes::install_version("shinyFiles",upgrade="never", version = "0.9.1")'
RUN Rscript -e 'remotes::install_version("shinybusy",upgrade="never", version = "0.2.2")'
RUN Rscript -e 'remotes::install_version("golem",upgrade="never", version = "0.3.1")'
RUN Rscript -e 'remotes::install_version("DT",upgrade="never", version = "0.20")'
RUN mkdir /build_zone
ADD . /build_zone
WORKDIR /build_zone
RUN Rscript -e 'package.dir <- system.file(package = "Wimtrap");utils::download.file(url = "https://github.com/RiviereQuentin/carepat/archive/main.zip",destfile = paste0(package.dir, "/carepat.zip" ),timeout = 600,quiet = FALSE);utils::unzip(zipfile = paste0(package.dir, "/carepat.zip"),exdir = package.dir)'
RUN Rscript -e 'options(repos = BiocManager::repositories());BiocManager::install("RiviereQuentin/WimtrapWeb", dependencies = TRUE, build_vignettes = FALSE, force = TRUE)'
RUN rm -rf /build_zone
CMD R -e "options(shiny.host='0.0.0.0');library(WimtrapWeb);run_app()"
