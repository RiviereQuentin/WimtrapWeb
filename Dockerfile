FROM rocker/r-ver:4.1.2
RUN apt-get update && apt-get install -y  git-core libcurl4-openssl-dev libgit2-dev libicu-dev libssl-dev libxml2-dev make pandoc pandoc-citeproc zlib1g-dev && rm -rf /var/lib/apt/lists/*
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" >> /usr/local/lib/R/etc/Rprofile.site
RUN R -e 'install.packages("remotes")'
RUN R -e 'install.packages("utils")'
RUN R -e 'install.packages("curl")'
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install(version = "3.14")'
RUN Rscript -e 'remotes::install_github("RiviereQuentin/Wimtrap@c8204ad804b1eeb41cc1c3a41fc2b7eebe91722d")'
RUN Rscript -e 'BiocManager::install(pkgs = "GenomicRanges", update = FALSE,  ask = FALSE)'
RUN Rscript -e 'remotes::install_version("glue",upgrade="never", version = "1.6.1")'
RUN Rscript -e 'remotes::install_version("processx",upgrade="never", version = "3.5.2")'
RUN Rscript -e 'remotes::install_version("htmltools",upgrade="never", version = "0.5.2")'
RUN Rscript -e 'remotes::install_version("data.table",upgrade="never", version = "1.14.2")'
RUN Rscript -e 'remotes::install_version("bslib",upgrade="never", version = "0.3.1")'
RUN Rscript -e 'remotes::install_version("pROC",upgrade="never", version = "1.18.0")'
RUN Rscript -e 'remotes::install_version("shiny",upgrade="never", version = "1.7.1")'
RUN Rscript -e 'remotes::install_version("config",upgrade="never", version = "0.3.1")'
RUN Rscript -e 'remotes::install_version("attempt",upgrade="never", version = "0.3.1")'
RUN Rscript -e 'remotes::install_version("testthat",upgrade="never", version = "3.1.2")'
RUN Rscript -e 'remotes::install_version("spelling",upgrade="never", version = "2.2")'
RUN Rscript -e 'remotes::install_version("thematic",upgrade="never", version = "0.1.2.1")'
RUN Rscript -e 'remotes::install_version("xgboost",upgrade="never", version = "1.5.0.2")'
RUN Rscript -e 'remotes::install_version("caret",upgrade="never", version = "6.0-90")'
RUN Rscript -e 'remotes::install_version("shinyFiles",upgrade="never", version = "0.9.1")'
RUN Rscript -e 'remotes::install_version("shinybusy",upgrade="never", version = "0.2.2")'
RUN Rscript -e 'remotes::install_version("golem",upgrade="never", version = "0.3.1")'
RUN Rscript -e 'remotes::install_version("DT",upgrade="never", version = "0.20")'
RUN mkdir /build_zone
ADD . /build_zone
WORKDIR /build_zone
RUN R -e 'remotes::install_local(upgrade="never")'
RUN rm -rf /build_zone
CMD R -e "options('shiny.port'=$PORT,shiny.host='0.0.0.0');PlantCarepat::run_app()"