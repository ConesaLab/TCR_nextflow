FROM r-base:4.1.2
LABEL maintainer=" Susanna Marquez [smargar1@upv.edu.es]" \
    description="Environment and dependencies for TCR Workflow analysis"

# Install ps, for Nextflow. https://www.nextflow.io/docs/latest/tracing.html
RUN apt-get update && \
    apt-get install -y procps \
    pandoc \
    libcurl4-openssl-dev \
    r-cran-curl

# Install required R packages
ARG R_DEPS="c( \
    'BiocManager', \
    'bookdown', \
    'dplyr', \
    'ggplot2', \
    'ggpubr', \
    'grid', \
    'gridExtra', \
    'immunarch', \
    'microbiome', \
    'openxlsx', \
    'plyr', \
    'purrr', \
    'RcmdrMisc', \
    'RColorBrewer', \
    'reshape2', \
    'stringdist', \
    'tibble', \
    'tidyr', \
    'vegan', \
    'viridis' \
    )"
ARG R_BIOC_DEPS="c( \
    'Biobase', \
    'microbiome')"

RUN Rscript -e "install.packages(${R_DEPS}, clean=TRUE)" && \
    Rscript -e "BiocManager::install(${R_BIOC_DEPS})"  && \
    Rscript -e "install.packages('NMF', clean=TRUE)"
CMD ["R"]
