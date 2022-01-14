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
ARG R_DEPS="c('bookdown', \
    'devtools', \
    'dplyr', \
    'ggplot2', \
    'ggpubr', \
    'grid', \
    'gridExtra', \
    'immunarch', \
    'microbiome', \
    'NMF', \
    'openxlsx', \
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
ARG R_BIOC_DEPS="c('microbiome')"

RUN Rscript -e "install.packages(${R_DEPS}, clean=TRUE)" && \
    Rscript -e "BiocManager::install(${R_BIOC_DEPS})"

CMD ["R"]
