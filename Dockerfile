FROM r-base:4.1.2
LABEL maintainer=" Susanna Marquez [smargar1@upv.edu.es]" \
    description="Environment and dependencies for TCR Workflow analysis"

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
RUN Rscript -e "install.packages(${R_DEPS}, clean=TRUE)"
