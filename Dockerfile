FROM bioconductor/bioconductor_docker:RELEASE_3_22

RUN R -e "BiocManager::install(
    c('tximport', 'DESeq2'),
    ask = FALSE,
    update = FALSE
)"

RUN R -e "install.packages(
    'ggplot2',
    repos = 'https://cloud.r-project.org'
)"

RUN R -e "library(tximport);
          library(DESeq2);
          library(ggplot2)"
