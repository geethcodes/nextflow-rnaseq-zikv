FROM bioconductor/bioconductor_docker:RELEASE_3_22

RUN R -e "BiocManager::install(c('tximport','DESeq2','apeglm'), ask=FALSE, update=FALSE)"
RUN R -e "install.packages(c('ggplot2','pheatmap'), repos='https://cloud.r-project.org')"

# Build fails here if anything didn't install
RUN R -e "library(tximport); library(DESeq2); library(apeglm); library(ggplot2); library(pheatmap); cat('R packages OK\n')"
