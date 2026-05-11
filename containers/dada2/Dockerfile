FROM rocker/r-ver:4.3.3

LABEL org.opencontainers.image.source="https://github.com/UKCEH-MolecularEcology/nf-honeypi"
LABEL org.opencontainers.image.description="R 4.3.3 + DADA2 for nf-honeypi"
LABEL org.opencontainers.image.licenses="MIT"

RUN Rscript -e " \
    install.packages('BiocManager', repos='https://cloud.r-project.org', quiet=TRUE); \
    BiocManager::install('dada2', ask=FALSE, update=FALSE, version='3.18', quiet=TRUE); \
    if (!requireNamespace('dada2', quietly=TRUE)) stop('dada2 failed to install'); \
    cat('dada2 version:', as.character(packageVersion('dada2')), '\n')"
