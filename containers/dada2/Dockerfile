FROM bioconductor/bioconductor_docker:RELEASE_3_18

LABEL org.opencontainers.image.source="https://github.com/UKCEH-MolecularEcology/nf-honeypi"
LABEL org.opencontainers.image.description="R 4.3 + DADA2 1.30 for nf-honeypi"
LABEL org.opencontainers.image.licenses="MIT"

RUN Rscript -e " \
    BiocManager::install('dada2', ask=FALSE, update=FALSE); \
    stopifnot(requireNamespace('dada2', quietly=TRUE)); \
    cat('dada2 version:', as.character(packageVersion('dada2')), '\n')"
