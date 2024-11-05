# shiny-server / CD 26062024
# Uses bioconductor image as base
FROM docker.io/bioconductor/bioconductor_docker:RELEASE_3_19

# Set noninteractive to avoid prompts during build
ENV DEBIAN_FRONTEND noninteractive

# Update system and install packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends apt-utils gdebi-core && \
    apt-get full-upgrade -y && \
    curl "https://download3.rstudio.org/ubuntu-18.04/x86_64/shiny-server-$(curl https://download3.rstudio.org/ubuntu-18.04/x86_64/VERSION)-amd64.deb" -o ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    apt-get clean && \
    rm -f ss-latest.deb && \
    rm -rf /var/lib/apt/lists/* /etc/services.d/rstudio /etc/init.d/rstudio-server /etc/cont-init.d

# Install rrvgo package and all its dependencies
RUN R --slave -e 'install.packages(c("devtools", "shiny", "shinydashboard"))' && \
    R --slave -e 'BiocManager::install(c("GOSemSim", "AnnotationDbi", "GO.db", "pheatmap", "ggplot2", "ggrepel", "treemap", "tm", "wordcloud", "shiny", "grDevices", "grid", "stats", "methods", "umap", "knitr", "rmarkdown", "BiocStyle", "DT", "plotly", "heatmaply", "magrittr", "utils", "clusterProfiler", "DOSE", "slam"))' && \
    R --slave -e 'BiocManager::install(c("org.Ag.eg.db", "org.At.tair.db", "org.Bt.eg.db", "org.Ce.eg.db", "org.Cf.eg.db", "org.Dm.eg.db", "org.Dr.eg.db", "org.EcK12.eg.db", "org.EcSakai.eg.db", "org.Gg.eg.db", "org.Hs.eg.db", "org.Mm.eg.db", "org.Mmu.eg.db", "org.Pt.eg.db", "org.Rn.eg.db", "org.Sc.sgd.db", "org.Ss.eg.db", "org.Xl.eg.db"))' && \
    R --slave -e 'BiocManager::install("rrvgo")' && \
    R --slave -e 'file.copy(system.file("shiny_rrvgo", package="rrvgo"), "/srv/shiny-server/", recursive=TRUE)' && \
    sed -i 's|#{{Impressum-placeholder}}|, p(style="text-align: right;", a("Institute of Molecular Biology gGmbH", href="https://imb.de/", target="_blank")), p(style="text-align: right;", a("Impressum - Imprint", href="https://imb.de/impressum-imprint", target="_blank"))|' /srv/shiny-server/shiny_rrvgo/app.R

# Setup permissions
RUN chown -R shiny:shiny /var/lib/shiny-server && \
    sed -i 's/directory_index on/directory_index off/' /etc/shiny-server/shiny-server.conf && \
    rm -rf /srv/shiny-server/index.html /srv/shiny-server/sample-apps

# Set user and run shiny-server
USER shiny
CMD ["shiny-server"]

