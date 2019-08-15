FROM kelshmo/covarr
LABEL maintainer='Kelsey Montgomery'
USER root
COPY . ${HOME}
RUN echo ${NB_USER}
RUN echo ${HOME}
RUN chown -R ${NB_USER} ${HOME}
USER ${NB_USER}
RUN wget https://github.com/kelshmo/nature-functional-data/raw/master/DESCRIPTION && R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-06-13/')); devtools::install_deps()"
