FROM bioconductor/devel_proteomics
MAINTAINER Christian Panse <cp@fgcz.ethz.ch>

RUN sudo apt-get install git -y 

RUN echo "install.packages(c('testthat', 'knitr'))" | R --no-save

RUN cd /tmp && git clone https://github.com/protViz/msqc1.git \
        && R CMD build msqc1 \
	&& R CMD INSTALL msqc1_*.tar.gz \
	&& echo "library(msqc1); packageVersion('msqc1')" \
	| R --no-save 
