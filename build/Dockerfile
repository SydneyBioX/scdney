# FROM bioconductor/bioconductor_docker:devel
FROM eddelbuettel/r2u:22.04
LABEL name='scdney'

RUN apt-get -y update
RUN apt-get -y upgrade
RUN apt-get install -y build-essential

RUN apt-get install -y curl libcurl4-openssl-dev libssl-dev

WORKDIR /install
ADD install.R .
RUN Rscript install.R

ENTRYPOINT [ "bash" ]
