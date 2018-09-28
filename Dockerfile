# image: afioregartland/rnaseq_prep
FROM ubuntu:18.04
# FROM amazonlinux:latest
LABEL maintainer="agartlan@fredhutch.org"

ENV PACKAGES git gcc make g++ cmake libboost-all-dev liblzma-dev libbz2-dev \
    ca-certificates zlib1g-dev curl unzip autoconf trimmomatic default-jre gnupg \
    ed less locales vim-tiny nano wget fonts-texgyre python3.6 python-pip python-dev build-essential

ENV SALMON_VERSION 0.11.3
ENV R_BASE_VERSION 3.5.1
ENV DEBIAN_FRONTEND noninteractive

WORKDIR /home

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen
RUN locale-gen en_US.utf8
RUN /usr/sbin/update-locale LANG=en_US.UTF-8

## Now install R and littler, and create a link for littler in /usr/local/bin
## Also set a default CRAN repo, and make sure littler knows about it too
## Also install stringr to make dococt install (from source) easier
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 8B48AD6246925553
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 7638D0442B90D010
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 04EE7237B7D453EC

## Use Debian unstable via pinning -- new style via APT::Default-Release
RUN echo "deb http://http.debian.net/debian sid main" > /etc/apt/sources.list.d/debian-unstable.list \
        && echo 'APT::Default-Release "testing";' > /etc/apt/apt.conf.d/default 

RUN apt-get update && \
    apt-get install -t unstable -y --no-install-recommends \
        littler \
        r-cran-littler \
        r-cran-stringr \
        r-base=${R_BASE_VERSION}-* \
        r-base-dev=${R_BASE_VERSION}-* \
        r-recommended=${R_BASE_VERSION}-*
RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"))' >> /etc/R/Rprofile.site
RUN echo 'source("/etc/R/Rprofile.site")' >> /etc/littler.r
RUN ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
    && ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
    && ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
    && ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r
RUN install.r docopt \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
    && rm -rf /var/lib/apt/lists/*

RUN echo 'local({r <- getOption("repos"); r["CRAN"] <- "http://cran.r-project.org"; options(repos=r)})' > ~/.Rprofile
RUN R -e 'source("http://bioconductor.org/biocLite.R"); biocLite("DESeq2"); biocLite("tximport"); biocLite("readr");'

# Install salmon from source
RUN curl -k -L https://github.com/COMBINE-lab/salmon/archive/v${SALMON_VERSION}.tar.gz -o salmon-v${SALMON_VERSION}.tar.gz && \
    tar xzf salmon-v${SALMON_VERSION}.tar.gz && \
    cd salmon-${SALMON_VERSION} && \
    mkdir build && \
    cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local && make && make install

# Salmon alternative?
# ADD https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/Salmon-${SALMON_VERSION}_linux_x86_64.tar.gz /opt/Salmon-${SALMON_VERSION}_linux_x86_64.tar.gz
# RUN cd /opt && tar -zxvf Salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && cp -p /opt/Salmon-*/bin/salmon /usr/local/bin && cp -p /opt/Salmon-*/lib/* /usr/local/lib && cd /opt && rm -rf Salmon*

# For dev versions
#RUN git clone https://github.com/COMBINE-lab/salmon.git && \
#    cd salmon && \
#    git checkout develop && \
#    mkdir build && \
#    cd build && \
#    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local && make && make install

# Install fastqc
RUN curl -O https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip && \
    unzip fastqc_v0.11.7.zip && \
    rm fastqc_v0.11.7.zip && \
    chmod 755 FastQC/fastqc && \
    ln -s /FastQC/fastqc /bin/fastqc

RUN curl -o /usr/local/bin/faToTwoBit http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit /usr/local/bin/faToTwoBit
RUN chmod 755 /usr/local/bin/faToTwoBit

RUN curl -o /usr/local/bin/twoBitToFa http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa /usr/local/bin/twoBitToFa
RUN chmod 755 /usr/local/bin/twoBitToFa

RUN pip install --upgrade pip
RUN pip install pybedtools pysam biopython numpy pandas scipy matplotlib sckitbio jupyter

ENV PATH /home/salmon-${SALMON_VERSION}/bin:${PATH}
