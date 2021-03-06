# image: afioregartland/rnaseq_prep
FROM ubuntu:18.04
# FROM amazonlinux:latest
LABEL maintainer="agartlan@fredhutch.org"

ENV PACKAGES git gcc make g++ cmake libboost-all-dev liblzma-dev libbz2-dev \
    ca-certificates zlib1g-dev curl unzip autoconf trimmomatic default-jre gnupg \
    ed less locales vim-tiny nano wget fonts-texgyre python3.6 python3.6-dev build-essential \
    openjdk-8-jdk wget screen

ENV SALMON_VERSION 0.13.1
ENV R_BASE_VERSION 3.6.0
ENV DEBIAN_FRONTEND noninteractive

WORKDIR /home

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN wget https://bootstrap.pypa.io/get-pip.py \
    && python3.6 get-pip.py
RUN ln -s /usr/bin/python3.6 /usr/local/bin/python3

RUN yes w | pip3 install --upgrade pip
RUN yes w | pip3 install setuptools numpy
RUN yes w | pip3 install 'matplotlib<3.0.0,>=2.1.1'
RUN yes w | pip3 install pybedtools pysam biopython pandas scipy scikit-bio jupyter feather-format awscli boto3 botocore multiqc

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
RUN R -e 'source("http://bioconductor.org/biocLite.R"); biocLite("DESeq2"); biocLite("tximport"); biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene"); biocLite("readr"); install.packages("feather", repos="http://cran.r-project.org");'

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
RUN curl -O https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip && \
    unzip fastqc_v0.11.8.zip && \
    rm fastqc_v0.11.8.zip && \
    chmod 755 FastQC/fastqc && \
    ln -s /home/FastQC/fastqc /bin/fastqc

RUN curl -o /usr/local/bin/faToTwoBit http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
RUN chmod 755 /usr/local/bin/faToTwoBit

RUN curl -o /usr/local/bin/twoBitToFa http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
RUN chmod 755 /usr/local/bin/twoBitToFa

RUN curl -k -L https://raw.githubusercontent.com/FredHutch/url-fetch-and-run/master/fetch-and-run/fetch_and_run.sh -o /usr/local/bin/fetch_and_run.sh
RUN chmod +x /usr/local/bin/fetch_and_run.sh

# Download & extract STAR - Repo includes binaries for linux
RUN curl -sSL https://github.com/alexdobin/STAR/blob/master/bin/Linux_x86_64_static/STAR?raw=true -o /usr/local/bin/star
RUN chmod +x /usr/local/bin/star

# RNA-seQC
RUN curl -sSL https://github.com/broadinstitute/rnaseqc/releases/download/v2.1.0/rnaseqc.v2.1.0.linux.gz -o /usr/local/bin/rnaseqc
RUN chmod +x /usr/local/bin/rnaseqc

# Picard tools
RUN mkdir /home/picard-tools && \
    wget --no-check-certificate -P /home/picard-tools/ https://github.com/broadinstitute/picard/releases/download/2.20.0/picard.jar
    
# SAM tools
RUN mkdir ~/src && \
    cd ~/src && \
    git clone https://github.com/samtools/htslib && \
    git clone https://github.com/samtools/samtools && \
    cd samtools && \
    make && \
    cp samtools /usr/bin

ENV PATH /home/salmon-${SALMON_VERSION}/bin:/home/picard-tools:/usr/bin:$PATH

RUN cd /home
# ENTRYPOINT ["/bin/bash"]
