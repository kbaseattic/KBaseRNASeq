FROM kbase/kbase:sdkbase2.latest
MAINTAINER KBase Developer
# -----------------------------------------
# Insert apt-get instructions here to install
# any required dependencies for your module.
# -----------------------------------------
RUN apt-get update && apt-get install -y unzip gcc bzip2 ncurses-dev sysstat
RUN pip install mpipe

WORKDIR /kb/module
COPY ./deps /kb/deps
RUN \
  sh /kb/deps/kb_tophat/install-tophat.sh && \
  sh /kb/deps/kb_bowtie/install-bowtie2.sh && \
  sh /kb/deps/kb_hisat2/install-hisat.sh && \
  sh /kb/deps/kb_cufflinks/install-cufflinks.sh && \
  sh /kb/deps/kb_stringTie/install-stringtie.sh && \
  sh /kb/deps/kb_tableMaker/install-tablemaker.sh && \
  sh /kb/deps/kb_ballgown/install-ballgown.sh

RUN \
  sh /kb/deps/pylib.sh

COPY ./ /kb/module
RUN \
  . /kb/dev_container/user-env.sh && \
  cd /kb/module && \
  kb-sdk install AssemblyUtil && \
  kb-sdk install GenomeFileUtil
RUN \
  cd /kb/module && \
  make && make deploy
#  rm -rf narrative_method_store \
ENV PATH=$PATH:/kb/dev_container/modules/kb_sdk/bin
WORKDIR /kb/module
RUN mkdir -p /kb/module/work
RUN pip install --upgrade requests==2.7.0
RUN pip freeze | grep requests
ENTRYPOINT [ "./scripts/entrypoint.sh" ]
CMD [ ]
