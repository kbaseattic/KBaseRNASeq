FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# Install the SDK (should go away eventually)
RUN pip install --upgrade virtualenv
RUN \
  . /kb/dev_container/user-env.sh && \
  cd /kb/dev_container/modules && \
  rm -rf jars && \
  git clone https://github.com/kbase/jars && \
  rm -rf kb_sdk && \
  git clone https://github.com/kbase/kb_sdk -b develop && \
  rm -rf handle_service && \
  git clone https://github.com/kbase/handle_service && \
  cd /kb/dev_container/modules/jars && \
  make deploy && \
  cd /kb/dev_container/modules/kb_sdk && \
  make && make deploy && \
  cd /kb/dev_container/modules/handle_service && \
  make && make deploy 
RUN \
  . /kb/dev_container/user-env.sh && \
  cd /kb/dev_container/modules && \
  rm -rf data_api && \
  git clone https://github.com/kbase/data_api -b 0.3.0-dev && \
  pip install --upgrade /kb/dev_container/modules/data_api
####END OF KBASE #############################
#apt-get update && apt-get install -y ant && \
# -----------------------------------------
# Insert apt-get instructions here to install
# any required dependencies for your module.
# -----------------------------------------
RUN apt-get update && apt-get install -y unzip gcc bzip2 ncurses-dev
RUN pip install mpipe
WORKDIR /kb/module
COPY ./deps /kb/deps
RUN \
  sh /kb/deps/kb_tophat/install-tophat.sh && \
  sh /kb/deps/kb_bowtie/install-bowtie2.sh && \
  sh /kb/deps/kb_hisat2/install-hisat.sh && \
  sh /kb/deps/kb_cufflinks/install-cufflinks.sh && \
  sh /kb/deps/kb_stringTie/install-stringtie.sh && \
  sh /kb/deps/kb_tableMaker/install-tablemaker.sh

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
