FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# Install the SDK (should go away eventually)
RUN \
  . /kb/dev_container/user-env.sh && \
  cd /kb/dev_container/modules && \
  rm -rf jars && \
  git clone https://github.com/kbase/jars && \
  rm -rf kb_sdk && \
  git clone https://github.com/kbase/kb_sdk -b develop && \
  cd /kb/dev_container/modules/jars && \
  make deploy && \
  cd /kb/dev_container/modules/kb_sdk && \
  make && make deploy
RUN \
  . /kb/dev_container/user-env.sh && \
  cd /kb/dev_container/modules && \
  rm -rf data_api && \
  git clone https://github.com/kbase/data_api -b develop && \
  pip install /kb/dev_container/modules/data_api
RUN \
  . /kb/dev_container/user-env.sh && \
  cd /kb/dev_container/modules && \
  rm -rf genome_util && \
  git clone https://github.com/kbase/genome_util && \
  cd /kb/dev_container/modules/genome_util && \
  make && make deploy

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
  sh /kb/deps/kb_cufflinks/install-cufflinks.sh 

COPY ./ /kb/module
RUN \
  cd /kb/module && \
  make && make deploy
#  rm -rf narrative_method_store \
ENV PATH=$PATH:/kb/dev_container/modules/kb_sdk/bin
WORKDIR /kb/module
RUN mkdir -p /kb/module/work
RUN pip freeze | grep requests
ENTRYPOINT [ "./scripts/entrypoint.sh" ]
CMD [ ]
