FROM kbase/depl:latest
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
  make
RUN \
  . /kb/dev_container/user-env.sh && \
  cd /kb/dev_container/modules && \
  rm -rf genome_util && \
  git clone https://github.com/kbase/genome_util && \
  cd /kb/dev_container/modules/genome_util && \
  make && make deploy
CMD \
  cd /kb/module ; \
  rm -rf KBaseRNASeq ;\
  git clone https://github.com/kbase/KBaseRNASeq ;\
  cd /kb/module/KBaseRNASeq ;\
  exec deps/kb_tophat/install-tophat.sh
  exec deps/kb_bowtie/install-bowtie2.sh

####END OF KBASE #############################
#apt-get update && apt-get install -y ant && \
# -----------------------------------------
# Insert apt-get instructions here to install
# any required dependencies for your module.
# -----------------------------------------
COPY ./ /kb/module
ENV PATH=$PATH:/kb/module/bin:/kb/dev_container/modules/kb_sdk/bin
WORKDIR /kb/module
RUN make
#RUN make deploy
RUN mkdir -p /kb/module/work
ENTRYPOINT [ "./scripts/entrypoint.sh" ]
CMD [ ]
