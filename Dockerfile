FROM kbase/depl:latest
MAINTAINER KBase Developer
# Install the SDK (should go away eventually)
RUN \
apt-get update && apt-get install -y ant && \
cd /kb/dev_container/modules && \
git clone https://github.com/kbaseIncubator/KBaseRNASeq.git -b develop && \
cd /kb/dev_container/modules/KBaseRNASeq && \
make
# -----------------------------------------
# Insert apt-get instructions here to install
# any required dependencies for your module.
# -----------------------------------------
COPY ./ /kb/module
ENV PATH=$PATH:/kb/dev_container/modules/KBaseRNASeq/bin
WORKDIR /kb/module
RUN make
RUN make deploy
RUN mkdir -p /kb/module/work
ENTRYPOINT [ "./scripts/entrypoint.sh" ]
CMD [ ]
