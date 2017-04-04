FROM kbase/kbase:sdkbase.latest
MAINTAINER Michael Sneddon (mwsneddon@lbl.gov)
# -----------------------------------------

RUN pip install coverage

# update security libraries in the base image
RUN pip install cffi --upgrade \
    && pip install pyopenssl --upgrade \
    && pip install ndg-httpsclient --upgrade \
    && pip install pyasn1 --upgrade \
    && pip install requests --upgrade \
    && pip install 'requests[security]' --upgrade

# update / install biopython
RUN pip install --upgrade setuptools pip \
    && pip install --upgrade biopython==1.66


# Install KBase Data API Library + dependencies
RUN mkdir -p /kb/module && cd /kb/module && git clone https://github.com/kbase/data_api && \
    cd data_api && git checkout 0.4.0-dev && cd /kb/module && \
    mkdir -p lib/ && cp -a data_api/lib/doekbase lib/ && \
    pip install -r /kb/module/data_api/requirements.txt


# Install KBase Transform Scripts + dependencies
# Note: may not always be safe to copy things to /kb/deployment/lib
#RUN mkdir -p /kb/module && cd /kb/module && git clone https://github.com/kbase/transform && \
#    cd transform && git checkout 53633ed && cd /kb/module && \
#    mkdir -p lib/ && cp -a transform/lib/biokbase /kb/deployment/lib/ && \
#    pip install -r /kb/module/data_api/requirements.txt


# Copy module files to image
COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
