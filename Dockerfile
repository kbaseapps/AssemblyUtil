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

# Copy module files to image
COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
