FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------

# pin biopython version
RUN pip install --upgrade biopython==1.70
RUN pip install --upgrade pytest==7.0.1 coverage==6.2 pytest-cov==4.0.0 python-dateutil==2.8.2 dill==0.3.4

# Copy module files to image
COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
