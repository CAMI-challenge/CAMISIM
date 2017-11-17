FROM ubuntu:16.04

RUN apt update
RUN apt install -y python2.7 python-pip perl
RUN perl -MCPAN -e 'install XML::Simple'
ADD requirements.txt /requirements.txt
RUN cat requirements.txt | xargs -n 1 pip install
ADD *.py /usr/local/bin/
ADD scripts /usr/local/bin/scripts
ADD tools /usr/local/bin/tools
ADD defaults /usr/local/bin/defaults
WORKDIR /usr/local/bin
ENTRYPOINT ["python"]
