FROM ubuntu:bionic
RUN apt-get update && \
apt-get install -y time vim git
RUN apt-get -y install python3-pip
RUN pip3 install cgatcore 
RUN pip3 install pyyaml
RUN pip3 install ruffus
RUN pip3 install gevent
RUN pip3 install paramiko
RUN pip3 install pyfiglet==0.7
RUN pip3 install pep8
RUN pip3 install pytest
RUN pip3 install pytest-pep8
RUN pip3 install drmaa
RUN pip3 install setuptools
RUN pip3 install six
RUN pip3 install sqlalchemy
RUN pip3 install apsw
RUN pip3 install pandas
RUN pip3 install numpy
RUN pip3 install pycoreutils
RUN pip3 install boto3
RUN pip3 install google-cloud-storage
RUN pip3 install google-cloud
RUN pip3 install ftputil
RUN pip3 install pysftp

RUN git clone https://github.com/flixofarrell/GQSLAP.git