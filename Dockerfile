#The base image is FreeSurfer's synthstrip package
FROM freesurfer/synthstrip:1.4

#Install pyants
#RUN python3 -m pip install antspyx
#RUN python3 -m pip install -v antspy

RUN apt-get -y update
RUN apt-get -y install git
RUN apt-get -y install build-essential cmake=3.26.1
RUN mkdir -p /code
RUN cd /code && git clone https://github.com/ANTsX/ANTsPy
RUN cd /code/ANTsPy && python3 ./setup.py install
