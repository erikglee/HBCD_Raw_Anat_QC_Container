#The base image is FreeSurfer's synthstrip package
FROM freesurfer/synthstrip:1.4

#Install pyants
#RUN python3 -m pip install antspyx
#RUN python3 -m pip install -v antspy

RUN apt-get -y update
RUN apt-get -y install git
RUN apt-get install build-essential cmake --no-install-recommends
RUN mkdir -p /code
RUN cd /code && git clone https://github.com/ANTsX/ANTsPy
RUN cd /code/ANTsPy && python3 ./setup.py install
