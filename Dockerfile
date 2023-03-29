#The base image is FreeSurfer's synthstrip package
FROM freesurfer/synthstrip:1.4

#Install pyants
#RUN python3 -m pip install antspyx
#RUN python3 -m pip install -v antspy

RUN mkdir -p /code
RUN cd /code & git clone https://github.com/ANTsX/ANTsPy
RUN cd /code & python3 setup.py install
