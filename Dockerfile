#The base image is FreeSurfer's synthstrip package
FROM freesurfer/synthstrip:1.4

#Install pyants
RUN python3 -m pip install antspyx