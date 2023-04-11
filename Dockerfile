#The base image is FreeSurfer's synthstrip package
FROM freesurfer/synthstrip:1.4

#Install relevant python packages
RUN python3 -m pip install nibabel==3.2.2
RUN python3 -m pip install dipy==1.6.0
RUN python3 -m pip install matplotlib==3.3.4

#Make code and data directory
RUN mkdir /hbcd_code && mkdir /image_templates

#Copy over images
ADD image_templates/tpl-MNIInfant_cohort-1_res-1_mask-applied_T1w.nii.gz /image_templates/tpl-MNIInfant_cohort-1_res-1_mask-applied_T1w.nii.gz
ADD image_templates/tpl-MNIInfant_cohort-1_res-1_mask-applied_T2w.nii.gz /image_templates/tpl-MNIInfant_cohort-1_res-1_mask-applied_T2w.nii.gz

#Copy code, assign permissions
ADD run.py /hbcd_code/run.py
RUN chmod 555 -R /hbcd_code
ENV PATH="${PATH}:/hbcd_code"
RUN pipeline_name=hbcd_qc && cp /hbcd_code/run.py /hbcd_code/$pipeline_name

#Define entrypoint
ENTRYPOINT ["hbcd_qc"]
