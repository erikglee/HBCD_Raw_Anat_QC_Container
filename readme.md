This repo contains the Dockerfile for building a simple utility to help visualize T1w, T2w, and QALAS
MRI data on infants. The pipeline assumes BIDS formatting of an input dataset and follows the argument
structure of a BIDS application. With that in mind it can be run as follows:

>> hbcd_qc /bids_dir /output_dir participant

Processing can also be restricted to a specific subject or session by adding --participant_label or
--session_id flags.

The pipeline looks for T1w/T2w/QALAS scans and tries to make snapshot png images for each of them. Once
an image is found, it is skull stripped by SynthStrip, then it is aligned to infant MNI space using a
12 DOF transformation in DIPY. Then a png is generated that displays three axial, coronal, and sagittal
views of the original MRI image that has been linearly transformed to infant MNI space. The color intensities
in the plot default to having a lower/upper bound on intensity of 0.3 and 1.7 times the modal signal intensity
within the brain mask.

If QALAS data is available in the BIDS dataset, png snapshots will only be created if there is one QALAS volume
with naming following the pattern ".../anat/...inv-2_QALAS.ni...". If a file satisfiying this naming is present,
the image will be registered to the T1w infant MNI template under the image_templates folder, same as would
occur if the utility instead found a T1w image. T2w images are instead registered to the T2w infant MNI template.



Related links:

SynthStrip - https://surfer.nmr.mgh.harvard.edu/docs/synthstrip/ 

DIPY - https://dipy.org/ 

TemplateFlow - https://www.templateflow.org/browse/ 

Infant Atlases - https://nist.mni.mcgill.ca/infant-atlases-0-4-5-years/ 

nibabel - https://nipy.org/nibabel/ 
