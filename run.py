#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import os, glob
import dipy
import nibabel as nib
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator
from dipy import align
import dipy
import argparse
from dipy.align.imaffine import (AffineMap,
                                 MutualInformationMetric,
                                 AffineRegistration)
from dipy.align.transforms import (TranslationTransform3D,
                                   RigidTransform3D)

def register_images(input_file_path, output_destination, reuse_existing_output_file = True):
    '''
    
    '''
    
    #Input file path = full path to either T1w or T2w image
    #Output Destination = the top directory for outputs (i.e. same for all subjects)
    
    input_path_split = input_file_path.split('/')
    if 'ses-' in input_path_split[-3]:
        partial_sub_path = '/'.join(input_path_split[-4:])
        beginning_path = '/'.join(input_path_split[:-4])
    else:
        partial_sub_path = '/'.join(input_path_split[-3:])
        beginning_path = '/'.join(input_path_split[:-3])

    anat_out_dir = os.path.join(output_destination, os.path.dirname(partial_sub_path))
    if os.path.exists(anat_out_dir) == False:
        os.makedirs(anat_out_dir)
        
    if 'T1w' in partial_sub_path:
        contrast = 'T1w'
        stripped_out_file = os.path.join(output_destination, partial_sub_path).replace('T1w.nii.gz', 'masked-brain_T1w.nii.gz')
        final_registered_out_file = os.path.join(output_destination, partial_sub_path).replace('T1w.nii.gz', 'reg-MNIInfant_T1w.nii.gz')
    elif 'T2w' in partial_sub_path:
        contrast = 'T2w'
        stripped_out_file = os.path.join(output_destination, partial_sub_path).replace('T2w.nii.gz', 'masked-brain_T2w.nii.gz')
        final_registered_out_file = os.path.join(output_destination, partial_sub_path).replace('T2w.nii.gz', 'reg-MNIInfant_T2w.nii.gz')
        
    if os.path.exists(final_registered_out_file):
        print('Using already existing registered out file with name: {}'.format(final_registered_out_file))
        return final_registered_out_file
    
    os.system('python3 /freesurfer/mri_synthstrip -i {input_path} -o {output_skull_stripped_path}'.format(
        input_path = input_file_path, output_skull_stripped_path = stripped_out_file))
    
    print('Attempting Native to MNI Infant Registration using DIPY: ')
    template_image_path = '/image_templates/tpl-MNIInfant_cohort-1_res-1_mask-applied_{}.nii.gz'.format(contrast)
    template_image = dipy.io.image.load_nifti(template_image_path)
    original_image = dipy.io.image.load_nifti(input_file_path)
    registered_img = align.affine_registration(stripped_out_file, template_image[0], static_affine=template_image[1])
    registered_out_path = stripped_out_file.replace('{}.nii.gz'.format(contrast), 'reg-MNIInfant_{}.nii.gz'.format(contrast))
    dipy.io.image.save_nifti(registered_out_path,
                         registered_img[0], template_image[1])
    
    affine_map = AffineMap(registered_img[1],
                       original_image[0].shape, original_image[1],
                       original_image[0].shape, original_image[1])
    resampled = affine_map.transform(original_image[0])
    dipy.io.image.save_nifti(final_registered_out_file,
                     resampled, original_image[1])

    
    
    
    return final_registered_out_file, stripped_out_file

def make_slices_image(image_nifti_path, slice_info_dict, output_img_name, close_plot = True,
                     upsample_factor = 2, mask_path = None):
    '''Takes a nifti and plots slices of the nifti according to slices_info_dict
    
    Parameters
    ----------
    image_nifti_path : str
        Path to nifti image to make plot with
    slice_info_dict : dict
        Dictionary that formats how the picture
        will be formatted. See example below.
    output_img_name : str
        The name/full path of the image
        to be created by this function
    close_plot : bool, default True
        Whether to close the plot after it
        is rendered
    mask_path : str or None, default None
        A masked version of the brain (with signal intensities).
        This will can be used to help with image contrast.
        
    Example slice_info_dict. The first entry in each key's
    list dictates which plane is being imaged. The second
    entry indicates where (in RAS) the center of the plane
    should be placed. And the third and fourth entries dictate
    the range of voxels to be included in the slice. For example,
    100 would mean that 200 units are included. Larger values
    will make larger field of views. 
    
    slice_info_dict = {'coronal_1' : [0, -25, 100, 100],
                   'coronal_2' : [0, 0, 100, 100],
                   'coronal_3' : [0, 25, 100, 100],
                   'sagittal_1' : [1, -50, 100, 100],
                   'sagittal_2' : [1, 0, 100, 100],
                   'sagittal_3' : [1, 30, 100, 100],
                   'axial_1' : [2, -50, 100, 100],
                   'axial_2' : [2, 0, 100, 100],
                   'axial_3' : [2, 50, 100, 100]}
    
    '''
    
    #Load the nifti image
    nifti_image = nib.load(image_nifti_path)
    
    #Load nifti mask (assume voxels > 0.5 are good)
    if type(None) == type(mask_path):
        vmin = None
        vmax = None
    else:
        mask_data = nib.load(mask_path).get_fdata()
        mask_vals = mask_data[mask_data > 0.5]
        #vmin = np.percentile(mask_vals, 1)
        #vmax = np.percentile(mask_vals, 95)
        hist_results = np.histogram(mask_vals, bins = 100)
        modal_value = hist_results[1][np.argmax(hist_results[0])]
        vmin = modal_value*.3
        vmax = modal_value*1.7
    
    #Grab data + affine
    full_data = nifti_image.get_fdata()
    full_affine = nifti_image.affine
    
    #Setup interpolator in scipy so we can
    #resample the image in RAS units instead
    #of voxel units
    i = np.arange(0,full_data.shape[0])
    j = np.arange(0,full_data.shape[1])
    k = np.arange(0,full_data.shape[2])
    interp = RegularGridInterpolator((i, j, k), full_data, method = 'linear', bounds_error = False)
    
    inv_affine = np.linalg.inv(full_affine) #To get to RAS
    imgs = [] #List to store all of the individual slice pixel intensities
    
    #Make each of the slice images
    for temp_img in slice_info_dict.keys():
        temp_setup = slice_info_dict[temp_img]
        temp_slice = []
        
        #Upsample by a factor of 2
        for i in range(-1*temp_setup[2]*upsample_factor,temp_setup[2]*upsample_factor):
            for j in range(-1*temp_setup[3]*upsample_factor, temp_setup[3]*upsample_factor):
                i_hat = i/upsample_factor
                j_hat = j/upsample_factor
                if temp_setup[0] == 0:
                    temp_slice.append(np.matmul(inv_affine, np.array([temp_setup[1],i_hat,j_hat,1])))
                elif temp_setup[0] == 1:
                    temp_slice.append(np.matmul(inv_affine, np.array([i_hat,temp_setup[1],j_hat, 1])))
                elif temp_setup[0] == 2:
                    temp_slice.append(np.matmul(inv_affine, np.array([i_hat,j_hat,temp_setup[1],1])))
                else:
                    raise ValueError('Error: the second entry must be 0,1,2 to indicate slicing axis')
        vals = interp(np.array(temp_slice)[:,0:3])
        imgs.append(np.rot90(vals.reshape((temp_setup[2]*2*upsample_factor, temp_setup[3]*2*upsample_factor))))

    dim1 = imgs[0].shape[0]
    dim2 = imgs[0].shape[1]
    full_img_panel = np.zeros((dim1*3, dim2*3))
    for i, temp_img in enumerate(imgs):
        y = np.mod(i, 3)
        x = np.floor(i/3)
        full_img_panel[int(x*dim1):int((1+x)*dim1),int(y*dim2):int((1+y)*dim2)] = temp_img
    full_img_panel[np.where(np.isnan(full_img_panel))] = 0


    fig = plt.figure(dpi = 250)
    plt.imshow(full_img_panel, cmap = 'gist_gray', interpolation='nearest', vmin=vmin, vmax=vmax)
    plt.xticks([])
    plt.yticks([])
    plt.axis('off')
    plt.savefig(output_img_name, bbox_inches='tight')
    if close_plot:
        plt.close()
    
    return


#Configure the commands that can be fed to the command line
parser = argparse.ArgumentParser()
parser.add_argument("bids_dir", help="The path to the BIDS directory for your study (this is the same for all subjects)", type=str)
parser.add_argument("output_dir", help="The path to the folder where outputs will be stored (this is the same for all subjects)", type=str)
parser.add_argument("analysis_level", help="Should always be participant", type=str)

parser.add_argument('--participant_label', '--participant-label', help="The name/label of the subject to be processed (i.e. sub-01 or 01)", type=str)
parser.add_argument('--session_id', '--session-id', help="OPTIONAL: the name of a specific session to be processed (i.e. ses-01)", type=str)
args = parser.parse_args()


#Get cwd in case relative paths are given
cwd = os.getcwd()

#reassign variables to command line input
bids_dir = args.bids_dir
if os.path.isabs(bids_dir) == False:
    bids_dir = os.path.join(cwd, bids_dir)
output_dir = args.output_dir
if os.path.isabs(output_dir) == False:
    output_dir = os.path.join(cwd, output_dir)
analysis_level = args.analysis_level
if analysis_level != 'participant':
    raise ValueError('Error: analysis level must be participant, but program received: ' + analysis_level)


#Set session label
if args.session_id:
    session_label = args.session_id
    if 'ses-' not in session_label:
        session_label = 'ses-' + session_label
else:
    session_label = None
    
#Find participants to try running
if args.participant_label:
    participant_split = args.participant_label.split(' ')
    participants = []
    for temp_participant in participant_split:
        if 'sub-' not in temp_participant:
            participants.append('sub-' + temp_participant)
        else:
            participants.append(temp_participant)
else:
    os.chdir(bids_dir)
    participants = glob.glob('sub-*')
    
#Dictionary for making slices
slice_info_dict = {'coronal_1' : [0, -25, 100, 100],
                   'coronal_2' : [0, 0, 100, 100],
                   'coronal_3' : [0, 25, 100, 100],
                   'sagittal_1' : [1, -50, 100, 100],
                   'sagittal_2' : [1, 0, 100, 100],
                   'sagittal_3' : [1, 30, 100, 100],
                   'axial_1' : [2, -50, 100, 100],
                   'axial_2' : [2, 0, 100, 100],
                   'axial_3' : [2, 50, 100, 100]}
    
#Iterate through all participants
for temp_participant in participants:
    
    #Check that participant exists at expected path
    subject_path = os.path.join(bids_dir, temp_participant)
    if os.path.exists(subject_path):
        os.chdir(subject_path)
    else:
        raise AttributeError('Error: no directory found at: ' + subject_path)
    
    #Find session/sessions
    if session_label == None:
        sessions = glob.glob('ses*')
        if len(sessions) < 1:
            sessions = ['']
    elif os.path.exists(session_label):
        sessions = [session_label]
    else:
        raise AttributeError('Error: session with name ' + session_label + ' does not exist at ' + subject_path)

    #Iterate through sessions
    for temp_session in sessions:

        #If there is no session structure, this will go to the subject path
        session_path = os.path.join(subject_path, temp_session)

        #Grab T1w file
        anats_dict = {}
        t1_anats = glob.glob(os.path.join(session_path,'anat/*T1w.ni*'))
        anats_dict['T1w_images'] = t1_anats
        
        #Grab T2w file
        t2_anats = glob.glob(os.path.join(session_path,'anat/*T2w.ni*'))
        anats_dict['T2w_images'] = t2_anats
        
        for temp_t1w in anats_dict['T1w_images']:
            registered_nii_for_slice_img, masked_image = register_images(temp_t1w, output_dir)
            slice_img_path = registered_nii_for_slice_img.replace('T1w.nii', 'T1w_image-slice.png')
            slice_img_path = slice_img_path.replace('slice.png.gz', 'slice.png') #For case when nifti is compressed
            make_slices_image(registered_nii_for_slice_img, slice_info_dict, slice_img_path, close_plot = True,
                     upsample_factor = 2, mask_path = masked_image)
            
            
        for temp_t2w in anats_dict['T2w_images']:
            registered_nii_for_slice_img, masked_image = register_images(temp_t2w, output_dir)
            slice_img_path = registered_nii_for_slice_img.replace('T2w.nii', 'T2w_image-slice.png')
            slice_img_path = slice_img_path.replace('slice.png.gz', 'slice.png') #For case when nifti is compressed
            make_slices_image(registered_nii_for_slice_img, slice_info_dict, slice_img_path, close_plot = True,
                     upsample_factor = 2, mask_path = masked_image)

        print('Finished with: {}'.format(session_path))