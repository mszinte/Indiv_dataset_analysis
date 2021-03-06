"""
-----------------------------------------------------------------------------------------
post_pp_roi.py
-----------------------------------------------------------------------------------------
Goal of the script:
Determine roi h5 files
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject number
sys.argv[2]: fit model ('gauss','css')
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
python post_fit/post_pp_roi.py sub-01 gauss
python post_fit/post_pp_roi.py sub-02 gauss
-----------------------------------------------------------------------------------------
"""

# Stop warnings
# -------------
import warnings
warnings.filterwarnings("ignore")

# General imports
# ---------------
import os
import sys
import json
import glob
import numpy as np
import matplotlib.pyplot as pl
import ipdb
import platform
import h5py
import scipy.io
opj = os.path.join
deb = ipdb.set_trace

# MRI imports
# -----------
import nibabel as nb
import cortex
from nipype.interfaces.freesurfer import SurfaceTransform

# Functions import
# ----------------
from plot_class import PlotOperator
from utils import set_pycortex_config_file, mask_gii_2_hdf5

# Bokeh imports
# ---------------
from bokeh.io import output_notebook, show, save, output_file, export_png, export_svgs
from bokeh.layouts import row, column, gridplot

# Popeye imports
import popeye.utilities as utils
from popeye.visual_stimulus import VisualStimulus
import popeye.css as css
import popeye.og as og

# Get inputs
# ----------
subject = sys.argv[1]
fit_model = sys.argv[2]

# Define analysis parameters
# --------------------------
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define cluster/server specific parameters
# -----------------------------------------
if 'aeneas' in platform.uname()[1]:
    base_dir = analysis_info['aeneas_base_folder'] 
elif 'lisa' in platform.uname()[1]:
    base_dir = analysis_info['lisa_base_folder'] 

deriv_dir = opj(base_dir,'pp_data',subject,fit_model,'deriv')
roi_masks_dir = opj(base_dir,'pp_data',subject,fit_model,'roi_masks')
h5_dir = opj(base_dir,'pp_data',subject,fit_model,'h5')
try: os.makedirs(roi_masks_dir)
except OSError: pass

# Check system
# ------------
sys.exit('Popeye error with Python 2. Use Python 3 Aborting.') if sys.version_info[0] < 3 else None

# Change cortex database folder
# -----------------------------
pycortex_folder     =   opj(base_dir,'pp_data','cortex')
set_pycortex_config_file(   project_folder = pycortex_folder)

# Create mask from overlay.svg
# ----------------------------
# make fsaverage roi masks
print('creating roi masks from overlay.svg')
masks = cortex.utils.get_roi_verts( subject = 'fsaverage', 
                                    roi = analysis_info['rois'], 
                                    mask = True)
mat_masks = []
for roi in analysis_info['rois']:
    mat_masks.append(masks[roi])
mat_masks = np.vstack(mat_masks)
mat_masks = mat_masks.astype('float32')

prf_deriv_L_all_fsaverage = nb.load(opj(deriv_dir,'all','prf_deriv_L_all_fsaverage.gii'))
mat_masks_L = mat_masks[:,0:163842]
darrays = [nb.gifti.gifti.GiftiDataArray(d) for d in mat_masks_L]
gii_out = nb.gifti.gifti.GiftiImage(header = prf_deriv_L_all_fsaverage.header,
                                    extra = prf_deriv_L_all_fsaverage.extra,
                                    darrays = darrays)
nb.save(gii_out,opj(roi_masks_dir,"masks_L_fsaverage.gii"))

prf_deriv_R_all_fsaverage = nb.load(opj(deriv_dir,'all','prf_deriv_R_all_fsaverage.gii'))
mat_masks_R = mat_masks[:,163842:327684]
darrays = [nb.gifti.gifti.GiftiDataArray(d) for d in mat_masks_R]
gii_out = nb.gifti.gifti.GiftiImage(header = prf_deriv_R_all_fsaverage.header, 
                                    extra = prf_deriv_R_all_fsaverage.extra, 
                                    darrays = darrays)
nb.save(gii_out,opj(roi_masks_dir,"masks_R_fsaverage.gii"))

# convert roi masks to fsaverage6
print('converting roi masks files to fsaverage6')
sxfm = SurfaceTransform()
sxfm.inputs.source_subject = "fsaverage"
sxfm.inputs.target_subject = "fsaverage6"
sxfm.terminal_output = 'none'

for hemi in ['L','R']:
    sxfm.inputs.subjects_dir = opj(base_dir,'derivatives','freesurfer')
    if hemi == 'L': sxfm.inputs.hemi = "lh"
    elif hemi == 'R': sxfm.inputs.hemi = "rh"
        
    sxfm.inputs.source_file = opj(roi_masks_dir,"masks_{hemi}_fsaverage.gii".format(hemi = hemi))
    sxfm.inputs.out_file = opj(roi_masks_dir,"masks_{hemi}_fsaverage6.gii".format(hemi = hemi))
    print(sxfm.inputs.out_file)
    sxfm.run()

# Save ROIS data in hdf5
# ----------------------
print('creating h5 files')
for roi_num, roi in enumerate(analysis_info['rois']):
    try: os.makedirs(h5_dir)
    except OSError: pass

    h5_file = opj(h5_dir,'{roi}.h5'.format(roi = roi))
    try: os.system('rm '+ h5_file)
    except: pass

    for hemi in ['L','R']:

        mask_file = opj(roi_masks_dir,"masks_{hemi}_fsaverage6.gii".format(hemi = hemi))
        
        for mask_dir in ['all','pos','neg']:
            
            in_file = opj(deriv_dir,mask_dir,"prf_deriv_{hemi}_{mask_dir}.gii".format(hemi = hemi, mask_dir = mask_dir))
            folder_alias = '{hemi}_{mask_dir}'.format(hemi = hemi,mask_dir = mask_dir)
            
            mask_gii_2_hdf5(in_file = in_file,
                            mask_file = mask_file,
                            hdf5_file = h5_file,
                            folder_alias = folder_alias,
                            roi_num = roi_num)

