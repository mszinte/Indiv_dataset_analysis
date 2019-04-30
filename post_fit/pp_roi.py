"""
-----------------------------------------------------------------------------------------
pp_roi.py
-----------------------------------------------------------------------------------------
Goal of the script:
Region of interests pre-processing
Compute pRF derivatives and plot on pycortex overlay.svg to determine visual ROI
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject number
sys.argv[2]: fit model ('gauss','css')
sys.argv[3]: voxels per fit (e.g 2500)
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
python post_fit/pp_roi.py sub-01 gauss 2500
python post_fit/pp_roi.py sub-02 gauss 2500
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
opj = os.path.join
deb = ipdb.set_trace

# MRI imports
# -----------
import nibabel as nb
import cortex
from nipype.interfaces.freesurfer import SurfaceTransform

# Functions import
# ----------------
from utils import set_pycortex_config_file, convert_fit_results

# Check system
# ------------
sys.exit('Popeye error with Python 2. Use Python 3 Aborting.') if sys.version_info[0] < 3 else None

# Get inputs
# ----------
subject = sys.argv[1]
fit_model = sys.argv[2]
job_vox = float(sys.argv[3])
if fit_model == 'gauss': fit_val = 6
elif fit_model == 'css': fit_val = 7
base_file_name = "{sub}_task-prf_space-fsaverage6.func_sg_psc".format(sub = subject)

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

# determine number of vertex and time_serie
data = []
data_file  =  sorted(glob.glob(opj(base_dir,'raw_data',subject,"{bfn}.gii".format(bfn = base_file_name))))
data_file_load = nb.load(data_file[0])
data.append(np.array([data_file_load.darrays[i].data for i in range(len(data_file_load.darrays))]))
data = np.vstack(data) 
ts_num,vox_num = data.shape[0],data.shape[1]

# Check if all slices are present
# -------------------------------
start_idx =  np.arange(0,vox_num,job_vox)
end_idx = start_idx+job_vox
end_idx[-1] = vox_num
num_miss_part = 0

fit_est_files_job = []
for iter_job in np.arange(0,start_idx.shape[0],1):
    fit_est_file = opj(base_dir,'pp_data',subject,fit_model,'fit', "{bfn}_est_{start_idx}_to_{end_idx}.gii".format(bfn = base_file_name,
                                                                                                                   start_idx = str(int(start_idx[iter_job])),
                                                                                                                   end_idx = str(int(end_idx[iter_job]))))
    if os.path.isfile(fit_est_file):
        if os.path.getsize(fit_est_file) == 0:
            num_miss_part += 1 
        else:
            fit_est_files_job.append(fit_est_file)
    else:
        num_miss_part += 1

if num_miss_part != 0:
    sys.exit('%i missing files, analysis stopped'%num_miss_part)

# Combine fit files
# -----------------
print('combining fit files')

data = np.zeros((fit_val,vox_num))

for fit_filename in fit_est_files_job:
    data_fit = []
    data_fit_file = nb.load(fit_filename)
    data_fit.append(np.array([data_fit_file.darrays[i].data for i in range(len(data_fit_file.darrays))]))
    data_fit = np.vstack(data_fit)
    data = data + data_fit

# Seperate hemi files
# -------------------
data_L = data[:,0:int(vox_num/2)]
data_R = data[:,int(vox_num/2):vox_num]

vox_num = int(vox_num/2.0)
for hemi in ['L','R']:
    exec("darrays_est_{hemi} = [nb.gifti.gifti.GiftiDataArray(d) for d in data_{hemi}]".format(hemi = hemi))
    exec("gii_out_{hemi} = nb.gifti.gifti.GiftiImage(header = data_fit_file.header, extra = data_fit_file.extra, darrays = darrays_est_{hemi})".format(hemi = hemi))
    exec("prf_filename_{hemi} = opj(base_dir,'pp_data',subject,fit_model,'fit','{bfn}_est_{hemi}.gii')".format(bfn =base_file_name, hemi = hemi))
    exec("nb.save(gii_out_{hemi}, prf_filename_{hemi})".format(hemi = hemi))

# Compute derived measures from prfs
# ----------------------------------
print('extracting pRF derivatives')
for hemi in ['L','R']:
    exec("prf_filename = prf_filename_{hemi}".format(hemi = hemi)) 
    convert_fit_results(prf_filename = prf_filename,
                        output_dir = deriv_dir,
                        hemi = hemi,
                        stim_width = analysis_info['stim_width'],
                        stim_height = analysis_info['stim_height'],
                        fit_model = fit_model)

# Resample gii to fsaverage
# -------------------------
print('converting derivative files to fsaverage')
sxfm = SurfaceTransform()
sxfm.inputs.source_subject = "fsaverage6"
sxfm.inputs.target_subject = "fsaverage"
sxfm.terminal_output = 'none'

for hemi in ['L','R']:
    sxfm.inputs.subjects_dir = opj(base_dir,'derivatives','freesurfer')
    if hemi == 'L': sxfm.inputs.hemi = "lh"
    elif hemi == 'R': sxfm.inputs.hemi = "rh"
        
    for mask_dir in ['all','pos','neg']:
        sxfm.inputs.source_file = opj(deriv_dir,mask_dir,"prf_deriv_{hemi}_{mask_dir}.gii".format(hemi = hemi, mask_dir = mask_dir))
        sxfm.inputs.out_file = opj(deriv_dir,mask_dir,"prf_deriv_{hemi}_{mask_dir}_fsaverage.gii".format(hemi = hemi, mask_dir = mask_dir))
        print(sxfm.inputs.out_file)
        sxfm.run()