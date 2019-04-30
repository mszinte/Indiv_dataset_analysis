## Indiv_dataset_analyis

- code for analysing Individual subject dataset
- subject data can be found here: https://osf.io/5bmdn/
- PRF output can be found here: https://osf.io/5bmdn/

## Analysis specifics
---------------------
- individual subjects are first preprocessed using FMRIPREP on LISA, see pre_fit/fmriprep_lisa.py
- Well registrered runs were kept for analysis, see select_block.json
- data are sg filter, psc and averaged within subjects, see pre_fit/pre_fit.py
- pRF parameters are extracted using fit/submit_fit_jobs.py on Lisa for gaussian model
- pRF derivatives of each subjects are analysed and drawn on fsaverage pycortex map of nprf_hcp using pp_roi.py
- each ROI masks are and saved in hdf5 files using post_pp_roi.py
- make Figure S1 with post_fit/notebooks/MakeFigureS1.ipynb
- make Figure S2 with post_fit/notebooks/MakeFigureS2.ipynb
