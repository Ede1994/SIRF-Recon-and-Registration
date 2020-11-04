# SIRF-Recon and Registration
 - RTA_EPI: Reconstruct-Transform-Average EPI-based
 - RTA_NAC: Reconstruct-Transform-Average NAC-based
 - Reg_Nifty: Registration with NiftyReg -> outputs TM and absolute displacement values
 
 - crerate_epi_nii.sh: convert .dcm files to .nii (needs ODIN)
 
 - motion_simu.py includes the possibility to perform a motion
 - reg_spm_NAC.py uses SPM for registration of NACs
 - reg_spm_NAC_EPI.py registration of NACs in EPI space
 - reg_spm_NAC_UTE.py registration of NACs in UTE space
 - reg_spm_epi.py uses SPM for registration of EPIs
 - sirf_recon_resample.py includes NiftyReg: resampling based on external transformation matrices
 - sirf_recon_spm_rta.py includes SPM for registration, RTA: Reconstruct-Transform-Add
