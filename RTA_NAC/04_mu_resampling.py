#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright University College London Hospital, 2020
Author: Eric Einspaenner, Institute of Nuclear Medicine
For internal research only.
"""
import os
import numpy
from art import tprint
import sirf.STIR as Pet
import sirf.Reg as Reg
import sirf.Reg as Eng_ref
import sirf.Reg as Eng_flo

#%% Definition of different path

py_path = os.getcwd()

# avoid cluttering of files, delete working-folder and start from scratch
working_folder = py_path + '/working'

# change the current working directory to the given path
os.chdir(working_folder)


#%% Create folders for results

path_NAC = working_folder + '/recon/NAC/'
path_tm = working_folder + '/tm/'
path_mu = working_folder + '/mu/'

if not os.path.exists(path_mu):
    os.makedirs(path_mu, mode=0o770)
    print('Create Folder: {}'.format(path_mu))


#%% Files
    
attn_file = '/home/eric/Dokumente/PythonProjects/SIRF/UKL_data/mu_Map/stir_mu_map.hv'  # .nii possible, requires ITK
print('mu-Map: {}'.format(attn_file))

attn_image = Pet.ImageData(attn_file)

#%% Do registration of mu map to NAC_ref
# =============================================================================
# 
# tprint(attn_file)
# 
# ref_file = path_NAC + 'NAC_0.nii'
# ref = Eng_ref.ImageData(ref_file)
# 
# 
# # settings for attn resampler
# reg_attn2nac = Reg.NiftyAladinSym()
# reg_attn2nac.set_reference_image(ref)
# reg_attn2nac.set_floating_image(attn_image)
# reg_attn2nac.process()
# 
# tm_attn2nac = reg_attn2nac.get_transformation_matrix_forward()
# =============================================================================


#%% resample mu-Map into each NAC space

ref_file = path_NAC + 'NAC_0.nii'
ref = Eng_ref.ImageData(ref_file)
flo = Eng_flo.ImageData(attn_image)

tprint('Start Resampling')

for i in range(30):
    print('Begin resampling mu-Maps: {}, with {}'.format(i, path_tm + 'tm_nac_inv_' + str(i) + '.txt'))
    tm = numpy.loadtxt(path_tm + 'tm_nac_inv_' + str(i) + '.txt')
    tm_fwd = Reg.AffineTransformation(tm)
    #tm_attn2ithNAC = tm_attn2nac * tm_fwd

    resampler = Reg.NiftyResample()
    resampler.set_reference_image(ref)
    resampler.set_floating_image(flo)
    resampler.add_transformation(tm_fwd)
    resampler.set_padding_value(0)
    resampler.set_interpolation_type_to_linear()
    mu_map_in_ith_NAC_space = resampler.forward(flo)

    mu_map_in_ith_NAC_space.write(path_mu + 'stir_mu_map_in_recon_space_' + str(i))
    Reg.ImageData(mu_map_in_ith_NAC_space).write(path_mu + 'mu_' + str(i) + '.nii')

    print('Finish resampling mu-Maps: {}'.format(i))


tprint('Finish Resampling')