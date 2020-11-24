#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright University College London Hospital, 2020
Author: Eric Einspaenner, Institute of Nuclear Medicine
For internal research only.
"""
import os
import re
from art import tprint
import sirf.Reg as Reg
import sirf.Reg as Eng_ref
import sirf.Reg as Eng_flo


#%% Functions

# function for sorting lists, containing strings with numbers (https://stackoverflow.com/a/48413757)
def sorted_alphanumeric(data):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(data, key=alphanum_key, reverse=False)


# registration of EPI images
def reg_nac(ref, flo_path, list_NACs):
    i = 0
    nac_reg = Reg.NiftyAladinSym()
    nac_reg.set_reference_image(ref)
    for image in sorted_alphanumeric(list_NACs):
        flo_file = flo_path + image
        print(flo_file)
        flo = Eng_flo.ImageData(flo_file)
        nac_reg.set_floating_image(flo)
        nac_reg.process()
        tm_nac = nac_reg.get_transformation_matrix_forward()
        tm_nac_inv = nac_reg.get_transformation_matrix_inverse()
        tm_nac.write(path_tm + 'tm_nac_' + str(i))
        tm_nac_inv.write(path_tm + 'tm_nac_inv_' + str(i))
        i += 1


#%% Definition of different path

py_path = os.getcwd()

# avoid cluttering of files, delete working-folder and start from scratch
working_folder = py_path + '/working'

# change the current working directory to the given path
os.chdir(working_folder)


#%% Create folders for results

path_NAC = working_folder + '/recon/NAC/'
path_smooth = working_folder + '/recon/SMOOTH/'
path_tm = working_folder + '/tm/'

if not os.path.exists(path_tm):
    os.makedirs(path_tm, mode=0o770)
    print('Create Folder: {}'.format(path_tm))


#%% Registration NAC, delivers transformation matrices 
# define reference image (first image) and float-path

tprint('Start Registration of NACs')

# refernce file
ref_file = path_smooth + 'smooth_0.nii'
ref = Eng_ref.ImageData(ref_file)

# float files
flo_path = path_smooth
list_smooth = [f for f in os.listdir(flo_path) if f.endswith(".nii")]

# Niftyreg with EPI images
reg_nac(ref, flo_path, list_smooth)

tprint('Finish Registration')