#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright University College London Hospital, 2020
Author: Eric Einspaenner, Institute of Nuclear Medicine
For internal research only.
"""
import os
import shutil
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
def reg_epi(ref, flo_path):
    i = 0
    epi_reg = Reg.NiftyAladinSym()
    epi_reg.set_reference_image(ref)
    for image in sorted_alphanumeric(os.listdir(flo_path)):
        flo_file = flo_path + image
        print(flo_file)
        flo = Eng_flo.ImageData(flo_file)
        epi_reg.set_floating_image(flo)
        epi_reg.process()
        tm_epi = epi_reg.get_transformation_matrix_forward()
        tm_epi_inv = epi_reg.get_transformation_matrix_inverse()
        tm_epi.write(path_EPI + 'tm_epi_' + str(i))
        tm_epi_inv.write(path_EPI + 'tm_epi_inv_' + str(i))
        i += 1


#%% Definition of different path

py_path = os.getcwd()

# avoid cluttering of files, delete working-folder and start from scratch
working_folder = py_path + '/working'
if os.path.exists(working_folder):
    shutil.rmtree(working_folder)
if not os.path.exists(working_folder):
    os.makedirs(working_folder, mode=0o770)

# change the current working directory to the given path
os.chdir(working_folder)

path_EPI = working_folder + '/EPI/'
if not os.path.exists(path_EPI):
    os.makedirs(path_EPI, mode=0o770)
    print('Create Folder: {}'.format(path_EPI))


#%% Registration EPI, delivers transformation matrices 
# define reference image (first image) and float-path

tprint('Start Registration of EPIs')

# refernce file
epi_data_path = py_path + '/UKL_data/EPI/1/'
ref_file = epi_data_path + 'epi_0.nii'
ref = Eng_ref.ImageData(ref_file)

# float files
flo_path = epi_data_path

# Niftyreg with EPI images
reg_epi(ref, flo_path)

tprint('Finish Registration')