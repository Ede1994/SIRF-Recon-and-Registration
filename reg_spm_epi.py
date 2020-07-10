# Copyright University College London Hospital, 2020
# Author: Eric Einspaenner, Institute of Nuclear Medicine
# For internal research only.

#%% import all necessary modules
import os
import shutil
import sys
import glob
import re
import numpy
import math
import matplotlib.pyplot as plt
from art import *
import sirf.STIR as Pet
import sirf.Reg as Reg
import sirf.Reg as Eng_ref
import sirf.Reg as Eng_flo


#%% functions
# function for sorting lists, containing strings with numbers (https://stackoverflow.com/a/48413757)
def sorted_alphanumeric(data):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(data, key=alphanum_key, reverse=False)

# generate a list, contains the time intervalls (sum of frames)
def sum_list(lists):
    sum_liste = []
    add = int(0)
    sum_liste.append(add)
    for frame in lists[0:]:
        add += frame
        sum_liste.append(add)
    return sum_liste

# Set the STIR verbosity
Pet.set_verbosity(1)


#%% data path and set filenames
# data path EPI files
data_path_EPI = '/home/rich/Documents/ReCo/UKL_data/Processed/MedPhys_MoCo_Test_20190703_33/10_Head_ep2d_bold_moco/'
# avoid cluttering of files, delete working-folder and start from scratch
working_folder = '/home/rich/Documents/ReCo/working_EPI'
# if os.path.exists(working_folder):
#     shutil.rmtree(working_folder)
if not os.path.exists(working_folder):
    os.makedirs(working_folder, mode=0o770)

# change the current working directory to the given path
os.chdir(working_folder)


#%% redirect STIR messages to some files
# you can check these if things go wrong
msg_red = Pet.MessageRedirector('info.txt', 'warn.txt')


#%% create folders for results
path_epi = working_folder + '/EPI/'
path_moco = working_folder + '/moco'
if not os.path.exists(path_epi):
    os.makedirs(path_epi, mode=0o770)
    print('Create Folder: {}'.format(path_epi))
if not os.path.exists(path_moco):
    os.makedirs(path_moco, mode=0o770)
    print('Create Folder: {}'.format(path_moco))


#%% read EPI images as dcm and convert it to nii
# for i, file in zip(range(len(os.listdir(data_path_EPI))), os.listdir(data_path_EPI)):
#     print(data_path_EPI + file)
#     epi_image = Pet.ImageData(data_path_EPI + file)
#     temp = Reg.ImageData(epi_image)
#     temp.write(path_epi + 'EPI_' + str(i))


#%% SPM registration NAC
# define reference image (first image) and float-path, NAC
ref_file = path_epi + '20191016_131125Headep2dboldmocos010a001_00001.nii'
flo_path = path_epi

tprint('Start Reg for EPI')

# SPM with EPI images
spm_reg = Reg.SPMRegistration()
spm_reg.set_reference_image_filename(ref_file)
for image in sorted_alphanumeric(os.listdir(flo_path)):
    flo_file = flo_path + image
    print(flo_file)
    spm_reg.add_floating_image_filename(flo_file)
spm_reg.set_working_folder(path_moco)
spm_reg.set_working_folder_file_overwrite(True)
spm_reg.set_delete_temp_files(True)
spm_reg.process()

tprint('Finish Reg for EPI')


#%% get TM from Reg and calculate mean translation amplitude for every EPI volume
translation_values = []
rotation_values = []
# for i in range(len(os.listdir(flo_path))):
#     tm_fwd = spm_reg.get_transformation_matrix_forward(i)
#     tm_fwd.write(path_moco + '/tm' + str(i))
#     tm_fwD_arr = tm_fwd.as_array()
#     translation = numpy.sqrt((tm_fwD_arr[0][3])**2 + (tm_fwD_arr[1][3])**2 + (tm_fwD_arr[2][3])**2)
#     translation_values.append(translation)

for file in sorted_alphanumeric(os.listdir(path_moco)):
    print(path_moco + '/' + file)
    tm_fwD_arr = numpy.loadtxt(path_moco + '/' + file)
    translation = numpy.sqrt((tm_fwD_arr[0][3]) ** 2 + (tm_fwD_arr[1][3]) ** 2 + (tm_fwD_arr[2][3]) ** 2)
    translation_values.append(translation)
    rotation = numpy.sqrt((math.atan2(tm_fwD_arr[2][1], tm_fwD_arr[2][2])) ** 2 + (math.atan2(-tm_fwD_arr[2][0], numpy.sqrt(tm_fwD_arr[2][1] ** 2 + tm_fwD_arr[2][2] ** 2))) ** 2 + (math.atan2(tm_fwD_arr[1][0], tm_fwD_arr[0][0])) ** 2)
    rotation_values.append(rotation)

#%% save list of translation amplitude values as txt-file
with open(path_moco + '/translation.txt', 'w') as f:
    for item in translation_values:
        f.write("%s\n" % item)

with open(path_moco + '/rotation.txt', 'w') as f:
    for item in rotation_values:
        f.write("%s\n" % item)

tprint('DONE!')
