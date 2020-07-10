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
import matplotlib.pyplot as plt
from art import *
from sirf import SIRF
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
# data path to list-mode and normalisation files
data_path_LM = '/home/rich/Documents/ReCo/UKL_data/LM/20191009/'
data_path_EPI_image = '/home/rich/Documents/ReCo/working_NAC/EPI_01.nii'

# avoid cluttering of files, delete working-folder and start from scratch
working_folder = '/home/rich/Documents/ReCo/working_NAC'
# if os.path.exists(working_folder):
#     shutil.rmtree(working_folder)
if not os.path.exists(working_folder):
    os.makedirs(working_folder, mode=0o770)

# change the current working directory to the given path
os.chdir(working_folder)

# input files
list_file = data_path_LM + '1.3.12.2.1107.5.2.38.51008.30000019100905461135500000006.l.hdr'
print('LM data: {}'.format(list_file))
norm_file = data_path_LM + '1.3.12.2.1107.5.2.38.51008.30000019100905461135500000008.n.hdr'
print('Norm data: {}'.format(norm_file))


# output filename prefixes
sino_file = 'sino'


#%% define time frames (read frames.txt) and time intervals
# path to folder with frames.txt
frames_path = '/home/rich/tmp/processed/home/rich/Documents/ReCo/UKL_data/Processed/MedPhys_MoCo_Test_20191009_19/40008__new2_4moco_PRR_AC_Images/'

with open(frames_path + 'frames.txt', 'r') as f:
    time_frames = [int(line.rstrip()) for line in f]

print('Frames: {}'.format(time_frames))

# define time intervals
time_intervals = sum_list(time_frames)
print('Time intervals: {}'.format(time_intervals))


#%% redirect STIR messages to some files
# you can check these if things go wrong
msg_red = Pet.MessageRedirector('info.txt', 'warn.txt')


#%% create template and set lm2sino converter
# template for acq_data
template_acq_data = Pet.AcquisitionData('Siemens_mMR', span=11, max_ring_diff=16, view_mash_factor=2)
template_acq_data.write('template.hs')


#%% create listmode-to-sinograms converter object
lm2sino = Pet.ListmodeToSinograms()

# set input, output and template files
lm2sino.set_input(list_file)
lm2sino.set_output_prefix(sino_file)
lm2sino.set_template('template.hs')


#%% create folders for results
path_sino = working_folder + '/sino'
path_rando = working_folder + '/rando'
path_recon = working_folder + '/recon'
path_NAC = working_folder + '/recon/NAC'
path_moco_NAC = working_folder + '/moco/NAC'
if not os.path.exists(path_sino):
    os.makedirs(path_sino, mode=0o770)
    print('Create Folder: {}'.format(path_sino))
if not os.path.exists(path_rando):
    os.makedirs(path_rando, mode=0o770)
    print('Create Folder: {}'.format(path_rando))
if not os.path.exists(path_recon):
    os.makedirs(path_recon, mode=0o770)
    print('Create Folder: {}'.format(path_recon))
if not os.path.exists(path_NAC):
    os.makedirs(path_NAC, mode=0o770)
    print('Create Folder: {}'.format(path_NAC))
if not os.path.exists(path_moco_NAC):
    os.makedirs(path_moco_NAC, mode=0o770)
    print('Create Folder: {}'.format(path_moco_NAC))


#%% register EPI to PET with NiftyReg
ref_file = path_NAC + '/nii/' + 'NAC_0.nii'
flo_file = data_path_EPI_image
ref = Eng_ref.ImageData(ref_file)
flo = Eng_ref.ImageData(flo_file)

register = Reg.NiftyAladinSym()
register.set_reference_image(ref)
register.set_floating_image(flo)
register.process()
tm_fwd = register.get_transformation_matrix_forward()
tm_fwd.write(working_folder + '/tm_fwd.txt')
tm_inv = register.get_transformation_matrix_inverse()
tm_inv.write(working_folder + '/tm_inv.txt')

# idy = tm_inv * tm_fwd
# idy.write(working_folder + '/idy.txt')

matrix = numpy.loadtxt('tm_test_epi.txt')
inverse = numpy.linalg.inv(matrix)
tm_fwd2 = Reg.AffineTransformation(matrix)
tm_inv2 = Reg.AffineTransformation(inverse)
numpy.savetxt('tm_test_inv_epi.txt', inverse)

resampler = Reg.NiftyResample()
resampler.set_reference_image(ref)
resampler.set_floating_image(ref)
resampler.add_transformation(tm_inv)
resampler.set_padding_value(0)
resampler.set_interpolation_type_to_linear()
resampler.process()
temp = resampler.get_output()
temp.write(working_folder + '/new_pet')

ref = Eng_ref.ImageData(working_folder + '/new_pet.nii')

resampler = Reg.NiftyResample()
resampler.set_reference_image(ref)
resampler.set_floating_image(ref)
resampler.add_transformation(tm_fwd)
resampler.set_padding_value(0)
resampler.set_interpolation_type_to_linear()
#resampler.process()
temp = resampler.forward(ref)
temp.write(working_folder + '/new_pet2')

'''
ref = Eng_ref.ImageData(working_folder + '/new_pet2.nii')
matrix = numpy.loadtxt('tm_test.txt')
inverse = numpy.linalg.inv(matrix)
tm_fwd2 = Reg.AffineTransformation(matrix)
tm_inv2 = Reg.AffineTransformation(inverse)
numpy.savetxt('tm_test_inv.txt', inverse)

resampler = Reg.NiftyResample()
resampler.set_reference_image(ref)
resampler.set_floating_image(ref)
resampler.add_transformation(tm_fwd2)
resampler.set_padding_value(0)
resampler.set_interpolation_type_to_linear()
resampler.process()
temp = resampler.get_output()
temp.write(working_folder + '/new_pet3')
'''

i = 0
for image in sorted_alphanumeric(os.listdir('/home/rich/Documents/ReCo/working_NAC/recon/NAC/nii')):
    print(image)
    ref = Eng_ref.ImageData(working_folder + '/new_pet.nii')
    resampler.set_reference_image(ref)
    flo = Eng_ref.ImageData('/home/rich/Documents/ReCo/working_NAC/recon/NAC/nii/' + image)
    resampler.set_floating_image(flo)
    resampler.set_padding_value(0)
    resampler.set_interpolation_type_to_linear()
    resampler.process()
    new_NAC = resampler.get_output()
    new_NAC.write(path_recon + '/new_NAC/new_NAC' + str(i))
    i += 1


#%% SPM registration NAC
#define reference image (first image) and float-path
ref_file = path_recon + '/new_NAC/new_NAC0.nii'
flo_path = path_recon + '/new_NAC/'

tprint('Start Reg for NAC')

# SPM with NAC images
spm_reg = Reg.SPMRegistration()
spm_reg.set_reference_image_filename(ref_file)
for image in sorted_alphanumeric(os.listdir(flo_path)):
    flo_file = flo_path + image
    print(flo_file)
    spm_reg.add_floating_image_filename(flo_file)
spm_reg.set_working_folder(path_moco_NAC)
spm_reg.set_working_folder_file_overwrite(True)
spm_reg.set_delete_temp_files(True)
spm_reg.process()

tprint('Finish Reg for NAC')


#%% get TM from Reg and calculate mean translation amplitude for every EPI volume
translation_values = []
for i in range(len(os.listdir(flo_path))):
    tm_fwd = spm_reg.get_transformation_matrix_forward(i)
    tm_fwd.write(path_moco_NAC + '/tm' + str(i))
    # tm_fwD_arr = tm_fwd.as_array()
    # translation = numpy.sqrt((tm_fwD_arr[0][3])**2 + (tm_fwD_arr[1][3])**2 + (tm_fwD_arr[2][3])**2)
    # translation_values.append(translation)

for file in sorted_alphanumeric(os.listdir(path_moco_NAC)):
    print(path_moco_NAC + '/' + file)
    tm_fwD_arr = numpy.loadtxt(path_moco_NAC + '/' + file)
    translation = numpy.sqrt((tm_fwD_arr[0][3]) ** 2 + (tm_fwD_arr[1][3]) ** 2 + (tm_fwD_arr[2][3]) ** 2)
    translation_values.append(translation)


#%% save list of translation amplitude values as txt-file
with open(path_moco_NAC + '/translation.txt', 'w') as f:
    for item in translation_values:
        f.write("%s\n" % item)

tprint('DONE!')