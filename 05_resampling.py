#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright University College London Hospital, 2020
Author: Eric Einspaenner, Institute of Nuclear Medicine
For internal research only.
"""
import os
import re
import numpy
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

# generate a list, contains the time intervalls (sum of frames)
def sum_list(lists):
    sum_liste = []
    add = int(0)
    sum_liste.append(add)
    for frame in lists[0:]:
        add += frame
        sum_liste.append(add)
    return sum_liste


#%% Definition of different path

py_path = os.getcwd()
print(py_path)

# data path to list-mode and normalisation files
data_path_LM = py_path + '/UKL_data/LM/20190723/20190723/'

# avoid cluttering of files, delete working-folder and start from scratch
working_folder = py_path + '/working'

# change the current working directory to the given path
os.chdir(working_folder)


#%% create folders for results
path_NAC = working_folder + '/recon/NAC/'
path_EPI = working_folder + '/EPI/'
path_moco = working_folder + '/moco/'
if not os.path.exists(path_NAC):
    os.makedirs(path_NAC, mode=0o770)
    print('Create Folder: {}'.format(path_NAC))
if not os.path.exists(path_EPI):
    os.makedirs(path_EPI, mode=0o770)
    print('Create Folder: {}'.format(path_EPI))
if not os.path.exists(path_moco):
    os.makedirs(path_moco, mode=0o770)
    print('Create Folder: {}'.format(path_moco))


#%% Define time frames (read frames.txt) and time intervals

# path to folder with frames.txt
frames_path = py_path + '/UKL_data/frames/'

# read frames.txt
with open(frames_path + 'frames.txt', 'r') as f:
    time_frames = [int(line.rstrip()) for line in f]
print('Frames: {}'.format(time_frames))

# number of motion steps
num_motion_steps = len(time_frames)

# define time intervals
time_intervals = sum_list(time_frames)
print('Time intervals: {}'.format(time_intervals))


# define the tm in the middle of the motion states
num_tm = []
for time, dur in zip(time_intervals, time_frames):
    t = int(time/2 + dur/4)
    num_tm.append(t)
print('Correct TM: {}'.format(num_tm))


#%% convert a array to a SIRF transformation matrix and then resample the float image

# list of all NACs
list_NACs = [f for f in os.listdir(path_NAC) if f.endswith(".nii")]

# define space matrices
tm_fwd = numpy.loadtxt(py_path + '/UKL_data/tm_epi/reg_NAC_EPI.txt')
tm_inv = numpy.loadtxt(py_path + '/UKL_data/tm_epi/reg_NAC_EPI_inv.txt')

# define reference image and float-path
ref_file = path_NAC + 'NAC_0.nii'
ref = Eng_ref.ImageData(ref_file)
flo_path = path_NAC

# settings for image resampler
resampler_im = Reg.NiftyResample()
resampler_im.set_reference_image(ref)
resampler_im.set_padding_value(0)
resampler_im.set_interpolation_type_to_linear()

tprint('Start Resampling')

i = 0

# for loop, simultaneous matrices and images
for num, image in zip(num_tm, sorted_alphanumeric(list_NACs)):
    print('TM: {}, Float-Image: {}'.format('tm_epi_' + str(num) + '.txt', image))

    flo = Eng_flo.ImageData(path_NAC + image)

    # read tm-matrix as numpy array and read load float image
    matrix = numpy.loadtxt(path_EPI + 'tm_epi_' + str(num) + '.txt')
    matrix2 = numpy.loadtxt(path_EPI + 'tm_epi_inv_' + str(num) + '.txt')
    
    # tm space transformation: EPI -> NAC
    # transform tm into PET space
    matrix2 = tm_inv * matrix2 * tm_fwd

    # create affine transformation from numpy array
    tm = Reg.AffineTransformation(matrix2)

    # motion information from another source and resample an image with this informations
    print('Begin resampling: {}'.format(image))
    resampler_im.clear_transformations()
    resampler_im.set_floating_image(flo)
    resampler_im.add_transformation(tm)
    new_im = resampler_im.forward(flo)
    new_im.write(working_folder + '/moco/moco_'+str(i))

    print('Resampling successful: {}'.format(image))

    i += 1

tprint('Finish Resampling')


#%% define RTA method

# define initial image (first image, first frame)
initial = Reg.NiftiImageData(working_folder + '/moco/moco_0.nii')
initial_array = initial.as_array()

# sum over all images (as array)
for image in sorted_alphanumeric(os.listdir(path_moco))[1:]:
    print(image)
    array = Reg.NiftiImageData(path_moco + image).as_array()
    initial_array += array

# create image
final_image = Reg.NiftiImageData(working_folder + '/moco/moco_0.nii')
final_image.fill(initial_array)
final_image.write('final_image_RTA.nii')

tprint('DONE!')

