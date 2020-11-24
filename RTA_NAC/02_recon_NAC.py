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
import sirf.STIR as Pet
import sirf.Reg as Reg


#%% Functions

# function for sorting lists, containing strings with numbers (https://stackoverflow.com/a/48413757)
def sorted_alphanumeric(data):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(data, key=alphanum_key, reverse=False)


#%% Definition of different path

py_path = os.getcwd()
print(py_path)

# data path to list-mode and normalisation files
data_path_LM = '/home/eric/Dokumente/PythonProjects/SIRF/UKL_data/LM/20190723/20190723/'

# avoid cluttering of files, delete working-folder and start from scratch
working_folder = py_path + '/working'

# change the current working directory to the given path
os.chdir(working_folder)
norm_file = data_path_LM + '1.3.12.2.1107.5.2.38.51008.30000019072305114947500000008.n.hdr'
print('Norm data: {}'.format(norm_file))


#%% Create folders for results

path_NAC = working_folder + '/recon/NAC/'
path_smooth = working_folder + '/recon/SMOOTH/'
path_sino = working_folder + '/sino/'
path_rando = working_folder + '/rando/'

if not os.path.exists(path_NAC):
    os.makedirs(path_NAC, mode=0o770)
    print('Create Folder: {}'.format(path_NAC))
if not os.path.exists(path_smooth):
    os.makedirs(path_smooth, mode=0o770)
    print('Create Folder: {}'.format(path_smooth))

#%% Settings for reconstruction

# set acq_model
acq_model = Pet.AcquisitionModelUsingRayTracingMatrix()
acq_model.set_num_tangential_LORs(5)

# set recon, OSEM
recon = Pet.OSMAPOSLReconstructor()
num_subsets = 21
num_subiterations = 168
recon.set_num_subsets(num_subsets)
recon.set_num_subiterations(num_subiterations)

# definitions for detector sensitivity modelling
asm_norm = Pet.AcquisitionSensitivityModel(norm_file)
acq_model.set_acquisition_sensitivity(asm_norm)

# dimensions
nxny = (256, 256)

#%% redirect STIR messages to some files
# you can check these if things go wrong
msg_red = Pet.MessageRedirector('info.txt', 'warn.txt')


#%% List of sino and rand
list_sino = [f for f in os.listdir(working_folder + '/sino/') if f.endswith(".hs")]
list_rando = [f for f in os.listdir(working_folder + '/rando/') if f.endswith(".hs")]


#%% NAC reconstruction

tprint('Start NAC Recon')

for i, sino, random in zip(range(len(path_sino)), sorted_alphanumeric(list_sino), sorted_alphanumeric(list_rando)):

    sino_pet = Pet.AcquisitionData(path_sino + sino)
    print(sino)
    randoms_pet = Pet.AcquisitionData(path_rando + random)
    print(random)

    # reconstruct the data (without mu-map)
    obj_fun = Pet.make_Poisson_loglikelihood(sino_pet)
    acq_model.set_background_term(randoms_pet)
    recon.set_objective_function(obj_fun)
    initial_image = sino_pet.create_uniform_image(1.0, nxny)
    image = initial_image
    recon.set_up(image)

    recon.set_current_estimate(image)
    recon.process()

    # save recon images
    recon_image = Reg.NiftiImageData(recon.get_output())
    recon_image.write(path_NAC + 'NAC_' + str(i))
    
    image = recon.get_output()
    
    # apply gaussian filter with 3mm fwhm
    gaussian_filter = Pet.SeparableGaussianImageFilter()
    gaussian_filter.set_fwhms((3, 3, 3))
    #gaussian_filter.set_max_kernel_sizes((10, 10, 2))
    gaussian_filter.set_normalise()
    gaussian_filter.set_up(image)
    gaussian_filter.apply(image)    

    # save Image as .nii
    smoothed_image = Reg.NiftiImageData(image)
    smoothed_image.write(path_smooth + 'smooth_' + str(i))

    print('Reconstruction successful: Frame {}'.format(i))

tprint('Finish NAC Recon')