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


#%% settings for reconstruction
# set acq_model
acq_model = Pet.AcquisitionModelUsingRayTracingMatrix()
acq_model.set_num_tangential_LORs(5)

# set recon, OSEM
recon = Pet.OSMAPOSLReconstructor()
num_subsets = 7
num_subiterations = 4
recon.set_num_subsets(num_subsets)
recon.set_num_subiterations(num_subiterations)

# definitions for detector sensitivity modelling
asm_norm = Pet.AcquisitionSensitivityModel(norm_file)
acq_model.set_acquisition_sensitivity(asm_norm)


#%% for-loop, reconstruction of time intervals
# NAC reconstruction
tprint('Start Recon NAC')
for i in range(len(time_intervals)-1):
    print('Begin reconstruction: Frame {}'.format(i))

    # listmode-to-sinogram
    lm2sino.set_time_interval(time_intervals[i], time_intervals[i+1])
    lm2sino.set_up()
    lm2sino.process()
    acq_data = lm2sino.get_output()
    acq_data.write(path_sino + '/sino' + str(i))

    # randoms estimate
    randoms = lm2sino.estimate_randoms()
    randoms.write(path_rando + '/rando' + str(i))

    # reconstruct the data (without mu-map)
    obj_fun = Pet.make_Poisson_loglikelihood(acq_data)
    recon.set_objective_function(obj_fun)
    initial_image = acq_data.create_uniform_image(1.0)
    image = initial_image
    recon.set_up(image)

    recon.set_current_estimate(image)
    recon.process()

    # save recon images
    recon_image = recon.get_output()
    recon_image.write(path_NAC + '/NAC_' + str(i))

    # save Image as .nii
    recon_image = Reg.NiftiImageData(recon.get_output())
    recon_image.write(path_NAC + '/nii/NAC_' + str(i))

    print('Reconstruction successful: Frame {}'.format(i))

tprint('Finish Recon NAC')


#%% SPM registration NAC
# define reference image (first image) and float-path, NAC
ref_file = path_NAC + '/nii/' + 'NAC_0.nii'
flo_path = path_NAC + '/nii/'
#ref = Eng_ref.ImageData(ref_file)
#flo = Eng_ref.ImageData(flo_file)

tprint('Start Reg for NAC')

# SPM with NAC images
spm_reg = Reg.SPMRegistration()
spm_reg.set_reference_image_filename(ref_file)
for image in sorted_alphanumeric(os.listdir(flo_path)):
    flo_file = flo_path + image
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
    tm_fwD_arr = tm_fwd.as_array()
    translation = numpy.sqrt((tm_fwD_arr[0][3])**2 + (tm_fwD_arr[1][3])**2 + (tm_fwD_arr[2][3])**2)
    translation_values.append(translation)


#%% save list of translation amplitude values as txt-file
with open(path_moco_NAC + '/translation.txt', 'w') as f:
    for item in translation_values:
        f.write("%s\n" % item)

tprint('DONE!')