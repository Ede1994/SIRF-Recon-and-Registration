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
# function to display an image
def imshow(image, limits, title=''):
    """Display an image with a colourbar, returning the plot handle. 

    Arguments:
    image -- a 2D array of numbers
    limits -- colourscale limits as [min,max]. An empty [] uses the full range
    title -- a string for the title of the plot (default "")
    """
    plt.title(title)
    bitmap = plt.imshow(image)
    if len(limits) == 0:
        limits = [image.min(), image.max()]

    plt.clim(limits[0], limits[1])
    plt.colorbar(shrink=.6)
    plt.axis('on')
    return bitmap

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
data_path_LM = '/home/rich/Documents/ReCo/UKL_data/LM/20190417/'
data_path_MR = '/home/rich/Documents/ReCo/UKL_data/Processed/MedPhys_MoCo_Test_20190417_7/'
# avoid cluttering of files, delete working-folder and start from scratch
working_folder = '/home/rich/Documents/ReCo/working'
if os.path.exists(working_folder):
    shutil.rmtree(working_folder)
if not os.path.exists(working_folder):
    os.makedirs(working_folder, mode=0o770)

# change the current working directory to the given path
os.chdir(working_folder)

# input files
list_file = data_path_LM + '1.3.12.2.1107.5.2.38.51008.30000019041705474367000000006.l.hdr'
print('LM data: {}'.format(list_file))
norm_file = data_path_LM + '1.3.12.2.1107.5.2.38.51008.30000019041705474367000000008.n.hdr'
print('Norm data: {}'.format(norm_file))
attn_file = data_path_MR + 'stir_mu_map.hv'  # .nii possible, requires ITK
print('mu-Map: {}'.format(attn_file))

# output filename prefixes
sino_file = 'sino'


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


#%% from template sinogram, ensure that mu-map has spacing/offset
# that matches the reconstructed image
attn_image = Pet.ImageData(attn_file)
template_image = template_acq_data.create_uniform_image(1.0)

resampler = Reg.NiftyResample()
resampler.set_reference_image(template_image)
resampler.set_floating_image(attn_image)
resampler.set_padding_value(0)
resampler.set_interpolation_type_to_linear()
attn_image = resampler.forward(attn_image)
attn_image.write(working_folder + '/stir_mu_map_in_recon_space')


#%% create folders for results
path_sino = working_folder + '/sino'
path_rando = working_folder + '/rando'
path_recon = working_folder + '/recon'
if not os.path.exists(path_sino):
    os.makedirs(path_sino, mode=0o770)
    print('Create Folder: {}'.format(path_sino))
if not os.path.exists(path_rando):
    os.makedirs(path_rando, mode=0o770)
    print('Create Folder: {}'.format(path_rando))
if not os.path.exists(path_recon):
    os.makedirs(path_recon, mode=0o770)
    print('Create Folder: {}'.format(path_recon))


#%% define time frames (read frames.txt) and time intervals
# path to folder with frames.txt
frames_path = '/home/rich/tmp/processed/home/rich/Documents/ReCo/UKL_data/MedPhys_MoCo_Test_20190417_7/40004__30x20sec_4moco_PRR_AC_Images/'

with open(frames_path + 'frames.txt', 'r') as f:
    time_frames = [int(line.rstrip()) for line in f]

print('Frames: {}'.format(time_frames))

# define time intervals
time_intervals = sum_list(time_frames)
print('Time intervals: {}'.format(time_intervals))


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

# definitions for attenuation
#attn_image = Pet.ImageData(attn_file)
attn_acq_model = Pet.AcquisitionModelUsingRayTracingMatrix()
asm_attn = Pet.AcquisitionSensitivityModel(attn_image, attn_acq_model)

# definitions for detector sensitivity modelling
asm_norm = Pet.AcquisitionSensitivityModel(norm_file)
acq_model.set_acquisition_sensitivity(asm_norm)


#%% for-loop, reconstruction of time intervals
tprint('Start Recon')
for i in range(len(time_intervals)-1):
    print('Begin reconstruction: Frame {}'.format(i))

    # listmode-to-sinogram
    lm2sino.set_time_interval(time_intervals[i], time_intervals[i+1])
    lm2sino.set_up()
    lm2sino.process()
    acq_data = lm2sino.get_output()
    acq_data.write(working_folder + '/sino/sino'+str(i))
    initial_image = acq_data.create_uniform_image(1.0)
    initial_image.write(working_folder + '/uniform')

    # randoms estimate
    randoms = lm2sino.estimate_randoms()
    randoms.write(working_folder + '/rando/rando'+str(i))

    # reconstruct the data (includes all)
    obj_fun = Pet.make_Poisson_loglikelihood(acq_data)
    asm_attn.set_up(acq_data)
    attn_factors = Pet.AcquisitionData(acq_data)
    attn_factors.fill(1.0)
    asm_attn.unnormalise(attn_factors)
    asm_attn = Pet.AcquisitionSensitivityModel(attn_factors)
    asm = Pet.AcquisitionSensitivityModel(asm_norm, asm_attn)
    acq_model.set_acquisition_sensitivity(asm)
    acq_model.set_background_term(randoms)
    obj_fun.set_acquisition_model(acq_model)
    recon.set_objective_function(obj_fun)
    initial_image = acq_data.create_uniform_image(1.0)
    image = initial_image
    recon.set_up(image)

    recon.set_current_estimate(image)
    recon.process()

    # save recon images
    recon_image = recon.get_output()
    recon_image.write(working_folder + '/recon/recon' + str(i))

    # save Image as .nii
    recon_image = Reg.NiftiImageData(recon.get_output())
    recon_image.write(working_folder + '/floates/recon'+str(i))

    print('Reconstruction successful: Frame {}'.format(i))

tprint('Finish Recon')


#%% create folder for motion corrected images
path_moco = working_folder + '/moco/'
if not os.path.exists(path_moco):
    os.makedirs(path_moco, mode=0o770)
    print('Create Folder: {}'.format(path_moco))


#%% convert a array to a SIRF transformation matrix and then resample the float image
# define reference image and float-path
ref_file = working_folder + '/floates/recon0.nii'
ref = Eng_ref.ImageData(ref_file)
flo_path = working_folder + '/floates/'

# path to folder, contain the MCFLIRT matrices
path_mat_mcflirt = frames_path + "epi_frames.mat/"

tprint('Start Reg')

# for loop, simultaneous matrices and images
for mat, image in zip(sorted(os.listdir(path_mat_mcflirt)), sorted_alphanumeric(os.listdir(flo_path))):
    print('TM: {}, Float-Image: {}'.format(mat, image))

    # read tm-matrix as numpy array and read load float image
    matrix = numpy.loadtxt(path_mat_mcflirt + mat)
    flo = Eng_flo.ImageData(flo_path + image)

    # create affine transformation from numpy array
    tm = Reg.AffineTransformation(matrix)

    # motion information from another source and resample an image with this informations
    print('Begin resampling: {}'.format(image))
    resampler = Reg.NiftyResample()
    resampler.set_reference_image(ref)
    resampler.set_floating_image(flo)
    resampler.add_transformation(tm)
    resampler.set_padding_value(0)
    resampler.set_interpolation_type_to_linear()
    output = resampler.forward(flo)
    output.write(working_folder + '/moco/moco_'+str(mat))

    print('Resampling successful: {}'.format(image))

tprint('Finish Reg')


#%% define RTA method
# define initial image (first image, first frame)
initial = Reg.NiftiImageData(working_folder + '/moco/moco_MAT_0000.nii')
initial_array = initial.as_array()

# sum over all images (as array)
for image in sorted_alphanumeric(os.listdir(path_moco))[1:]:
    print(image)
    array = Reg.NiftiImageData(path_moco + image).as_array()
    initial_array += array

print('Final array (RTA): {}'.format(initial_array))

# create image
final_image = Reg.NiftiImageData(working_folder + '/moco/moco_MAT_0000.nii')
print(type(final_image))
final_image.fill(initial_array)
print(type(final_image))
final_image.write(working_folder + '/final_image_RTA.nii')

tprint('DONE!')


#%% plot image
# dataset must be an array
slice_num = initial_array.shape[0]//2

# plot figure
plt.figure()
imshow(initial_array[slice_num, :, :, ], [], 'MoCo - RTA')