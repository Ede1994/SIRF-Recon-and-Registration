#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright University College London Hospital, 2020
Author: Eric Einspaenner, Institute of Nuclear Medicine
For internal research only.
"""
import os
import re
import shutil
import numpy
from art import tprint
import sirf.STIR as Pet
import sirf.Reg as Reg
import sirf.Reg as Eng_ref
import sirf.Reg as Eng_flo


#%% Functions

# generate a list, contains the time intervalls (sum of frames)
def sum_list(lists):
    sum_liste = []
    add = int(0)
    sum_liste.append(add)
    for frame in lists[0:]:
        add += frame
        sum_liste.append(add)
    return sum_liste

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


#%% Settings for reconstruction

# set acq_model
acq_model = Pet.AcquisitionModelUsingRayTracingMatrix()
acq_model.set_num_tangential_LORs(5)

# set recon, OSEM
recon = Pet.OSMAPOSLReconstructor()
num_subsets = 21
num_subiterations = 63
recon.set_num_subsets(num_subsets)
recon.set_num_subiterations(num_subiterations)


#%% redirect STIR messages to some files
# you can check these if things go wrong

msg_red = Pet.MessageRedirector('info.txt', 'warn.txt')


#%% main
data_path = '/home/eric/Dokumente/PythonProjects/SIRF/UKL_data/LM/20190723/20190723/'

list_LM_files = [f for f in os.listdir(data_path) if f.endswith(".l.hdr")]
list_norm_files = [f for f in os.listdir(data_path) if f.endswith(".n.hdr")]

print(sorted_alphanumeric(list_LM_files))
print(sorted_alphanumeric(list_norm_files))

j = 0
for list_file, norm_file in zip(sorted_alphanumeric(list_LM_files), sorted_alphanumeric(list_norm_files)):

    # Definition of different path
    py_path = os.getcwd()
    print(py_path)

    # avoid cluttering of files, delete working-folder and start from scratch
    working_folder = py_path + '/working'
    if os.path.exists(working_folder):
        shutil.rmtree(working_folder)
    if not os.path.exists(working_folder):
        os.makedirs(working_folder, mode=0o770)

    # change the current working directory to the given path
    os.chdir(working_folder)

    # input files
    list_file = data_path + list_file
    norm_file = data_path + norm_file
    print('LM data: {}'.format(list_file))
    print('Norm data: {}'.format(norm_file))
    
    attn_file = '/home/eric/Dokumente/PythonProjects/SIRF/UKL_data/mu_Map/stir_mu_map.hv'
    print('mu-Map: {}'.format(attn_file))

    attn_image = Pet.ImageData(attn_file)

    # output filename prefixes
    sino_file = 'sino'


    ### Create folders for results ###
    path_sino = working_folder + '/sino/'
    path_rando = working_folder + '/rando/'
    path_NAC = working_folder + '/recon/NAC/'
    path_tm = working_folder + '/tm/'
    path_mu = working_folder + '/mu/'
    path_AC = working_folder + '/recon/AC/'
    path_moco = working_folder + '/moco/'

    if not os.path.exists(path_sino):
        os.makedirs(path_sino, mode=0o770)
        print('Create Folder: {}'.format(path_sino))
    if not os.path.exists(path_rando):
        os.makedirs(path_rando, mode=0o770)
        print('Create Folder: {}'.format(path_rando))
    if not os.path.exists(path_NAC):
        os.makedirs(path_NAC, mode=0o770)
        print('Create Folder: {}'.format(path_NAC))
    if not os.path.exists(path_tm):
        os.makedirs(path_tm, mode=0o770)
        print('Create Folder: {}'.format(path_tm))
    if not os.path.exists(path_mu):
        os.makedirs(path_mu, mode=0o770)
        print('Create Folder: {}'.format(path_mu))
    if not os.path.exists(path_AC):
        os.makedirs(path_AC, mode=0o770)
        print('Create Folder: {}'.format(path_AC))
    if not os.path.exists(path_moco):
        os.makedirs(path_moco, mode=0o770)
        print('Create Folder: {}'.format(path_moco))
    
    ### create template and set lm2sino converter ###
    # template for acq_data
    template_acq_data = Pet.AcquisitionData('Siemens_mMR', span=11, max_ring_diff=16, view_mash_factor=1)
    template_acq_data.write('template.hs')
   
    ### Create listmode-to-sinograms converter object ###
    lm2sino = Pet.ListmodeToSinograms()

    ### set input, output and template files ###
    lm2sino.set_input(list_file)
    lm2sino.set_output_prefix(sino_file)
    lm2sino.set_template('template.hs')


    ### Define time frames (read frames.txt) and time intervals ###    
    # path to folder with frames.txt
    frames_path = '/home/eric/Dokumente/PythonProjects/SIRF/UKL_data/frames/'

    # read frames.txt
    with open(frames_path + 'frames.txt', 'r') as f:
        time_frames = [int(line.rstrip()) for line in f]
        print('Frames: {}'.format(time_frames))

    # number of motion steps
    num_motion_steps = len(time_frames)

    # define time intervals
    time_intervals = sum_list(time_frames)
    print('Time intervals: {}'.format(time_intervals))


    ### NAC reconstruction ###   
    # definitions for detector sensitivity modelling
    asm_norm = Pet.AcquisitionSensitivityModel(norm_file)
    acq_model.set_acquisition_sensitivity(asm_norm)

    tprint('Start NAC Recon')

    for i in range(len(time_intervals)-1):
        print('Begin: Frame {}'.format(i))
        print('Time interval: {} - {}'.format(time_intervals[i], time_intervals[i+1]))

        # listmode-to-sinogram
        lm2sino.set_time_interval(time_intervals[i], time_intervals[i+1])
        lm2sino.set_up()
        lm2sino.process()
        acq_data = lm2sino.get_output()
        acq_data.write(path_sino + 'sino'+str(i))

        # randoms estimate
        randoms = lm2sino.estimate_randoms()
        randoms.write(path_rando + 'rando' + str(i))
        
        # reconstruct the data (without mu-map)
        obj_fun = Pet.make_Poisson_loglikelihood(acq_data)
        acq_model.set_background_term(randoms)
        recon.set_objective_function(obj_fun)
        initial_image = acq_data.create_uniform_image(1.0)
        image = initial_image
        recon.set_up(image)

        recon.set_current_estimate(image)
        recon.process()

        # save recon images
        recon_image = recon.get_output()
        recon_image.write(path_NAC + 'NAC_' + str(i))

        # save Image as .nii
        recon_image = Reg.NiftiImageData(recon.get_output())
        recon_image.write(path_NAC + 'NAC_' + str(i))

        print('Reconstruction successful: Frame {}'.format(i))

    tprint('Finish NAC Recon')
    
    ### NAC registration ###
    # Registration NAC, delivers transformation matrices 
    # define reference image (first image) and float-path

    tprint('Start Registration of NACs')

    # refernce file
    ref_file = path_NAC + 'NAC_0.nii'
    ref = Eng_ref.ImageData(ref_file)

    # float files
    flo_path = path_NAC
    list_NACs = [f for f in os.listdir(path_NAC) if f.endswith(".nii")]

    # Niftyreg with NAC images
    reg_nac(ref, flo_path, list_NACs)

    tprint('Finish Registration')
    
    ### resample mu-Map into each NAC space ###
    ref_file = path_NAC + 'NAC_0.nii'
    ref = Eng_ref.ImageData(ref_file)
    flo = Eng_flo.ImageData(attn_image)

    tprint('Start Resampling')

    for i in range(30):
        print('Begin resampling mu-Maps: {}, with {}'.format(i, path_tm + 'tm_nac_inv_' + str(i) + '.txt'))
        tm = numpy.loadtxt(path_tm + 'tm_nac_inv_' + str(i) + '.txt')
        tm_fwd = Reg.AffineTransformation(tm)

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
    
    ### List of sino and rand ###
    list_sino = [f for f in os.listdir(working_folder + '/sino/') if f.endswith(".hs")]
    list_rando = [f for f in os.listdir(working_folder + '/rando/') if f.endswith(".hs")]
    list_mu = [f for f in os.listdir(path_mu) if f.endswith(".nii")]


    ### AC reconstruction ###
    tprint('Start AC Recon')

    for i, sino, random, mu in zip(range(len(path_sino)), sorted_alphanumeric(list_sino), sorted_alphanumeric(list_rando), sorted_alphanumeric(list_mu)):

        sino_pet = Pet.AcquisitionData(path_sino + sino)
        print(sino)
        randoms_pet = Pet.AcquisitionData(path_rando + random)
        print(random)
        mu_pet = Pet.ImageData(path_mu + mu)
        print(mu)

        # definitions for attenuation
        attn_acq_model = Pet.AcquisitionModelUsingRayTracingMatrix()
        asm_attn = Pet.AcquisitionSensitivityModel(mu_pet, attn_acq_model)

        # reconstruct the data (includes all)
        obj_fun = Pet.make_Poisson_loglikelihood(sino_pet)
        asm_attn.set_up(sino_pet)
        attn_factors = Pet.AcquisitionData(sino_pet)
        attn_factors.fill(1.0)
        asm_attn.unnormalise(attn_factors)
        asm_attn = Pet.AcquisitionSensitivityModel(attn_factors)
        asm = Pet.AcquisitionSensitivityModel(asm_norm, asm_attn)
        acq_model.set_acquisition_sensitivity(asm)
        acq_model.set_background_term(randoms_pet)
        obj_fun.set_acquisition_model(acq_model)
        recon.set_objective_function(obj_fun)
        initial_image = sino_pet.create_uniform_image(1.0)
        image = initial_image
        recon.set_up(image)

        recon.set_current_estimate(image)
        recon.process()

        # save recon images
        recon_image = recon.get_output()
        recon_image.write(path_AC + 'AC_' + str(i))

        # save Image as .nii
        recon_image = Reg.NiftiImageData(recon.get_output())
        recon_image.write(path_AC + 'AC_' + str(i))

        print('Reconstruction successful: Frame {}'.format(i))

    tprint('Finish AC Recon')
    
    ### resample the float images back to reference ###
    # list of all ACs
    list_ACs = [f for f in os.listdir(path_AC) if f.endswith(".nii")]

    # define reference image and float-path
    ref_file = path_AC + 'AC_0.nii'
    ref = Eng_ref.ImageData(ref_file)
    flo_path = path_AC

    # settings for image resampler
    resampler_im = Reg.NiftyResample()
    resampler_im.set_reference_image(ref)
    resampler_im.set_padding_value(0)
    resampler_im.set_interpolation_type_to_linear()

    tprint('Start Resampling')

    # for loop, simultaneous matrices and images
    for num, image in zip(range(len(list_ACs)), sorted_alphanumeric(list_ACs)):
        print('TM: {}, Float-Image: {}'.format('tm_nac_' + str(num) + '.txt', image))

        flo = Eng_flo.ImageData(path_AC + image)

        # read tm-matrix as numpy array
        matrix = numpy.loadtxt(path_tm + 'tm_nac_' + str(num) + '.txt')
  
        # create affine transformation from numpy array
        tm = Reg.AffineTransformation(matrix)

        # motion information from another source and resample an image with this informations
        print('Begin resampling: {}'.format(image))
        resampler_im.clear_transformations()
        resampler_im.set_floating_image(flo)
        resampler_im.add_transformation(tm)
        new_im = resampler_im.forward(flo)
        new_im.write(path_moco + 'moco_' + str(num))

        print('Resampling successful: {}'.format(image))

    tprint('Finish Resampling')


    ### define RTA method ##
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
    final_image.write('/home/eric/Dokumente/PythonProjects/SIRF/final_image_RTA' + str(j) + '.nii')

    tprint('DONE!')
    
    j += 1