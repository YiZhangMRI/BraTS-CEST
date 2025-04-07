%% Pre-process fastMRI dataset to generate CSM & combined img for phase map
% 2023.10.25 by Wang Yuyan
clc, clear, close all
addpath function
%% preset params
h5_folder = './demo_fastMRI_h5file'; % enter fastMRI h5 file folder dir
% you can download fastMRI data at https://fastmri.med.nyu.edu/
h5_folder_list = dir(fullfile(h5_folder,'*.h5'));
full_csm = "full"; % choose "full" or "partial" for region of csm
mat_folder = './Data/fastMRI'; % enter fastMRI mat file folder save dir
total_data_num = length(h5_folder_list);
data_num = 2; % set the number of data to be converted
data_num = min(data_num, total_data_num);
fprintf('load .h5 file from %s that has %s datas\n', [h5_folder, string(total_data_num)])
img_size = [96, 96]; % resized width & hight
slice_num = 10; % num of slices extracted from rawdata
coil_num = 16; % num of coils extracted from rawdata
%% read each data and convert to .mat
for d=1:data_num
    data_dir = fullfile(h5_folder_list(d).folder, h5_folder_list(d).name);
    ksp = h5read(data_dir,'/kspace');
    ksp = ksp.r + 1i*ksp.i;
    if size(ksp,3)<coil_num || size(ksp,4)<slice_num
        fprintf('{> data %s does not have enough coil or slice\n', string(d));
        continue
    else
        fprintf('{> data %s processing\n   ', string(d));tic;
    end
    %% resize img to preset sizes and calculate csm
    image = fft2c_mri(ksp);
    x_num = size(image,1);
    y_num = size(image,2);
    image = image(:, y_num*0.5-x_num*0.5+1:y_num*0.5+x_num*0.5, :, :); % skip sampling on ky-axis
    for s=1:slice_num
        fprintf('slice %s => ', string(s));
        image_slice = squeeze(image(:,:,:,s));
        img_coil_resize = zeros(img_size(1),img_size(2),coil_num);
        for c=1:coil_num
            img_coil = image_slice(:,:,c);
            img_coil = imresize(img_coil, img_size);
            img_coil = rot90(img_coil);
            img_coil_resize(:,:,c) = img_coil;
        end
        [sensitivities, reference] = cal_csm(img_coil_resize, full_csm);
        save_path = fullfile(mat_folder, num2str(d,"%04d"));
        if ~exist(save_path,'dir')
            mkdir(save_path);
        end
        save(fullfile(save_path, "espirit"+num2str(s,"%02d")+".mat"), "sensitivities", "reference");
    end
    fprintf('done\n')
    fprintf('   time %s sec\n', string(toc));
end
