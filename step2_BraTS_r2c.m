%% Simulate complex MRI data using BraTS2020 dataset with only amplitude images
% 2023.05.08 by Wang Yuyan
clc, clear, close all
addpath function
%% preset params
brats_folder = './demo_BraTS2020_niifile'; % enter BraTS file folder dir
% you can download BraTS data at http://braintumorsegmentation.org/
fastmri_folder = './Data/fastMRI'; % enter fastMRI file folder dir for phase map
cpx_folder = './Data/BraTS_complex'; % enter complex images file save dir
[folder_list, total_patient_num] = get_sub_folder(brats_folder);
fprintf('load .nii file from %s that has %s patients\n', [brats_folder, string(total_patient_num)])
patient_num = 5; % set the number of patient to be converted
patient_num = min(patient_num, total_patient_num);
lesion_size = 144; % minimum lesion size
brain_size = 1e4; % minimum brain size
slice_inter = 5; % interval to select frames
show_flag = 0; % wheter to print output figure
screen_size=get(0,'ScreenSize'); % set canvas size
screen_size(1:2)=screen_size(1:2)+0.05*screen_size(3:4);
screen_size(3:4)=screen_size(3:4)*0.8;
%% These data may cause problems
QC_flag = 1; % wether to skip 
QC_list = xlsread("BraTS2020_QC.xlsx","Sheet1","A:A")';
count = 1;
while count<patient_num+1 && QC_flag
    patient_id = split(folder_list(count),filesep);
    patient_id = patient_id(end);
    patient_no = split(patient_id,"_");
    patient_no = patient_no(end);
    if ismember(str2double(patient_no), QC_list) && QC_flag
        fprintf('{> patient %s has certain problem and skipped\n', patient_id);
        folder_list(count) = [];
        continue
    end
    count = count+1;
end
%% notice that BraTS2020's data contains 2019 & 2018, see name_mapping.csv in BraTS2020 folder for more details
for patient = 1:patient_num
    patient_loc = char(folder_list(patient));
    patient_id = split(folder_list(patient),filesep);
    patient_id = patient_id(end);
    % read nii from patient folder
    [T1_org, T1ce_org, T2_org, Flair_org, Tumor_mask_org] = read_nii(patient_loc); tic;
    if patient == 344; Tumor_mask_org(105:107,1:2,57:59) = 0; end  % 344's mask has a fault 
    %% area of lesion/tumor_core/enhanced/necrotic/brain, choose slice according to settings
    Tumor_mask_org = loc_brain(T2_org, Tumor_mask_org);
    Mask_area = cal_mask_area(Tumor_mask_org);
    Feasible_slice = find(Mask_area(1,:)>lesion_size & Mask_area(5,:)>brain_size);
    Feasible_slice = Feasible_slice(1:slice_inter:end);
    if isempty(Feasible_slice)
        fprintf('{> patient %s does not have enough slice\n', patient_id);
        continue
    else
        fprintf('{> patient %s processing with slice num %s\n', [patient_id, length(Feasible_slice)]);
    end
    T1_org = T1_org(:,:,Feasible_slice);
    T1ce_org = T1ce_org(:,:,Feasible_slice);
    T2_org = T2_org(:,:,Feasible_slice);
    Flair_org = Flair_org(:,:,Feasible_slice);
    Tumor_mask_org = Tumor_mask_org(:,:,Feasible_slice);
    Mask_area = cal_mask_area(Tumor_mask_org);
    show_slice = find(Mask_area(1,:) == max(Mask_area(1,:)), 1); % choose one slice with largest lesion to show
    %% create skin & skull with nature texture & record gap:6 skull:7 in mask
    fprintf('   create skull => ')
    [T1_sk, T1ce_sk, T2_sk, Flair_sk, Tumor_mask] = create_skull(T1_org, T1ce_org, T2_org, Flair_org, Tumor_mask_org);
    %% show amp image with added skin & skull
    if show_flag && mod(patient,10)==1
        figure('Name', patient_id+" amp image")
        amp_range=[0,1];
        for i=1:3
            s = min(max(show_slice+5*(i-2), 1), size(T2_sk,3));
            subplot(3,4,4*i-3), imshow(squeeze(T1ce_sk(:,:,s)),amp_range), colormap('gray'), title("T1ce slice "+s);
            subplot(3,4,4*i-2), imshow(squeeze(T2_sk(:,:,s)),amp_range), colormap('gray'), title("T2 slice "+s);
            subplot(3,4,4*i-1), imshow(squeeze(Flair_sk(:,:,s)),amp_range), colormap('gray'), title("Flair slice "+s);
            subplot(3,4,4*i), imshow(squeeze(Flair_sk(:,:,s)),[0,1]), colormap('gray');
            color_mask = cmp_brats(squeeze(Tumor_mask(:,:,s))); % tumor mask in color according to BraTS settings
            trans = ones([size(color_mask,1), size(color_mask,2)])*0.5;
            trans(sum(color_mask,3)==0) = 0;
            hold on; ov=imshow(color_mask);
            set(ov,'AlphaData',trans), title("Mask slice "+s);
            hold off;
        end
        set(gcf, 'position', screen_size);
    end
    %% add phase to get complex valued image
    fprintf('simulate phase => ')
    [T1, T1ce, T2, Flair] = real_to_complex(T1_sk, T1ce_sk, T2_sk, Flair_sk, Tumor_mask, fastmri_folder);
    %% show complex image with phase
    if show_flag && mod(patient,10)==1
        figure('Name',patient_id+" complex image")
        amp_range=[0,1];
        ri_range=[-0.5,0.5];
        T1ce_cpx = squeeze(T1ce(:,:,show_slice));
        subplot(3,4,1), imshow(abs(T1ce_cpx),amp_range),colorbar(),title("T1ce cpx amp");
        subplot(3,4,2), imshow(angle(T1ce_cpx),[-pi,pi]),colorbar(),title("T1ce cpx angle");
        subplot(3,4,3), imshow(real(T1ce_cpx),ri_range),colorbar(),title("T1ce cpx real");
        subplot(3,4,4), imshow(imag(T1ce_cpx),ri_range),colorbar(),title("T1ce cpx imag");
        T2_cpx = squeeze(T2(:,:,show_slice));
        subplot(3,4,5), imshow(abs(T2_cpx),amp_range),colorbar(),title("T2 cpx amp");
        subplot(3,4,6), imshow(angle(T2_cpx),[-pi,pi]),colorbar(),title("T2 cpx angle");
        subplot(3,4,7), imshow(real(T2_cpx),ri_range),colorbar(),title("T2 cpx real");
        subplot(3,4,8), imshow(imag(T2_cpx),ri_range),colorbar(),title("T2 cpx imag");
        Flair_cpx = squeeze(Flair(:,:,show_slice));
        subplot(3,4,9), imshow(abs(Flair_cpx),amp_range),colorbar(),title("Flair cpx amp");
        subplot(3,4,10), imshow(angle(Flair_cpx),[-pi,pi]),colorbar(),title("Flair cpx angle");
        subplot(3,4,11), imshow(real(Flair_cpx),ri_range),colorbar(),title("Flair cpx real");
        subplot(3,4,12), imshow(imag(Flair_cpx),ri_range),colorbar(),title("Flair cpx imag");
        set(gcf, 'position', screen_size);
    end
    %% save complex data to .mat
    if ~exist(fullfile(cpx_folder,patient_id),'dir')
        mkdir(fullfile(cpx_folder,patient_id))
    end
    save(fullfile(cpx_folder,patient_id,patient_id+"_t1.mat"),"T1")
    save(fullfile(cpx_folder,patient_id,patient_id+"_t1ce.mat"),"T1ce")
    save(fullfile(cpx_folder,patient_id,patient_id+"_t2.mat"),"T2")
    save(fullfile(cpx_folder,patient_id,patient_id+"_flair.mat"),"Flair")
    save(fullfile(cpx_folder,patient_id,patient_id+"_seg.mat"),"Tumor_mask")
    fprintf('done\n')
    fprintf('   time %s sec\n', string(toc));
end
