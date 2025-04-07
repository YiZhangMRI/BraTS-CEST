%% Simulate CEST images using complexed BraTS images, according to T1/T2w images & BM simulation
% 2023.06.05 by Wang Yuyan
addpath function
clc, clear, close all
%% preset params
cpx_folder='./Data/BraTS_complex'; % enter complex images file save dir
[cpx_folder_list, total_patient_num] = get_sub_folder(cpx_folder);
fprintf('load complex images from %s that has %s patients\n', [cpx_folder, string(total_patient_num)])
CEST_folder='./demo_BraTS_CEST'; % enter CEST images file save dir
full_csm = "full"; % choose "full" or "partial" for region of csm
fastmri_folder = './Data/fastMRI'; % enter fastMRI file folder dir for csm
[fastmri_folder_list,~] = get_sub_folder(fastmri_folder);
CEST_imgsize = [96, 96]; % size for generated CEST img
patient_num = 5; % set the number of patient to be converted
patient_num = min(patient_num, total_patient_num);
slice_inter = 1; % interval to select frames
enhance_flag = 1; % choose 1 or 0 for rigid body transform enhancement
shift_distance = 20; % max translation distance from the edge
T2blur_flag = 1; % choose 1 or 0 for TSE T2 blur
motion_flag = 0; % choose 1 or 0 for motion artifact
show_flag = 0; % choose 1 or 0 for fig exhibition
screen_size=get(0,'ScreenSize');
screen_size(1:2)=screen_size(1:2)+0.05*screen_size(3:4);
screen_size(3:4)=screen_size(3:4)*0.8;
%% read scan freq
scanned_offset = transpose(load('Offset frequencies for APT scans - 54pts.txt'));
scanned_offset = round(scanned_offset/1.2314)/100;
scanned_offset_sorted = sort(scanned_offset,'descend');
[scanned_offset_uniq, col] = unique(scanned_offset_sorted, "stable");
%% These data may cause problems
QC_flag = 1; % wether to skip 
QC_list = xlsread("BraTS2020_QC.xlsx","Sheet1","A:A")';
count = 1;
while count<patient_num+1 && QC_flag
    patient_id = split(cpx_folder_list(count),filesep);
    patient_id = patient_id(end);
    patient_no = split(patient_id,"_");
    patient_no = patient_no(end);
    if ismember(str2double(patient_no), QC_list) && QC_flag
        fprintf('{> patient %s has certain problem and skipped\n', patient_id);
        cpx_folder_list(count) = [];
        continue
    end
    count = count+1;
end
%% read patient data
for patient= 1:patient_num
    patient_loc = cpx_folder_list(patient);
    patient_id = split(cpx_folder_list(patient),filesep);
    patient_id = patient_id(end);
    %% load complex valuled data
    load(fullfile(patient_loc, patient_id+"_t1.mat"));
    load(fullfile(patient_loc, patient_id+"_t1ce.mat"));
    load(fullfile(patient_loc, patient_id+"_t2.mat")); 
    load(fullfile(patient_loc, patient_id+"_flair.mat"));
    load(fullfile(patient_loc, patient_id+"_seg.mat"));
    fprintf('{> patient %s converting to CEST img', patient_id); tic;
    %% select slices at intervals since the adjacent slices are similar
    Feasible_slice=1:slice_inter:size(Tumor_mask,3);
    fprintf('with slice num %s\n', string(length(Feasible_slice)));
    T1 = T1(:,:,Feasible_slice);
    T1ce = T1ce(:,:,Feasible_slice);
    T2 = T2(:,:,Feasible_slice);
    Flair = Flair(:,:,Feasible_slice);
    Tumor_mask = Tumor_mask(:,:,Feasible_slice);
    Tumor_area = cal_mask_area(Tumor_mask);
    show_slice = find(Tumor_area(1,:)==max(Tumor_area(1,:)), 1);
    %% show selected slices
    if show_flag && mod(patient,10)==1
        fig1=figure('Name',patient_id+" amp image");
        amp_range=[0,1];
        for i=1:3
            s = min(max(show_slice+5*(i-2),1),size(T2,3));
            subplot(3,4,4*i-3), imshow(abs(squeeze(T1(:,:,s))),amp_range), colormap('gray'),title("T1ce slice "+s);
            subplot(3,4,4*i-2), imshow(abs(squeeze(T2(:,:,s))),amp_range), colormap('gray'),title("T2 slice "+s);
            subplot(3,4,4*i-1), imshow(abs(squeeze(Flair(:,:,s))),amp_range), colormap('gray'),title("Flair slice "+s);
            subplot(3,4,4*i), imshow(abs(squeeze(Flair(:,:,s))),[0,1]), colormap('gray');
            color_mask = cmp_brats(squeeze(Tumor_mask(:,:,s)));
            trans = ones([size(color_mask,1),size(color_mask,2)])*0.5;
            trans(sum(color_mask,3)==0) = 0;
            hold on; ov=imshow(color_mask);
            set(ov,'AlphaData',trans),title("Mask slice "+s);
            hold off;
        end
        set(gcf, 'position', screen_size);
    end
    %% divide image into normal/tumor/edema/necrotic/skull/skin & calculate T1/T2 grade maps through discretization
    Division = zone_division(Tumor_mask);
    tissue = ["normal", "tumor", "edema", "necro", "skull", "skin"];
    t_level = [9, 9, 9, 4, 4, 4]; % num of relaxation time grades 
    t1_lb = [0.25, 0.1, 0.1, 0.05, 0.1, 0.1]; % lower bound of T1
    t1_ub = [0.9, 0.9, 0.95, 0.25, 0.95, 0.9]; % upper bound of T1
    t2_lb = [0.1, 0.1, 0.1, 0.8, 0.1, 0.1]; % lower bound of T2
    t2_ub = [0.8, 0.8, 0.9, 0.95, 0.95, 0.5]; % upper bound of T2
    t_edge = ["edge", "noedge", "noedge", "noedge", "edge", "edge"]; % whether to excluded edges
    T1_divzone = T1.*Division;
    T2_divzone = T2.*Division;
    T1_dismap = zeros(size(T1_divzone));
    T2_dismap = zeros(size(T2_divzone));
    fprintf('   T1/T2 mapping => ')
    for zone=1:length(tissue)
        for s=1:size(T1_divzone,3)
            % T1 discretization for estimated T1 mapping, lower amp correspondences to higher T1
            T1_zone = abs(squeeze(T1_divzone(:,:,s,zone)));
            T1_roi = squeeze(Division(:,:,s,zone));
            T1_dis_config.levels  = t_level(zone);
            T1_dis_config.lb      = t1_lb(zone);
            T1_dis_config.ub      = t1_ub(zone);
            T1_dis_config.sort    = "descend";
            T1_dis_config.edge    = t_edge(zone);
            T1_dismap(:,:,s,zone) = amp_discretization(T1_zone, T1_roi, T1_dis_config);
            % T2 discretization for estimated T2 mapping, higher amp correspondences to higher T2
            T2_zone = abs(squeeze(T2_divzone(:,:,s,zone)));
            T2_roi = squeeze(Division(:,:,s,zone));
            T2_dis_config.levels  = t_level(zone);
            T2_dis_config.lb      = t2_lb(zone);
            T2_dis_config.ub      = t2_ub(zone);
            T2_dis_config.sort    = "ascend";
            T2_dis_config.edge    = t_edge(zone);
            T2_dismap(:,:,s,zone) = amp_discretization(T2_zone, T2_roi, T2_dis_config);
        end
    end
    %% generate Z map for each region according to T1/T2 grade maps
    fprintf('generate Z map => ')
    Z_image_brain = zeros([size(T1),length(scanned_offset)]);
    Z_image_bg = zeros([size(T1),length(scanned_offset)]);
    parfor s=1:size(T1_divzone,3)
        Z_image_brain_slice = zeros([size(T1,1),size(T1,2),length(scanned_offset),length(tissue)]);
        for zone=1:length(tissue)
           T1_zone = squeeze(T1_dismap(:,:,s,zone));
           T2_zone = squeeze(T2_dismap(:,:,s,zone));
           roi = squeeze(Division(:,:,s,zone));
           Z_image_brain_slice(:,:,:,zone) = zsp_map(T1_zone, T2_zone, roi, t_level(zone), tissue(zone));
        end
        Z_image_brain_slice = sum(Z_image_brain_slice,4);
        Z_image_brain(:,:,s,:) = Z_image_brain_slice;
        T1_bg = abs(squeeze(T1_divzone(:,:,s,7)));
        T2_bg = abs(squeeze(T2_divzone(:,:,s,7)));
        roi_bg = squeeze(Division(:,:,s,7));
        Z_image_bg(:,:,s,:) = zsp_map_bg(T1_bg, T2_bg, roi_bg, Z_image_brain_slice);
    end
    % smooth junction between different tissues
    Z_image = Z_image_brain + Z_image_bg;
    Z_image = zimage_smooth(Z_image, Division, 2);
    if show_flag && mod(patient,10)==1
        fig2=figure('Name',patient_id+" Z map");
        for f=1:32
            subplot(4,8,f), imshow(abs(squeeze(Z_image(:,:,show_slice,f))),[0,1]),colormap("jet")
        end
        set(gcf, 'position', screen_size);
    end
    %% get CEST image for T1/T2/Flair
    fprintf('synthesize CEST image => ')
    CEST_image_T1 = T1.*Z_image;
    CEST_image_T2 = T2.*Z_image;
    CEST_image_Flair = Flair.*Z_image;
    %% randomly crop CEST img for enhancement at each slice
    if enhance_flag
        fprintf('rigid trans enhancement => ')
        brain_roi = sum(Division(:,:,:,1:length(tissue)),4);
        CEST_image_T1 = random_crop(CEST_image_T1, brain_roi, shift_distance);
        CEST_image_T2 = random_crop(CEST_image_T2, brain_roi, shift_distance);
        CEST_image_Flair = random_crop(CEST_image_Flair, brain_roi, shift_distance);
    end
    %% resize CEST img to preset size
    CEST_image_T1 = imresize(CEST_image_T1, CEST_imgsize, 'nearest');
    CEST_image_T2 = imresize(CEST_image_T2, CEST_imgsize, 'nearest');
    CEST_image_Flair = imresize(CEST_image_Flair, CEST_imgsize, 'nearest');
    %% add extra noise, T2 blur or motion effect
    fprintf('add noise => ')
    noise_level = 0.02*(1+(rand()-0.5)*0.5);
    T2decay_filter = TSE_pss(0.25*(1+(rand()-0.5)*0.2), CEST_imgsize);
    for s=1:size(T1_divzone,3)
        CEST_image_T1(:,:,s,:) = ksp_add_noise(squeeze(CEST_image_T1(:,:,s,:)), noise_level);
        CEST_image_T2(:,:,s,:) = ksp_add_noise(squeeze(CEST_image_T2(:,:,s,:)), noise_level);
        CEST_image_Flair(:,:,s,:) = ksp_add_noise(squeeze(CEST_image_Flair(:,:,s,:)), noise_level);
        if motion_flag
            % add alias and extra noise around 0ppm
            alias_level = 0.1*(1+(rand()-0.5)*0.5);
            CEST_image_T1(:,:,s,:) = ksp_add_alias(squeeze(CEST_image_T1(:,:,s,:)), alias_level);
            CEST_image_T2(:,:,s,:) = ksp_add_alias(squeeze(CEST_image_T2(:,:,s,:)), alias_level);
            CEST_image_Flair(:,:,s,:) = ksp_add_alias(squeeze(CEST_image_Flair(:,:,s,:)), alias_level);
        elseif T2blur_flag
            % add T2 decay filtering caused by TSE
            CEST_image_T1(:,:,s,:) = T2_blur(squeeze(CEST_image_T1(:,:,s,:)), T2decay_filter);
            CEST_image_T2(:,:,s,:) = T2_blur(squeeze(CEST_image_T2(:,:,s,:)), T2decay_filter);
            CEST_image_Flair(:,:,s,:) = T2_blur(squeeze(CEST_image_Flair(:,:,s,:)), T2decay_filter);
        end
    end
    %% check simulated CEST image
    if show_flag && mod(patient,10)==1
        fig3=figure('Name',patient_id+" CEST image");
        for f=1:6
            subplot(3,6,f), imshow(abs(squeeze(CEST_image_T1(:,:,show_slice,(f-1)*9+1))),[0,1]),colormap("jet")
            title(num2str(scanned_offset_sorted((f-1)*9+1),'%.1f')+" ppm")
            subplot(3,6,f+6), imshow(abs(squeeze(CEST_image_T2(:,:,show_slice,(f-1)*9+1))),[0,1]),colormap("jet")
            subplot(3,6,f+12), imshow(abs(squeeze(CEST_image_Flair(:,:,show_slice,(f-1)*9+1))),[0,0.8]),colormap("jet")
        end
        set(gcf, 'position', screen_size);
        CEST_image = squeeze(CEST_image_Flair(:,:,show_slice,:));
        ksp_CEST = fft2c_mri(squeeze(CEST_image));
        fig4=figure('Name',patient_id+" normalized CEST image");
        ratio = cal_ratio(CEST_image);
        for f=1:32
            slice = col(f); subplot(4,8,f), imshow(abs(squeeze(CEST_image(:,:,slice)))/ratio(slice),[0,0.8]),colormap("jet")
            title(num2str(scanned_offset_uniq(f),'%.1f')+" ppm")
        end
        set(gcf, 'position', screen_size);
        fig5=figure('Name',patient_id+" normalized CEST kspace");
        for f=1:32
            slice = col(f); subplot(4,8,f), imshow(abs(squeeze(ksp_CEST(:,:,slice)))/ratio(slice),[0,0.25]),colormap("jet")
            title(num2str(scanned_offset_uniq(f),'%.1f')+" ppm")
        end
        set(gcf, 'position', screen_size);
        fig6=figure('Name',patient_id+" final Z map");
        for f=1:32
            slice = col(f); subplot(4,8,f), imshow(abs(squeeze(CEST_image(:,:,slice)))./abs(CEST_image(:,:,1)),[0,1]),colormap("jet")
            title(num2str(scanned_offset_uniq(f),'%.1f')+" ppm")
        end
        set(gcf, 'position', screen_size);
        fig7=figure('Name',patient_id+" APTw image");
        apt_neg = squeeze(CEST_image(:,:,47))./squeeze(CEST_image(:,:,1));
        apt_pos = squeeze(CEST_image(:,:,9))./squeeze(CEST_image(:,:,1));
        rainbow_map = rainbow('rainbow_idl');
        subplot(1,3,1), imshow(abs(apt_pos),[0,1]),colormap("jet"),colorbar(),title("-3.5 ppm")
        subplot(1,3,2), imshow(abs(apt_neg),[0,1]),colormap("jet"),colorbar(),title("3.5 ppm")
        subplot(1,3,3), imshow(abs(apt_neg)-abs(apt_pos),[-0.05,0.05]),colormap(gca,rainbow_map),colorbar(),title("APTw")
        set(gcf, 'position', screen_size);
    end
    %% assign CSM from fastMRI dataset
    fprintf('assign csm => ')
    save_id = split(patient_id,"_");
    save_id = save_id(3);
    T1_dir = fullfile(CEST_folder, "T1_"+full_csm, save_id);
    T2_dir = fullfile(CEST_folder, "T2_"+full_csm, save_id);
    Flair_dir = fullfile(CEST_folder, "FLAIR_"+full_csm, save_id);
    if ~exist(T1_dir,'dir'); mkdir(T1_dir); mkdir(T2_dir); mkdir(Flair_dir); end
    for s=1:size(T1_divzone,3)
        csm_org = load_fastmri_csm(CEST_imgsize, fastmri_folder_list);
        [csm_T1, ref_T1] = assign_csm(squeeze(CEST_image_T1(:,:,s,1)), csm_org, full_csm);
        [csm_T2, ref_T2] = assign_csm(squeeze(CEST_image_T2(:,:,s,1)), csm_org, full_csm);
        [csm_Flair, ref_Flair] = assign_csm(squeeze(CEST_image_Flair(:,:,s,1)), csm_org, full_csm);
    %% save simulated CEST image & CSM
        APT_data = squeeze(CEST_image_T1(:,:,s,:)); save(fullfile(T1_dir, "CEST"+num2str(s,"%02d")+".mat"), "APT_data");
        APT_data = squeeze(CEST_image_T2(:,:,s,:)); save(fullfile(T2_dir, "CEST"+num2str(s,"%02d")+".mat"), "APT_data");
        APT_data = squeeze(CEST_image_Flair(:,:,s,:)); save(fullfile(Flair_dir, "CEST"+num2str(s,"%02d")+".mat"), "APT_data");
        sensitivities = csm_T1; reference = ref_T1;
        save(fullfile(T1_dir, "espirit"+num2str(s,"%02d")+".mat"), "sensitivities", "reference");
        sensitivities = csm_T2; reference = ref_T2;
        save(fullfile(T2_dir, "espirit"+num2str(s,"%02d")+".mat"), "sensitivities", "reference");
        sensitivities = csm_Flair; reference = ref_Flair;
        save(fullfile(Flair_dir, "espirit"+num2str(s,"%02d")+".mat"), "sensitivities", "reference");
        if show_flag && s==show_slice && mod(patient,10)==1
            fig8=figure('Name',patient_id+" csm");
            for c=1:4
                subplot(2,4,c),imshow(abs(squeeze(csm_org(:,:,c))),[]),title("org coil "+string(c))
                subplot(2,4,c+4),imshow(abs(squeeze(csm_Flair(:,:,c))),[]),title("new coil "+string(c))
            end
            set(gcf, 'position', screen_size);
        end   
    end
    fprintf('done\n')
    fprintf('   time %s sec\n', string(toc));
end
