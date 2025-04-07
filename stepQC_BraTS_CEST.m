%% Quality control of BraTS_CEST, check simulated results
% 2023.06.13 by Wang Yuyan
addpath function
clc, clear, close all
%% preset params
full_csm = "full";
CEST_folder='./demo_BraTS_CEST'; % enter CEST images file save dir
cpx_folder = './Data/BraTS_complex'; % enter complex images file save dir
T1_folder=fullfile(CEST_folder, char("T1_"+full_csm));
T2_folder=fullfile(CEST_folder, char("T2_"+full_csm));
Flair_folder=fullfile(CEST_folder, char("FLAIR_"+full_csm));
[T1_folder_list, ~] = get_sub_folder(T1_folder);
[T2_folder_list, ~] = get_sub_folder(T2_folder);
[Flair_folder_list, CEST_num] = get_sub_folder(Flair_folder);
fprintf('quality control for CEST images in %s that has %s patients\n', [CEST_folder, string(CEST_num)])
[cpx_folder_list, ~] = get_sub_folder(cpx_folder);
close_flag = 1; % whether to close previous fig
pause_time = 1; % pause time between figs
show_modality = "Flair"; % which modality for more ditails T1/T2/Flair
slice_inter = 1; % interval to select frames
screen_size=get(0,'ScreenSize');
screen_size(1:2)=screen_size(1:2)+0.05*screen_size(3:4);
screen_size(3:4)=screen_size(3:4)*0.8;
has_last_fig=0;
%% read scan freq
scanned_offset = transpose(load(fullfile('BM_simu_Zspetrum','Offset frequencies for APT scans - 54pts.txt')));
scanned_offset = round(scanned_offset/1.2314)/100;
scanned_offset_sorted = sort(scanned_offset,'descend');
[scanned_offset_uniq, col] = unique(scanned_offset_sorted, "stable");
%% read patient data
for d=1:CEST_num
    patient_num = split(T1_folder_list(d),filesep); patient_num = str2double(patient_num(end));
    patient_id = split(cpx_folder_list(d),filesep); patient_id = patient_id(end);
    %% locate the best view slice
    load(fullfile(cpx_folder_list(d), patient_id+"_seg.mat"));
    Tumor_area = cal_mask_area(Tumor_mask);
    show_slice = find(Tumor_area(1,:)==max(Tumor_area(1,:)), 1);
    CEST_slice = 1:slice_inter:size(Tumor_mask,3);
    show_slice_loc = abs(CEST_slice/show_slice-1);
    show_CEST_slice = find(show_slice_loc==min(show_slice_loc));
    %% load data
    load(fullfile(T1_folder_list(d), "CEST"+num2str(show_CEST_slice,"%02d")+".mat")); CEST_image_T1 = APT_data;
    load(fullfile(T2_folder_list(d), "CEST"+num2str(show_CEST_slice,"%02d")+".mat")); CEST_image_T2 = APT_data;
    load(fullfile(Flair_folder_list(d), "CEST"+num2str(show_CEST_slice,"%02d")+".mat")); CEST_image_Flair = APT_data;
    load(fullfile(Flair_folder_list(d), "espirit"+num2str(show_CEST_slice,"%02d")+".mat")); CSM = sensitivities;
    fprintf('{> patient %s slice %s simulated CEST image\n', [patient_id, string(show_CEST_slice)]);
    %% show simulated CEST image for T1/T2/Flair
    fig1=figure('Name',patient_id+" CEST image");
    color_mask=imresize(cmp_brats(squeeze(Tumor_mask(:,:,CEST_slice(show_CEST_slice)))), [size(CEST_image_T1,1), size(CEST_image_T1,2)]);
    trans=ones([size(color_mask,1),size(color_mask,2)])*0.5; trans(sum(color_mask,3)==0)=0;
    subplot(3,7,1), imshow(color_mask), title("T1");
    subplot(3,7,8), imshow(color_mask), title("T2");
    subplot(3,7,15), imshow(color_mask), title("Flair");
    for f=1:6
        subplot(3,7,f+1),imshow(abs(squeeze(CEST_image_T1(:,:,(f-1)*9+1))),[0,1]),colormap(gca,"jet")
        title(num2str(scanned_offset_sorted((f-1)*9+1),'%.1f')+" ppm")
        subplot(3,7,f+8),imshow(abs(squeeze(CEST_image_T2(:,:,(f-1)*9+1))),[0,1]),colormap(gca,"jet")
        subplot(3,7,f+15),imshow(abs(squeeze(CEST_image_Flair(:,:,(f-1)*9+1))),[0,0.8]),colormap(gca,"jet")
    end
    set(fig1, 'position', screen_size); if close_flag; pause(pause_time), if has_last_fig; close(fig6); end; end
    %% show CEST image/ksp with chosen modality
    eval(['CEST_image = CEST_image_',char(show_modality),';']);
    ksp_CEST = fft2c_mri(squeeze(CEST_image));
    fig2=figure('Name',patient_id+" normalized CEST image");
    ratio=cal_ratio(CEST_image);
    for f=1:32
        slice = col(f); subplot(4,8,f),imshow(abs(squeeze(CEST_image(:,:,slice)))/ratio(slice),[0,0.8]),colormap("jet")
        title(num2str(scanned_offset_uniq(f),'%.1f')+" ppm")
    end
    set(fig2, 'position', screen_size); if close_flag; pause(pause_time), close(fig1); end
    fig3=figure('Name',patient_id+" normalized CEST kspace");
    for f=1:32
        slice = col(f); subplot(4,8,f),imshow(abs(squeeze(ksp_CEST(:,:,slice)))/ratio(slice),[0,0.5]),colormap("jet")
        title(num2str(scanned_offset_uniq(f),'%.1f')+" ppm")
    end
    set(fig3, 'position', screen_size); if close_flag; pause(pause_time), close(fig2); end
    %% show Z map & APTw img
    fig4=figure('Name',patient_id+" Z map");
    for f=1:32
        slice = col(f); subplot(4,8,f),imshow(abs(squeeze(CEST_image(:,:,slice)))./abs(CEST_image(:,:,1)),[0,1]),colormap("jet")
        title(num2str(scanned_offset_uniq(f),'%.1f')+" ppm")
    end
    set(fig4, 'position', screen_size); if close_flag; pause(pause_time), close(fig3); end
    fig5=figure('Name',patient_id+" APTw image");
    apt_neg=squeeze(CEST_image(:,:,47))./squeeze(CEST_image(:,:,1));
    apt_pos=squeeze(CEST_image(:,:,9))./squeeze(CEST_image(:,:,1));
    apt_roi=imclose(imbinarize(abs(squeeze(CEST_image(:,:,1))),0.2),strel('disk',20));
    rainbow_map=rainbow('rainbow_idl');
    subplot(1,3,1),imshow(abs(apt_pos),[0,1]),colormap("jet"),colorbar(),title("-3.5 ppm")
    subplot(1,3,2),imshow(abs(apt_neg),[0,1]),colormap("jet"),colorbar(),title("3.5 ppm")
    subplot(1,3,3),imshow((abs(apt_neg)-abs(apt_pos)).*apt_roi-(1-apt_roi),[-0.05,0.05]),colormap(gca,rainbow_map),colorbar(),title("APTw")
    set(fig5, 'position', screen_size); if close_flag; pause(pause_time), close(fig4); end
    %% show csm
    fig6=figure('Name',patient_id+" csm");
    for c=1:8
        subplot(2,4,c),imshow(abs(squeeze(CSM(:,:,c))),[]),title("coil "+string(c))
    end
    set(fig6, 'position', screen_size); has_last_fig=1; if close_flag; pause(pause_time), close(fig5); end
end