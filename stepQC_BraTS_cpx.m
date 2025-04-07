%% Quality control of BraTS_r2c, check simulated results
% 2023.05.20 by Wang Yuyan
addpath function
clc, clear, close all
%% preset params
cpx_folder = './Data/BraTS_complex';
[cpx_folder_list, patient_num] = get_sub_folder(cpx_folder);
fprintf('quality control for complex images in %s that has %s patients\n', [cpx_folder, string(patient_num)])
close_flag = 1; % whether to close previous fig
pause_time = 1; % pause time between figs
has_last_fig=0;
screen_size=get(0,'ScreenSize');
screen_size(1:2)=screen_size(1:2)+0.05*screen_size(3:4);
screen_size(3:4)=screen_size(3:4)*0.8;
%% read patient data
for patient = 1:patient_num
    %% load data
    patient_loc = cpx_folder_list(patient);
    patient_id = split(cpx_folder_list(patient),filesep);
    patient_id = patient_id(end);
    load(fullfile(patient_loc,patient_id+"_t1.mat"));
    load(fullfile(patient_loc,patient_id+"_t1ce.mat"));
    load(fullfile(patient_loc,patient_id+"_t2.mat")); 
    load(fullfile(patient_loc,patient_id+"_flair.mat"));
    load(fullfile(patient_loc,patient_id+"_seg.mat"));
    fprintf('{> patient %s with slice num %s\n', [patient_id, string(size(Tumor_mask,3))]);
    %% locate the best view slice
    Tumor_area = cal_mask_area(Tumor_mask);
    show_slice = find(Tumor_area(1,:)==max(Tumor_area(1,:)), 1);
    %% show amp image
    fig1=figure('Name', patient_id+" amp image");
    amp_range=[0,1];
    for i=1:3
        s = min(max(show_slice+50*(i-2),1),size(T2,3));
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
    set(fig1, 'position', screen_size); if close_flag; pause(pause_time), if has_last_fig; close(fig2); end; end
    %% show complex image with phase
    fig2=figure('Name',patient_id+" complex image");
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
    set(fig2, 'position', screen_size); has_last_fig=1; if close_flag; pause(pause_time), close(fig1); end
end
