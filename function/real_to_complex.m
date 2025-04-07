%% create complex image with background & phase info
function [T1_cpx, T1ce_cpx, T2_cpx, Flair_cpx] = real_to_complex(T1_sk, T1ce_sk, T2_sk, Flair_sk, tumor_mask, fastmri_dir)
% background mask i.e. air & gap
bg_mask = zeros(size(tumor_mask));
bg_mask(tumor_mask==0) = 1;
bg_mask_phase = imerode(bg_mask, strel('disk',3));
gap_mask_phase = zeros(size(tumor_mask));
gap_mask_phase(tumor_mask==6) = 1;
gap_mask_phase = imerode(gap_mask_phase, strel('disk',3));
bg_mask_phase(gap_mask_phase==1) = 1;
bg_mask_amp = zeros(size(tumor_mask));
bg_mask_amp(tumor_mask==0 | tumor_mask==6) = 1;
bg_mask_amp_dilated = imdilate(bg_mask_amp, strel('disk',1));
% generate background amp & phase
[T1_bg_amp, T1_bg_phase] = create_background(T1_sk, bg_mask_amp_dilated, bg_mask_phase);
[T1ce_bg_amp, T1ce_bg_phase] = create_background(T1ce_sk, bg_mask_amp_dilated, bg_mask_phase);
[T2_bg_amp, T2_bg_phase] = create_background(T2_sk, bg_mask_amp_dilated, bg_mask_phase);
[Flair_bg_amp, Flair_bg_phase] = create_background(Flair_sk, bg_mask_amp_dilated, bg_mask_phase);
% add background to image
bg_mask_amp_edge = imfilter(edge_detection(bg_mask_amp_dilated, 2),fspecial('gaussian',3,3)','replicate');
T1_bg = T1_sk + T1_bg_amp.*bg_mask_amp; 
T1_bg = T1_bg.*(1-bg_mask_amp_edge) + imfilter(T1_sk,fspecial('gaussian',3,3)','replicate').*bg_mask_amp_edge;
T1ce_bg = T1ce_sk + T1ce_bg_amp.*bg_mask_amp;
T1ce_bg = T1ce_bg.*(1-bg_mask_amp_edge) + imfilter(T1ce_sk,fspecial('gaussian',3,3)','replicate').*bg_mask_amp_edge;
T2_bg = T2_sk + T2_bg_amp.*bg_mask_amp;
T2_bg = T2_bg.*(1-bg_mask_amp_edge) + imfilter(T2_sk,fspecial('gaussian',3,3)','replicate').*bg_mask_amp_edge;
Flair_bg = Flair_sk + Flair_bg_amp.*bg_mask_amp;
Flair_bg = Flair_bg.*(1-bg_mask_amp_edge) + imfilter(Flair_sk,fspecial('gaussian',3,3)','replicate').*bg_mask_amp_edge;
% figure()
% subplot(1,3,1),imshow(squeeze(T1ce_bg(:,:,1)),[0,1]), colormap('gray'), title("T1ce 0~1")
% subplot(1,3,2),imshow(squeeze(T1ce_bg(:,:,1)),[0,0.1]), colormap('gray'), title("T1ce 0~0.1")
% subplot(1,3,3),imshow(squeeze(T1ce_bg(:,:,1)),[0,0.01]), colormap('gray'), title("T1ce 0~0.01")
% generate phase for brain area
brain_mask_phase = ones(size(tumor_mask));
brain_mask_phase(bg_mask_phase==1) = 0;
[fastmri_folder_list,~] = get_sub_folder(fastmri_dir);
T1_br_phase = create_brain_phase_map(size(tumor_mask), fastmri_folder_list) .* brain_mask_phase;
T1ce_br_phase = create_brain_phase_map(size(tumor_mask), fastmri_folder_list) .* brain_mask_phase;
T2_br_phase = create_brain_phase_map(size(tumor_mask), fastmri_folder_list) .* brain_mask_phase;
Flair_br_phase = create_brain_phase_map(size(tumor_mask), fastmri_folder_list) .* brain_mask_phase;
% construct complex image
T1_cpx = r2c(T1_bg, imfilter(T1_bg_phase+T1_br_phase,fspecial('gaussian',3,3)','replicate'));
T1ce_cpx = r2c(T1ce_bg, imfilter(T1ce_bg_phase+T1ce_br_phase,fspecial('gaussian',3,3)','replicate'));
T2_cpx = r2c(T2_bg, imfilter(T2_bg_phase+T2_br_phase,fspecial('gaussian',3,3)','replicate'));
Flair_cpx = r2c(Flair_bg, imfilter(Flair_bg_phase+Flair_br_phase,fspecial('gaussian',3,3)','replicate'));
end
%% transfer real to complex according to cosine & sine
function img_cpx=r2c(img,phase)
img_real=img.*phase(:,:,:,1);
img_imag=img.*phase(:,:,:,2);
img_cpx=img_real+1j*img_imag;
end