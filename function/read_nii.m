%% read T1\T2\FLAIR\tumor_mask .nii file in src & normalized with max
function [T1, T1ce, T2, Flair, Tumor_mask] = read_nii(src)
addpath ./function/Tools_for_NIfTI_and_ANALYZE_image
loc=dir(fullfile(src,'*t1.nii'));
T1 = load_nii(fullfile(src,loc.name));
T1 = double(rot90(T1.img));
T1 = T1./max(T1(:));
% T1 = T1./(0.2*max(T1(:))); % for patient 193
% T1(T1>1) = 1;
loc=dir(fullfile(src,'*t1ce.nii'));
T1ce = load_nii(fullfile(src,loc.name));
T1ce = double(rot90(T1ce.img));
T1ce = T1ce./max(T1ce(:));
% T1ce = T1ce./(0.2*max(T1ce(:))); % for patient 193
% T1ce(T1ce>1) = 1;
loc=dir(fullfile(src,'*t2.nii'));
T2 = load_nii(fullfile(src,loc.name));
T2 = double(rot90(T2.img));
T2 = T2./max(T2(:));
loc=dir(fullfile(src,'*flair.nii'));
Flair = load_nii(fullfile(src,loc.name));
Flair = double(rot90(Flair.img));
Flair = Flair./max(Flair(:));
% Flair = Flair./(0.2*max(Flair(:))); % for patient 193
% Flair(Flair>1) = 1;
loc=dir(fullfile(src,'*seg.nii'));
Tumor_mask = load_nii(fullfile(src,loc.name));
Tumor_mask = rot90(Tumor_mask.img);
end