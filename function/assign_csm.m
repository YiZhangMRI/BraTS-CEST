%% assign new CSM for CEST image using ESPIRiT
function [csm_new, reference] = assign_csm(img, csm, full_flag)
addpath ./function/ESPIRiT_code
img_coil = img.*csm;
coil_num = size(csm, 3);
calib_size = [12, 12]; % use 12 calibration lines to compute compression
ksize = [3,3]; % kernel size
% threshold of eigen vector decomposition in image space. 
eigThresh_1 = 0.02; if full_flag=="full"; eigThresh_2 = 0.25; else; eigThresh_2 = 0.95; end
% crop a calibration area
ksp_coil = fft2c(img_coil);
calib = crop(ksp_coil, [calib_size,coil_num]);
[k,S] = dat2Kernel(calib, ksize);
idx = find(S >= S(1)*eigThresh_1, 1, 'last');
if idx<coil_num; idx=coil_num; end
[M,W] = kernelEig(k(:,:,:,1:idx),[size(img_coil, 1), size(img_coil, 2)]);
P = sum(repmat(img_coil,[1,1,1,coil_num]).*conj(M),3);
reference = squeeze(P(:,:,1,16)); % reference img
csm_new = M(:,:,:,end).*repmat(W(:,:,end)>eigThresh_2,[1,1,coil_num]);
end