%% create background with random air noise & slight aliasing pattern
function [bg_amp, bg_phase] = create_background(img, bg_mask_amp, bg_mask_phase)
img_size = size(img);
% subsample mask to create artifacts
sub_sample = ones(img_size);
alpha = rand(1)*20+10;
for f=1:img_size(3)
    sub_sample(:, randi([1,img_size(2)], [1,round(img_size(2)*(1/alpha))]), f) = 0;
end
sub_sample(:, round(img_size(2)*0.45):round(img_size(2)*0.55), :) = 1;
alias_ksp = fft2c_mri(img);
% randomly overwrite with ore ksp to simulate moving artifacts
parfor f=1:img_size(3)
    sub_sample_slice = squeeze(sub_sample(:,:,f));
    sub_sample_line = sum(1-sub_sample_slice(1,:));
    alias_ksp_slice = squeeze(alias_ksp(:,:,f));
    alias_ksp_slice(sub_sample_slice==0) = alias_ksp_slice(:, randi([1,img_size(1)*0.4], [1,sub_sample_line]));
    alias_ksp(:,:,f) = alias_ksp_slice;
end
alias_img=ifft2c_mri(alias_ksp);
% generate phase map with artifact & white noise
phase_alias = angle(alias_img) .* gaussian_window_centralized(img, 2.5, 2);
phase_alias = phase_alias .* bg_mask_phase;
phase_alias = (phase_alias/(max(phase_alias(:))*0.5))*pi;
phase_alias = (phase_alias + (rand(size(alias_img))-0.5)*2*pi.*(1-gaussian_window_centralized(img, 3, 1.5))) .* bg_mask_phase;
% phase_alias = imfilter(phase_alias, fspecial('gaussian',3,3)', 'replicate');
phase_alias(phase_alias>pi) = pi;
phase_alias(phase_alias<-pi) = -pi;
% generate texture for artifact & noise
nature_dir = '.\Data\Natural_image\';
nature_texture = gen_nature_texture([img_size(1),img_size(2)], nature_dir);
alias_img = abs(alias_img) .* nature_texture;
bg_mask_amp_edge = edge_detection(bg_mask_amp, 2);
% figure()
% subplot(1,4,1),imshow(squeeze(sub_sample(:,:,1)),[0,1]),title("sub sample mask")
% subplot(1,4,1),imshow(squeeze(abs(alias(:,:,1))),[0,0.1]),title("alias amp")
% amplitude normalization
alias_edge = alias_img / (max(alias_img(bg_mask_amp_edge==1))*20);
alias_body = alias_img .* (1-bg_mask_amp_edge)/max(alias_img(bg_mask_amp_edge==0));
alias_edge = imfilter(alias_edge + alias_body, fspecial('gaussian', 3, 3)', 'replicate') .* bg_mask_amp_edge;
alias_img = alias_body.*bg_mask_amp + alias_edge;
alias_img = alias_img / max(alias_img(:));
alias_window = gaussian_window_centralized(img, 1, 0.5);
% return phase in cosine & sine mode for simplify
bg_amp = (rand(img_size).*bg_mask_amp.*(1-alias_window)*0.02 + alias_img.*alias_window*0.5)*(0.5+rand(1)*1);
bg_phase = cat(4, cos(phase_alias).*bg_mask_phase, sin(phase_alias).*bg_mask_phase);
% subplot(1,4,2),imshow(squeeze(abs(alias(:,:,1))),[0,1]),title("alias amp")
% subplot(1,4,3),imshow(squeeze(phase_alias(:,:,1)),[-pi,pi]),title("bg phase")
% subplot(1,4,4),imshow(squeeze(abs(bg_amp(:,:,1))),[0,0.1]),title("bg amp")
end