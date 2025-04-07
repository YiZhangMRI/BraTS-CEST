%% T2 blur caused by TSE
function blur_image = T2_blur(image, T2decay_filter)
ksp = fft2c_mri(image);
blur_image = ifft2c_mri(ksp.*T2decay_filter);
end