%% add noise accrording to MTR
function img_noisy = ksp_add_noise(img, snr)
    ratio = cal_ratio(img);
    ratio = reshape(ratio, [1,1,length(ratio)]);
    noise_ratio = (1-ratio).^(4.*ratio);
%     figure()
%     subplot(2,1,1),plot(squeeze(ratio));
%     subplot(2,1,2),plot(squeeze(noise_ratio));
    high_freq_mask = gausswin(size(img,1), 0.5)*gausswin(size(img,2), 4)';
    noise = (randn(size(img))+randn(size(img))*1j)*snr;
    noise = noise .*high_freq_mask;
    noise(:, fix(normalize(randn([1,round(size(img,2)/5)]),'range')*(size(img,2)-1))+1) = 0;
    noise(:, round(size(img,2)*0.45):round(size(img,2)*0.55))=0;
%     figure()
%     subplot(1,4,1),imshow(high_freq_mask,[0,1])
%     subplot(1,4,2),imshow(abs(ifft2c_mri(high_freq_mask)),[0,1]),colormap("jet")
%     subplot(1,4,3),imshow(abs(squeeze(noise(:,:,1))),[0,snr])
%     subplot(1,4,4),imshow(abs(ifft2c_mri(squeeze(noise(:,:,1)))),[0,snr])
    ksp = fft2c_mri(img./ratio);
    ksp_noisy = ksp + noise.*noise_ratio;
    img_noisy = ifft2c_mri(ksp_noisy.*ratio);
%     figure()
%     subplot(2,2,1),imshow(abs(squeeze(img(:,:,1)./ratio(:,:,1))),[0,0.8])
%     subplot(2,2,2),imshow(abs(squeeze(img_noisy(:,:,1)./ratio(:,:,1))),[0,0.8])
%     subplot(2,2,3),imshow(abs(squeeze(img(:,:,27)./ratio(:,:,27))),[0,0.8])
%     subplot(2,2,4),imshow(abs(squeeze(img_noisy(:,:,27)./ratio(:,:,27))),[0,0.8]),colormap("jet")
end