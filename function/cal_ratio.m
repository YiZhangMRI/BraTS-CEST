%% estimate MTR of CEST img
function new_ratio=cal_ratio(img)
ksp=fft2c_mri(img);
mask=zeros(size(ksp));
mask(floor(size(img,1)/2)-5:floor(size(img,1)/2)+6,floor(size(img,2)/2)-5:floor(size(img,2)/2)+6,:)=1;
partial_ksp=ksp.*mask;
foggy_img=ifft2c_mri(partial_ksp);
center_img=squeeze(foggy_img(floor(size(img,1)/2)-19:floor(size(img,1)/2)+20,floor(size(img,2)/2)-19:floor(size(img,2)/2)+20,:));
center_img=reshape(center_img,[1600,size(img,3)]);
new_ratio=abs(mean(center_img,1));
new_ratio=new_ratio./new_ratio(1);
end