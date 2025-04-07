%% add motion aliasing accrording to MTR
function img_alias = ksp_add_alias(img, alpha)
    ratio = cal_ratio(img);
    ratio = reshape(ratio, [1,1,length(ratio)]);
    line_num = round(size(img,2)*alpha*(1-ratio).^(4.*ratio));
    ksp = fft2c_mri(img);
    ksp_base = zeros(size(ksp));
    ksp_motion = zeros(size(ksp));
    mask_motion = zeros(size(ksp));
    for f=1:size(ksp,3)
        pos=normalize(randn([1,line_num(1,1,f)]),'range');
        posnew=pos; posnew(pos>=0.5)=pos(pos>=0.5)-0.5; posnew(pos<0.5)=pos(pos<0.5)+0.5;
%         posnew=rand([1,line_num(1,1,f)]);
        pos=fix(posnew*(size(img,2)-1))+1;
        mask_motion(:, pos, f) = 1;
        mask_motion(:, round(size(img,2)*0.45):round(size(img,2)*0.55),f) = 0;
        line_pos=find(mask_motion(1,:,f)==1);
        ksp_base(:,:,f) = ksp(:,:,f).*(1-mask_motion(:,:,f));
        line_pos_moved=line_pos+round(rand([1,length(line_pos)]))*2-1;
        line_pos_moved=max(1, min(size(img,2), line_pos_moved));
        ksp_motion(:,line_pos,f) = ksp(:,line_pos_moved,f)*0.2 + ksp(:,line_pos,f)*0.8;
    end
%     figure()
%     subplot(1,4,1),imshow(abs(squeeze(ksp(:,:,27))),[0,0.1])
%     subplot(1,4,2),imshow(abs(squeeze(ksp_base(:,:,27))),[0,0.1])
%     subplot(1,4,3),imshow(abs(squeeze(ksp_motion(:,:,27))),[0,0.1])
    ksp=ksp_base+ksp_motion;
%     subplot(1,4,4),imshow(abs(squeeze(ksp(:,:,27))),[0,0.1]),colormap("jet")
    img_alias = ifft2c_mri(ksp);
end