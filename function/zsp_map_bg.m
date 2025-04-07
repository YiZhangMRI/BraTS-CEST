%% assign Zspetrum for background
function Z_image_bg = zsp_map_bg(T1, T2, roi, Z_image_brain)
alias_region = gaussian_window_centralized(reshape(1-roi, [size(T1),1]), 1, 0.5);
alias_region = imbinarize(alias_region, 0.5);
alias_region(roi==0)=0;
roi_edge=edge_detection(roi, 3);
max_T1=max(T1(roi_edge==0 & alias_region==1));
max_T2=max(T2(roi_edge==0 & alias_region==1));
alias_mask = zeros(size(T1));
alias_mask((T1>0.1*max_T1 | T2>0.1*max_T2) & alias_region==1)=1;
alias_mask=imclose(alias_mask,strel('disk',1));
bg_mask = 1-alias_mask;
bg_mask(roi~=1)=0;
% figure()
% subplot(1,3,1),imshow(alias_region,[0,1]);
% subplot(1,3,2),imshow(alias_mask,[0,1]);
% subplot(1,3,3),imshow(bg_mask,[0,1]);
T1_alias = T1; T1_alias(alias_mask~=1) = 0;
T2_alias = T2; T2_alias(alias_mask~=1) = 0;
T1_bg = T1; T1_bg(alias_mask==1) = 0;
T2_bg = T2; T2_bg(alias_mask==1) = 0;
alias_config.levels = 4;
alias_config.lb     = 0.05;
alias_config.ub     = 0.5;
alias_config.sort   = "descend";
alias_config.edge   = "edge";
alias_mask_dilated  = imdilate(alias_mask, strel('disk',10));
T1_alias_dismap     = amp_discretization(T1_alias, alias_mask_dilated, alias_config);
alias_config.sort   = "ascend";
T2_alias_dismap     = amp_discretization(T2_alias, alias_mask_dilated, alias_config);
bg_config.levels = 4;
bg_config.lb     = 0.05;
bg_config.ub     = 0.5;
bg_config.sort   = "descend";
bg_config.edge   = "noedge";
T1_bg_dismap     = amp_discretization(T1_bg, bg_mask, bg_config);
bg_config.sort   = "ascend";
T2_bg_dismap     = amp_discretization(T2_bg, bg_mask, bg_config);
Z_image_bg = zsp_map(T1_bg_dismap, T2_bg_dismap, bg_mask, 4, "bg");
% add random distortion about air noise
random_amp = (1+(2*rand(size(Z_image_bg))-1)*0.1).*(round(rand(size(Z_image_bg)))*2-1);
% figure()
% subplot(1,5,1),imshow(squeeze(Z_image_bg(:,:,2)),[0,1]),colormap("jet")
Z_image_bg = Z_image_bg .* random_amp;
% subplot(1,5,2),imshow(squeeze(random_amp(:,:,2)),[-1,1]),colormap("jet")
% subplot(1,5,3),imshow(squeeze(Z_image_bg(:,:,2)),[-1,1]),colormap("jet")
Z_image_alias = zsp_map(T1_alias_dismap, T2_alias_dismap, alias_mask_dilated, 4, "alias");
% subplot(1,5,4),imshow(squeeze(Z_image_alias(:,:,2)),[0,1]),colormap("jet")
Z_image_alias = Z_image_alias.*roi+Z_image_brain;
Z_image_alias = imfilter(Z_image_alias, fspecial('gaussian', 3, 3)','replicate').*alias_mask;
alias_random_mask = edge_detection(1-roi+alias_mask, 3);
Z_alias_random = imfilter(Z_image_alias+Z_image_brain+ones(size(Z_image_alias)).*(roi-alias_mask), fspecial('gaussian', 3, 3)','replicate');
Z_alias_random = Z_alias_random .* alias_random_mask .* random_amp;
Z_image_bg = (Z_image_bg + Z_image_alias.*(alias_mask-alias_random_mask) + Z_alias_random) .* roi;
% subplot(1,5,5),imshow(abs(squeeze(Z_image_bg(:,:,2))),[0,1]),colormap("jet")
end