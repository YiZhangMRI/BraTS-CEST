%% colormapping to Brats document
function color_mask = cmp_brats(tumor_mask)
% set 0:black||1:green||2:yellow||3:red||4:blue||5:magenta||6:cyan||7:orange
% 0:bg||1:necrotic||2:edema||3:no_enhance||4:enhanced_tumor||5:normal||6:skull||7:skin
color_mask_r = zeros(size(tumor_mask));
color_mask_r(tumor_mask==2 | tumor_mask==3 | tumor_mask==5 | tumor_mask==7)=1;
color_mask_g = zeros(size(tumor_mask));
color_mask_g(tumor_mask==1 | tumor_mask==2 | tumor_mask==6)=1;
color_mask_g(tumor_mask==7)=0.5;
color_mask_b = zeros(size(tumor_mask));
color_mask_b(tumor_mask==4 | tumor_mask==5 | tumor_mask==6)=1;
color_mask = cat(3,color_mask_r,color_mask_g,color_mask_b);
end