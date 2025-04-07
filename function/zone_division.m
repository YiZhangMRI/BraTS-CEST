%% divide image into different parts according to tumor_mask
function Division = zone_division(Tumor_mask)
% 0:bg||1:necrotic||2:edema||3:no_enhance||4:enhanced_tumor||5:normal||6:skull||7:skin
Normal = zeros(size(Tumor_mask));
Normal(Tumor_mask==5) = 1;
Tumor = zeros(size(Tumor_mask));
Tumor(Tumor_mask==3 | Tumor_mask==4) = 1;
Edema = zeros(size(Tumor_mask));
Edema(Tumor_mask==2) = 1;
Necro = zeros(size(Tumor_mask));
Necro(Tumor_mask==1) = 1;
Skull = zeros(size(Tumor_mask));
Skull(Tumor_mask==6) = 1;
Skin = zeros(size(Tumor_mask));
Skin(Tumor_mask==7) = 1;
Background = zeros(size(Tumor_mask));
Background(Tumor_mask==0) = 1;
Division = cat(4, Normal, Tumor, Edema, Necro, Skull, Skin, Background);
end