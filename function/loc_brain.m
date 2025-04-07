%% locate whole brain region & record to mask as label 5
function extend_tumor_mask = loc_brain(T2, tumor_mask)
extend_tumor_mask = uint8(zeros(size(T2)));
extend_tumor_mask(T2>0 & tumor_mask==0) = 5;
extend_tumor_mask = extend_tumor_mask + uint8(tumor_mask);
end