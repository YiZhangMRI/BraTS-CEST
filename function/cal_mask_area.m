%% calculate each label's area to find the most likely tumor slice
function area = cal_mask_area(tumor_mask)
lesion=zeros([1,size(tumor_mask,3)]);
tumor_core=zeros([1,size(tumor_mask,3)]);
enhance=zeros([1,size(tumor_mask,3)]);
necrotic=zeros([1,size(tumor_mask,3)]);
brain=zeros([1,size(tumor_mask,3)]);
for s=1:size(tumor_mask,3)
    tumor_slice=squeeze(tumor_mask(:,:,s));
    necrotic(s)=length(tumor_slice(tumor_slice==1));
    enhance(s)=length(tumor_slice(tumor_slice==4));
    tumor_core(s)=length(tumor_slice(tumor_slice==3))+necrotic(s)+enhance(s);
    lesion(s)=length(tumor_slice(tumor_slice==2))+tumor_core(s);
    brain(s)=length(tumor_slice(tumor_slice==5))+lesion(s);
end
area = cat(1,lesion,tumor_core,enhance,necrotic,brain);
end