%% return a gaussian_window centralized according to img content
function window = gaussian_window_centralized(img, alpha_row, alpha_column)
window=zeros(size(img));
for f=1:size(img,3)
    frame=squeeze(img(:,:,f));
    brain_region=zeros(size(frame));
    brain_region(abs(frame)>0.05*max(frame(:)))=1;
    brain_region = imclose(brain_region,strel('disk', 20));
    [corner_row_0, corner_column_0]=find_corner_point(brain_region, 0);
    [corner_row_45, corner_column_45]=find_corner_point(brain_region, 45);
    [corner_row_90, corner_column_90]=find_corner_point(brain_region, 90);
    corner_row=cat(1,corner_row_0,corner_row_45,corner_row_90);
    corner_column=cat(1,corner_column_0,corner_column_45,corner_column_90);
    corner_row_u=unique(corner_row);
    corner_column_u=unique(corner_column);
    center_row=0.5*(max(corner_row_u(1),2)+min(corner_row_u(length(corner_row_u)),size(frame,1)-2));
    center_column=0.5*(max(corner_column_u(1),2)+min(corner_column_u(length(corner_column_u)),size(frame,2)-2));
    gaussain_window=gausswin(size(frame,1),alpha_row)*gausswin(size(frame,2),alpha_column)';
    window(:,:,f)=imtranslate(gaussain_window,[center_column-size(frame,2)/2, center_row-size(frame,1)/2],'FillValues',0,'OutputView','same');
end
end