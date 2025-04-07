%% using bounding box to find the corner point of image content
function [corner_row, corner_column] = find_corner_point(slice, angle)
% rotation & transform matrix
slice_angle = imrotate(slice,angle,'crop');
rotM = [cosd(angle), sind(angle); -sind(angle), cosd(angle)];
% find corner at left & right
lr_search = find(slice_angle>0);
lr = [lr_search(1), lr_search(length(lr_search))];
lr_row = mod(lr, size(slice,1));
lr_column = [floor(lr(1)/size(slice,1)), ceil(lr(2)/size(slice,1))];
% find corner at up & down
ud_search = find(slice_angle'>0);
ud = [ud_search(1), ud_search(length(ud_search))];
ud_column = mod(ud,size(slice,2));
ud_row = [floor(ud(1)/size(slice,2)), ceil(ud(2)/size(slice,2))];
corner_row_angle = double([lr_row(1), ud_row(1), lr_row(2), ud_row(2)]);
corner_column_angle = double([lr_column(1), ud_column(1), lr_column(2), ud_column(2)]);
corner_angle_recover = round(rotM * [corner_row_angle-floor(size(slice,1)*0.5); corner_column_angle-floor(size(slice,2)*0.5)]);
corner_row = corner_angle_recover(1,:) + floor(size(slice,1)*0.5);
corner_column = corner_angle_recover(2,:) + floor(size(slice,2)*0.5);
end