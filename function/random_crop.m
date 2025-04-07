%% randomly crop image around context center
function img_crop = random_crop(img, roi, shift_range)
img_crop = zeros(size(img));
% figure()
for f=1:size(img, 3)
    brain_region = roi(:,:,f);
    brain_region = imclose(brain_region, strel('disk',20));
    % locate brain region bounding box for each slice
    [corner_row_0, corner_column_0] = find_corner_point(brain_region, 0);
    [corner_row_30, corner_column_30] = find_corner_point(brain_region, 30);
    [corner_row_60, corner_column_60] = find_corner_point(brain_region, 60);
    corner_row = cat(1, corner_row_0, corner_row_30, corner_row_60);
    corner_column = cat(1, corner_column_0, corner_column_30, corner_column_60);
    corner_row_u = unique(corner_row);
    corner_column_u = unique(corner_column);
    % specify context box
    contex_edge = 10;
    context_left = corner_column_u(1)-contex_edge;
    context_right = corner_column_u(length(corner_column_u))+contex_edge;
    context_up = corner_row_u(1)-contex_edge;
    context_down = corner_row_u(length(corner_row_u))+contex_edge;
    context_width = max(context_right-context_left, context_down-context_up);
    context_column = round((context_left+context_right)/2);
    context_row = round((context_up+context_down)/2);
    context_left = max(round(context_column - context_width/2), 2);
    context_right = min(round(context_column + context_width/2), size(img,2)-1);
    context_up = max(round(context_row - context_width/2), 2);
    context_down = min(round(context_row + context_width/2), size(img,1)-1);
    % specify crop box
    iter = 1;
    while iter<10
        crop_left = randi([1, min(context_left, shift_range)]);
        crop_up = randi([1, min(context_up, shift_range)]);
        crop_width_a = max(context_right-crop_left, context_down-crop_up);
        crop_width_b = min(size(img,2)-crop_left, size(img,1)-crop_up);
        if crop_width_a < crop_width_b
            crop_width = randi([crop_width_a, crop_width_b]);
            break
        elseif iter==9
            crop_width = crop_width_b;
        end
        iter = iter+1;
    end
    frame_crop = img(crop_up:crop_up+crop_width, crop_left:crop_left+crop_width,f,:);
    img_crop(:,:,f,:) = imresize(frame_crop, [size(img,1), size(img,2)], 'nearest');
%     % plot boxs
%     context_region = zeros([size(img,1),size(img,2)]);
%     crop_region = zeros([size(img,1),size(img,2)]);
%     context_region(context_up:context_down, context_left:context_right) = 1;
%     crop_region(crop_up:crop_up+crop_width, crop_left:crop_left+crop_width) = 1;
%     if f<=12; subplot(2,6,f),imshow(brain_region+context_region+crop_region, [0,3]); end
end
end