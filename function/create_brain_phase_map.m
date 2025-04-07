%% generate phase map for brain area from fastMRI data
function phase_map = create_brain_phase_map(map_size, fastmri_folder_list)
% randomly select a fastMRI data to copy phase map
image_dir = [];
while isempty(image_dir)
    data_dir = char(fastmri_folder_list(randi([1,length(fastmri_folder_list)])));
    image_dir = dir([data_dir,'\','espirit*.mat']);
end
data_id = randi([1,length(image_dir)]);
image_dir = [data_dir,'\',image_dir(data_id).name];
% load one APT_data and select S0 as phase img
fastmri_image = load(image_dir);
phase_image = fastmri_image.reference;
maxvalue = max(max(abs(phase_image)));
phase_image = phase_image./maxvalue;
% figure()
% subplot(2,2,1),imshow(abs(phase_image),[0,1]),colorbar(),title("amp")
% subplot(2,2,2),imshow(angle(phase_image),[-pi,pi]),colorbar(),title("angle")
% subplot(2,2,3),imshow(real(phase_image)./abs(phase_image),[-1,1]),colorbar(),title("cosine")
% subplot(2,2,4),imshow(imag(phase_image)./abs(phase_image),[-1,1]),colorbar(),title("sine")
% find region box for brain
max_abs = abs(phase_image);
max_abs = max(max_abs(:));
brain_region = imbinarize(abs(phase_image), 0.05*max_abs);
brain_region = imclose(brain_region, strel('disk',20));
brain_region = imdilate(brain_region, strel('disk',5));
[corner_row_0, corner_column_0] = find_corner_point(brain_region, 0);
[corner_row_45, corner_column_45] = find_corner_point(brain_region, 45);
[corner_row_90, corner_column_90] = find_corner_point(brain_region, 90);
corner_row=cat(1, corner_row_0, corner_row_45, corner_row_90);
corner_column=cat(1, corner_column_0, corner_column_45, corner_column_90);
% figure()
% subplot(1,2,1),imshow(abs(phase_image),[0,1]),title("fastmri")
% subplot(1,2,2),imshow(brain_region,[0,1]),title("phase region")
% hold on, plot(corner_column, corner_row, 'ro')
corner_row_u = unique(corner_row);
corner_column_u = unique(corner_column);
row_range = max(corner_row_u(1), 2) : min(corner_row_u(length(corner_row_u)), size(phase_image,1)-2);
column_range = max(corner_column_u(1), 2) : min(corner_column_u(length(corner_column_u)), size(phase_image,2)-2);
phase_image = phase_image(row_range, column_range);
phase_image = imresize(phase_image, [map_size(1),map_size(2)], 'nearest');
% smoothing the phase img and return phase map trough cosine & sine
phase_image = imfilter(phase_image, fspecial('gaussian',32,32)', 'replicate');
phase_cos = real(phase_image)./abs(phase_image);
phase_sin = imag(phase_image)./abs(phase_image);
% figure()
% subplot(1,4,1),imshow(abs(phase_image),[0,1]),title("amp");
% subplot(1,4,2),imshow(phase_cos,[-2,2]),title("cosine");
% subplot(1,4,3),imshow(phase_sin,[-2,2]),title("sine");
% subplot(1,4,4),imshow(angle(phase_image),[-pi,pi]),title("phase");
phase_cos = repmat(phase_cos, [1,1,map_size(3)]);
phase_sin = repmat(phase_sin, [1,1,map_size(3)]);
phase_map = cat(4, phase_cos, phase_sin);
end