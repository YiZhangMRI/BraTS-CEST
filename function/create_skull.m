%% create dummy skin & skull according to corner point & convex spline curve
function [T1, T1ce, T2, Flair, tumor_mask] = create_skull(T1, T1ce, T2, Flair, tumor_mask)
addpath ./function/interpclosed
% find brain region
brain_region = zeros(size(T2));
brain_region(tumor_mask~=0) = 1;
% edge smoothing & update brain area 5 in tumor_mask
brain_edge = imdilate(brain_region, strel('disk', 1)) - imerode(brain_region, strel('disk', 1));
tumor_mask(brain_edge==1 & tumor_mask==0) = 5;
brain_edge = imfilter(brain_edge, fspecial('gaussian',3,3)', 'replicate');
T1 = T1.*(1-brain_edge) + imfilter(T1, fspecial('gaussian',3,3)', 'replicate') .* brain_edge;
T1ce = T1ce.*(1-brain_edge) + imfilter(T1ce, fspecial('gaussian',3,3)', 'replicate') .* brain_edge;
T2 = T2.*(1-brain_edge) + imfilter(T2, fspecial('gaussian',3,3)', 'replicate') .* brain_edge;
Flair = Flair.*(1-brain_edge) + imfilter(Flair, fspecial('gaussian',3,3)', 'replicate') .* brain_edge;
% randomly select skin & skull thickness
org_skin_thick = randi([2,4]);
org_skull_think = randi([3,5]);
% remove sharp details to form a rounded brain edge for corner point detection 
brain_region(T2>0) = 1;
brain_region = imclose(brain_region,strel('disk', 40));
brain_region = imdilate(brain_region,strel('disk', org_skin_thick+org_skull_think));
% find corner point at each angle to locate brain bouding box
skin_region = zeros(size(tumor_mask));
skull_region = zeros(size(tumor_mask));
angle=-40:10:40; % from 4 o'clock to 2 o'clock
parfor f=1:size(tumor_mask,3)
    slice = squeeze(brain_region(:,:,f));
    corner_row = zeros([1,4*length(angle)]);
    corner_column = zeros([1,4*length(angle)]);
    for a=1:length(angle)
        [corner_row(4*a-3:4*a), corner_column(4*a-3:4*a)] = find_corner_point(slice, angle(a));
    end
    % sort corner point counterclockwise
    corner_point = zeros(2,length(corner_row));
    for a=1:length(angle)
        corner_point(:,(1:4)*length(angle)-(length(angle)-a))=[corner_row((a-1)*4+1:a*4);corner_column((a-1)*4+1:a*4)];
    end
    % sort left side to have decaying rows <¡ü>
    left_side = cat(2, corner_point(:, 4*length(angle)-1:4*length(angle)), corner_point(:, 1:length(angle)+2));
    left_side = sortrows(left_side',1,'descend')';
    corner_point(:, 4*length(angle)-1:4*length(angle)) = left_side(:, 1:2);
    corner_point(:, 1:length(angle)+2) = left_side(:, 3:length(angle)+4);
    % sort upper side to have increasing columns <¡ú>
    upper_side = corner_point(:, length(angle)-1:2*length(angle)+2);
    upper_side = sortrows(upper_side',2,'ascend')';
    corner_point(:, length(angle)-1:2*length(angle)+2) = upper_side;
    % sort right side to have increasing rows <¡ý>
    right_side = corner_point(:, 2*length(angle)-1:3*length(angle)+2);
    right_side = sortrows(right_side',1,'ascend')';
    corner_point(:, 2*length(angle)-1:3*length(angle)+2) = right_side;
    % sort bottom side to have decaying columns <¡û>
    bottom_side = cat(2,corner_point(:, 3*length(angle)-1:4*length(angle)), corner_point(:, 1:2));
    bottom_side = sortrows(bottom_side',2,'descend')';
    corner_point(:, 3*length(angle)-1:4*length(angle)) = bottom_side(:,1:length(angle)+2);
    corner_point(:, 1:2) = bottom_side(:, length(angle)+3:length(angle)+4);
    % spline curve interpolation
    spline_xy = interpclosed(corner_point(1,:), corner_point(2,:), 0:0.001:1, 'spline');
    spline_xy = unique(round(spline_xy)', "row", "stable")';
%     figure()
%     imshow(slice,[0,1])
%     hold on, plot(spline_xy(2,:),spline_xy(1,:),'g-','linewidth',1.0)
%     hold on, plot(corner_point(2,1:length(angle)),corner_point(1,1:length(angle)),'ro')
%     hold on, plot(corner_point(2,1+length(angle):2*length(angle)),corner_point(1,1+length(angle):2*length(angle)),'bo')
%     hold on, plot(corner_point(2,1+2*length(angle):3*length(angle)),corner_point(1,1+2*length(angle):3*length(angle)),'mo')
%     hold on, plot(corner_point(2,1+3*length(angle):4*length(angle)),corner_point(1,1+3*length(angle):4*length(angle)),'co')
%     break
    % create skin with random thickness along spline, also include meninges
    skin_slice = zeros(size(brain_region,1),size(brain_region,2));
    thickness = org_skin_thick;
    for p=1:size(spline_xy,2)
        currentpoint = zeros(size(brain_region,1), size(brain_region,2));
        currentpoint(spline_xy(1,p), spline_xy(2,p))=1;
        % momentum update skin thickness
        if mod(p,10)==0
            thickness = round(0.75*thickness + 0.25*(org_skin_thick*(1.6*rand(1)+0.2)));
        end
        currentpoint = imdilate(currentpoint, strel('disk', thickness));
        skin_slice(currentpoint==1) = 1;
    end
    % find skull region between brain & skin, also include CSF & sinuses
    skin_slice = imclose(skin_slice, strel('disk',20));
    brain_plus_skin = zeros(size(skin_slice));
    brain_plus_skin(squeeze(tumor_mask(:,:,f))>0 | skin_slice>0) = 1;
    skull_region(:,:,f) = imclose(brain_plus_skin, strel('disk', 20)) - brain_plus_skin;
    skin_slice = imclose(brain_plus_skin, strel('disk', 40));
    skin_slice(squeeze(tumor_mask(:,:,f))>0 | squeeze(skull_region(:,:,f))>0) = 0;
    skin_region(:,:,f) = skin_slice;
end
% add nature texture to skull & gap
nature_dir = '.\Data\Natural_image\';
skull_T1_T2 = skull_region .* gen_nature_texture([size(T2,1),size(T2,2)], nature_dir) * 0.02*(1+(rand(1)-0.5)*0.2);
skull_T1_T2 = imfilter(skull_T1_T2, fspecial('gaussian',2,2)', 'replicate');
skull_Flair = skull_region .* gen_nature_texture([size(T2,1),size(T2,2)], nature_dir) * 0.01*(1+(rand(1)-0.5)*0.2);
skull_Flair = imfilter(skull_Flair, fspecial('gaussian',2,2)', 'replicate');
fat_T1 = skin_region .* gen_nature_texture([size(T2,1),size(T2,2)], nature_dir) * 0.75*(1+(rand(1)-0.5)*0.25);
fat_T1 = imfilter(fat_T1, fspecial('gaussian',3,3)', 'replicate');
fat_T2 = skin_region .* gen_nature_texture([size(T2,1),size(T2,2)], nature_dir) * 0.8*(1+(rand(1)-0.5)*0.25);
fat_T2 = imfilter(fat_T2, fspecial('gaussian',3,3)', 'replicate');
fat_Flair = skin_region .* gen_nature_texture([size(T2,1),size(T2,2)], nature_dir) * 0.5*(1+(rand(1)-0.5)*0.25);
fat_Flair = imfilter(fat_Flair, fspecial('gaussian',3,3)', 'replicate');
% update tumor mask
fat_mask = uint8(zeros(size(T2)));
% dilate 1 for gaussian smooth & mark skull as 7
skin_region = imdilate(skin_region, strel('disk', 1));
fat_mask(skin_region>0 & tumor_mask==0) = 7;
tumor_mask = fat_mask + tumor_mask;
% dilate 1 for gaussian smooth & mark gap as 6
skull_region = imdilate(skull_region,strel('disk', 1));
skull_mask = uint8(zeros(size(T2)));
skull_mask(skull_region>0 & tumor_mask==0) = 6;
tumor_mask = skull_mask + tumor_mask;
T1 = T1 + skull_T1_T2 + fat_T1;
T1ce = T1ce + skull_T1_T2 + fat_T1;
T2 = T2 + skull_T1_T2 + fat_T2;
Flair = Flair + skull_Flair + fat_Flair;
% re-normalized by max
T1 = T1./max(T1(:));
T1ce = T1ce./max(T1ce(:));
T2 = T2./max(T2(:));
Flair = Flair./max(Flair(:));
% show result
% show_slice = round(size(T2,3)*0.5);
% figure()
% subplot(2,3,1),imshow(squeeze(fat_region(:,:,show_slice-5)),[]);
% subplot(2,3,2),imshow(squeeze(fat_region(:,:,show_slice)),[]);
% subplot(2,3,3),imshow(squeeze(fat_region(:,:,show_slice+5)),[]);
% subplot(2,3,4),imshow(squeeze(skull_region(:,:,show_slice-5)),[]);
% subplot(2,3,5),imshow(squeeze(skull_region(:,:,show_slice)),[]);
% subplot(2,3,6),imshow(squeeze(skull_region(:,:,show_slice+5)),[]);
% figure()
% subplot(1,3,1),imshow(squeeze(brain_region(:,:,show_slice)),[0,1])
% subplot(1,3,2),imshow(squeeze(skull_region(:,:,show_slice)),[0,1])
% subplot(1,3,3),imshow(squeeze(fat_region(:,:,show_slice)),[0,1])
% figure()
% subplot(1,4,1),imshow(squeeze(skull_T1_T2(:,:,show_slice)),[0,0.1])
% subplot(1,4,2),imshow(squeeze(skull_Flair(:,:,show_slice)),[0,0.1])
% subplot(1,4,3),imshow(squeeze(fat_T2(:,:,show_slice)),[0,1])
% subplot(1,4,4),imshow(squeeze(fat_Flair(:,:,show_slice)),[0,1])
end