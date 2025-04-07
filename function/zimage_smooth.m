%% smooth Zmap on the edge of regions
function Zmap_smoothed = zimage_smooth(Zmap, region, width)
edge = zeros(size(region));
for z=1:6
   edge(:,:,:,z) = edge_detection(squeeze(region(:,:,:,z)), width);
end
edge = sum(edge, 4);
edge(edge>0) = 1;
whole = sum(region(:,:,:,1:6), 4);
whole_edge = edge_detection(whole, width);
edge = edge - whole_edge;
edge = imfilter(edge, fspecial('gaussian', width*2, width*2)','replicate');
bg = squeeze(region(:,:,:,7));
Zmap_smoothed = imfilter(Zmap, fspecial('gaussian', 3, 3)','replicate');
Zmap_edge = imfilter(Zmap_smoothed, fspecial('gaussian', width*2, width*2)','replicate');
Zmap_smoothed = Zmap_smoothed.*(1-edge) + Zmap_edge.*edge;
Zmap_smoothed = Zmap_smoothed.*(whole-whole_edge) + Zmap.*(whole_edge+bg);
end