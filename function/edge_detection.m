%% detecte edge of area
function edge = edge_detection(area, width)
area(area>0)=1;
edge=area-imerode(area,strel('disk',width));
end