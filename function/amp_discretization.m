%% discretize the input into several levels by amplitude
function discret_map = amp_discretization(image, roi, config)
levels = config.levels;
lb     = config.lb;
ub     = config.ub;
sort   = config.sort;
edge   = config.edge;
discret_map=zeros(size(image));
% not conting edge for max & min
if edge=="edge"
    edge_roi = edge_detection(roi,3);
    amp_max=max(image(roi~=0 & edge_roi~=1));
    amp_min=min(image(roi~=0 & edge_roi~=1));
else
    amp_max=max(image(roi~=0));
    amp_min=min(image(roi~=0));
end
% linear discretization in lb ~ ub range
amp_softmax=amp_max - (1-ub)*(amp_max-amp_min);
amp_softmin=amp_min + lb*(amp_max-amp_min);
amp_dis=amp_softmin:(amp_softmax-amp_softmin)/levels:amp_softmax;
amp_max=max(image(roi~=0)); % recalculate max & min after you got amp_dis
amp_min=min(image(roi~=0));
if ~isempty(amp_dis)
    if sort=="ascend"
        for a=1:levels
            discret_map(roi~=0 & image>=amp_dis(a) & image<=amp_dis(a+1)) = a;
        end
        discret_map(roi~=0 & image>=amp_min & image<=amp_dis(1)) = 1;
        discret_map(roi~=0 & image>=amp_dis(levels+1) & image<=amp_max) = levels;
    elseif sort=="descend"
        for a=1:levels
            discret_map(roi~=0 & image>=amp_dis(a) & image<=amp_dis(a+1)) = levels+1-a;
        end
        discret_map(roi~=0 & image>=amp_min & image<=amp_dis(1)) = levels;
        discret_map(roi~=0 & image>=amp_dis(levels+1) & image<=amp_max) = 1;
    end
end
end