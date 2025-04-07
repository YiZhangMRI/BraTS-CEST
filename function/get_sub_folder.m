%% read all subfolder in src
function [s, len] = get_sub_folder(src)
    folder = dir(src);
    subfolder = folder(3:end); % first 2 are "." & ".."
    s = [];
    len = length(subfolder);
    for i = 1:len
        path = fullfile(src, subfolder(i).name);
        if isfolder(path)
            [subpath, sublen] = get_sub_folder(path); % recursion
            if sublen==0
                s = [s, string(path)]; % append s with self
            else
                s = [s, subpath]; % append s with new folder
            end
        elseif isfile(path)
            len = len-1;
        end
    end
end