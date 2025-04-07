%% get random nature texture
function texture = gen_nature_texture(texture_size, nature_dir)
dirs = dir([nature_dir,'*.jpg']);
filename = [nature_dir, dirs(randi([1,length(dirs)])).name];
nature_img = double(rgb2gray(imread(filename)))/255;
nature_img = abs(imresize(nature_img, texture_size, 'nearest'));
texture = imfilter(nature_img, fspecial('gaussian',15,15)', 'replicate');
texture = 1+0.5*(normalize(texture,'range')-0.5); % 0.75~1.25
end