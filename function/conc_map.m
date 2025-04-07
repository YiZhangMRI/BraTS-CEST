%% generate amide/gag/MT concentration map using natural texture
function C_map = conc_map(map_size, roi, levels, lb, ub)
nature_dir = '.\Data\Natural_image\';
nature_texture = gen_nature_texture(map_size, nature_dir);
config.levels = levels;
config.lb     = lb;
config.ub     = ub;
config.sort   = "ascend";
config.edge   = "noedge";
C_map         = amp_discretization(nature_texture, roi, config);
end