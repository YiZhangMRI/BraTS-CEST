%% assign Zspetrum based on T1/T2 discret map
function Z_image = zsp_map(T1_map, T2_map, roi, t_levels, tissue)
Z_image = zeros([size(T1_map), 54]);
% create concentration map for amide/gag pool & MT
Z_file = load("./Data/simuZ/z_"+tissue+"_T1_1_T2_1.mat");
eval(['Z_file=Z_file.z_',char(tissue),';']);
ag_level = length(fieldnames(Z_file));
mt_level = size(Z_file.z_1,1);
% fprintf(['out ',num2str(ag_level),' ',num2str(mt_level),'\n']);
ag_cmap = conc_map(size(T1_map), roi, ag_level, 0.1, 0.9);
mt_cmap = conc_map(size(T1_map), roi, mt_level, 0.1, 0.9);
for l_T1=1:t_levels
   for l_T2=1:t_levels
       Z_file=load("./Data/simuZ/z_"+tissue+"_T1_"+string(l_T1)+"_T2_"+string(l_T2)+".mat");
       eval(['Z_file=Z_file.z_',char(tissue),';']);
       loc = (T1_map==l_T1 & T2_map==l_T2);
       for row=1:size(T1_map, 1)
        for col=1:size(T1_map, 2)
         if loc(row,col)==1
           ag=char(num2str(ag_cmap(row,col)));
           mt=char(num2str(mt_cmap(row,col)));
           zone=char(num2str(fix(9*rand(1)+1)));
%            fprintf(['in ',ag,' ',mt,' ',zone,'\n']);
           eval(['Zspetrum = squeeze(Z_file.z_', ag,'(', mt,',', zone,',:));']);
           % add random distortion at Zspetrum by 1%
           Zspetrum = Zspetrum .* (1+(normalize(randn(size(Zspetrum)),'range')-0.5)*0.02);
           Zspetrum(1) = 1;
           Z_image(row,col,:)=Zspetrum;
         end    
        end
       end
   end
end
end