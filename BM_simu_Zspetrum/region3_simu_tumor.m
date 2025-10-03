%% simulation of tumor tissue using B-M equation
% 2023.06.05 by Wang Yuyan
clear, close all
%% input params
save_dir = fullfile('..','Data','simuZ');
if ~exist(save_dir,'dir'); mkdir(save_dir); end
scanned_offset = transpose(load(fullfile('BM_simu_Zspetrum','Offset frequencies for APT scans - 54pts.txt')));
% scanned_offset = transpose(load(fullfile('BM_simu_Zspetrum','Offset frequencies for APT scans - 32pts.txt')));
scanned_offset(1) = [];  % get rid of infinite value
scanned_offset_sorted = sort(scanned_offset,'descend');
basicConfig.paramStruct.B1 = 2; % 2¦ÌT
basicConfig.paramStruct.freq_selected = scanned_offset_sorted;  
basicConfig.paramStruct.gyr = 42.576; % Larmor freq
basicConfig.paramStruct.water_pc = 74.9; % water pool proton concentration
basicConfig.paramStruct.isGaussWave = true; % RF pulse shape, gauss=1, rect=0
basicConfig.paramStruct.isNOE = false; % include NOE or not
numStep = 64;  % divide gauss RF-pulse into parts
z_num = 10; % number of Z-spetrums at each param set
num_all = 0;
%% pool params
T1w_range=[1.2, 3.0];
T2w_range=[0.1, 0.5];
Mapt_range=[0.2, 0.48];
Mmt_range=[2.5, 5.5];
poolConfig.T2m = 4e-5;
poolConfig.Mnoe = 0.5;
%% discretization of params range
gray_scale=-1./log(0.1:0.1:0.9);
gray_scale=gray_scale/max(gray_scale);
T1w_range_d=[T1w_range(1), T1w_range(1)+gray_scale.*(T1w_range(2)-T1w_range(1))];
T2w_range_d=[T2w_range(1), T2w_range(1)+gray_scale.*(T2w_range(2)-T2w_range(1))];
Mapt_range_d=Mapt_range(1)+(0:0.2:1).*(Mapt_range(2)-Mapt_range(1));
Mmt_range_d=Mmt_range(1)+(0:0.2:1).*(Mmt_range(2)-Mmt_range(1));
%% simulate z-spectra using BM equation
for n_T1w=1:length(T1w_range_d)-1
    T1w_low  = T1w_range_d(n_T1w);
    T1w_high = T1w_range_d(n_T1w+1);
    for n_T2w=1:length(T2w_range_d)-1
        T2w_low  = T2w_range_d(n_T2w);
        T2w_high = T2w_range_d(n_T2w+1);
        disp(['gray scale: T1 ',num2str(T1w_low),' T2 ',num2str(T2w_low)])
        for n_Mapt=1:length(Mapt_range_d)-1
            Mapt_low  = Mapt_range_d(n_Mapt);
            Mapt_high = Mapt_range_d(n_Mapt+1);
            disp(['amide grade ',num2str(n_Mapt)])
            zSpectrum = ones(length(Mmt_range_d)-1, z_num, length(scanned_offset_sorted)+1);
            for n_Mmt=1:length(Mmt_range_d)-1
                Mmt_low  = Mmt_range_d(n_Mmt);
                Mmt_high = Mmt_range_d(n_Mmt+1);
                for n_z=1:z_num
                    lb = [T1w_low,  T2w_low,  Mapt_low,  Mmt_low];
                    ub = [T1w_high, T2w_high, Mapt_high, Mmt_high];
                    pool_params = lb + (ub - lb).*rand(size(lb));
                    poolConfig.T1w  = pool_params(1);
                    poolConfig.T2w  = pool_params(2);
                    poolConfig.Mapt = pool_params(3);
                    poolConfig.Mmt  = pool_params(4);
                    zSpectrum_sim = BM_simu_4pool_amide(basicConfig, poolConfig, numStep);
                    zSpectrum(n_Mmt, n_z, 2:length(scanned_offset)+1) = zSpectrum_sim;
                end
            end
            eval(['z_tumor.z_',num2str(n_Mapt),'=','zSpectrum',';']);
        end
        save([save_dir,filesep,'z_tumor_T1_',num2str(n_T1w),'_T2_',num2str(n_T2w),'.mat'],'z_tumor');
        num_all = num_all+1;
        disp(['*** num of z-spectra : ',num2str(num_all*(length(Mapt_range_d)-1)*(length(Mmt_range_d)-1)*z_num)])
    end
end
%% show Z & MTRa
% Z=squeeze(z_tumor.z_1(1,1,2:length(scanned_offset_sorted)+1));
% freq=scanned_offset_sorted/123.47;
% min(Z)
% MTRa=Z(length(Z):-1:floor(length(Z)/2)+2)-Z(1:floor(length(Z)/2));
% figure(),plot(freq, Z, '-o'),set(gca,'XDir','reverse')
% ylim=get(gca,'Ylim'); hold on;plot([3.5,3.5],ylim,'m--');hold on;plot([-3.5,-3.5],ylim,'m--');
% figure(),plot(freq(1:floor(length(Z)/2)), MTRa, '-o','linewidth',1.2),set(gca,'XDir','reverse')
% ylim=get(gca,'Ylim');hold on;plot([3.5,3.5],ylim,'m--','linewidth',1.2)



