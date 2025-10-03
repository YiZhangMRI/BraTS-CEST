function zspectrum_fit = BM_simu_4pool_gag(input, poolConfig, numStep)
% full Bloch-McConell equation simulation 
% 2021.12.17 created by Yong Xingwang, 2023.06.05 modified by Wang Yuyan
% 4-pool: water, semisolid MT, gag-OH and gag-NH 
%% parse input
B1                = input.paramStruct.B1;
freq_selected     = input.paramStruct.freq_selected;  
gyr               = input.paramStruct.gyr;
water_pc          = input.paramStruct.water_pc;
isGaussWave       = input.paramStruct.isGaussWave;

%% pool params
T1w               = poolConfig.T1w;
T2w               = poolConfig.T2w;
T2m               = poolConfig.T2m;
Mmt               = poolConfig.Mmt/water_pc;
Mnh               = poolConfig.Mnh;
Moh               = poolConfig.Moh;

%% pool specification 
% water pool
M_0a = 1;     % pool a is water pool
R_2a = 1/T2w;               
R_1a = 1/T1w;

% gag-NH, reference, DOI: 10.1007/s10334-020-00868-y
CESTfreq = 394;    % Hz, 394Hz==3.2ppm at 2.9T
M_0b = Mnh/water_pc;
R_2b = 1/0.01;
R_1b = 1/1;
k_ba = 50;
k_ab = k_ba * M_0b / M_0a;

% gag-OH, reference, DOI: 10.1007/s10334-020-00868-y
poolcfreq = 123.5;  % Hz, 123.5Hz==1ppm at 2.9T
M_0c = Moh/water_pc;
R_2c = 1/0.01;
R_1c = 1/1;
k_ca = 1000;
k_ac = k_ca * M_0c / M_0a;
    
% MT, i.e. semisolid pool
R_2s = 1/T2m;
R_1s = 1/1; 
M_0s = Mmt;
k_sa = 60;
k_as = k_sa * M_0s / M_0a;

M0   = [0; 0; M_0a; 0; 0; M_0b; 0; 0; M_0c; 0; 0; M_0s; 1]; % a=water, b=APT, c=GAG/NOE, s=MT

%% pulse specification
[omega_1, crush_counter_init, crush_counter_max, step_size, numPulse,tp,td] = get_RF_pulse_mod(input, numStep);

%% simulation
off_resonance_range = freq_selected*2*pi;    % Hz to rad/s
M_z_w = zeros(size(off_resonance_range));
parfor iter = 1:numel(off_resonance_range)
    off_resonance = off_resonance_range(iter);
    delta_omega_a = off_resonance;
    delta_omega_b = off_resonance - CESTfreq*2*pi;
    delta_omega_c = off_resonance - poolcfreq*2*pi;
    delta_omega_s = off_resonance;
    
    M = M0;
    if isGaussWave
        crush_counter = crush_counter_init;
        for time_iter = 1: length(omega_1)
            A = [-R_2a-k_ab-k_ac-k_as,    -delta_omega_a   ,          0          ,     k_ba     ,          0         ,         0         ,     k_ca     ,          0         ,         0         ,     k_sa     ,          0         ,         0         ,     0    ;
                    delta_omega_a    , -R_2a-k_ab-k_ac-k_as,  omega_1(time_iter) ,       0      ,        k_ba        ,         0         ,       0      ,        k_ca        ,         0         ,       0      ,        k_sa        ,         0         ,     0    ;
                          0          , -omega_1(time_iter) , -R_1a-k_ab-k_ac-k_as,       0      ,          0         ,        k_ba       ,       0      ,          0         ,        k_ca       ,       0      ,          0         ,        k_sa       , R_1a*M_0a;
                         k_ab        ,          0          ,          0          ,  -R_2b-k_ba  ,   -delta_omega_b   ,         0         ,       0      ,          0         ,         0         ,       0      ,          0         ,         0         ,     0    ;
                          0          ,         k_ab        ,          0          , delta_omega_b,     -R_2b-k_ba     , omega_1(time_iter),       0      ,          0         ,         0         ,       0      ,          0         ,         0         ,     0    ;
                          0          ,          0          ,         k_ab        ,       0      , -omega_1(time_iter),     -R_1b-k_ba    ,       0      ,          0         ,         0         ,       0      ,          0         ,         0         , R_1b*M_0b;
                         k_ac        ,          0          ,          0          ,       0      ,          0         ,         0         ,  -R_2c-k_ca  ,   -delta_omega_c   ,         0         ,       0      ,          0         ,         0         ,     0    ;
                          0          ,         k_ac        ,          0          ,       0      ,          0         ,         0         , delta_omega_c,     -R_2c-k_ca     , omega_1(time_iter),       0      ,          0         ,         0         ,     0    ;
                          0          ,          0          ,         k_ac        ,       0      ,          0         ,         0         ,       0      , -omega_1(time_iter),     -R_1c-k_ca    ,       0      ,          0         ,         0         , R_1c*M_0c;
                         k_as        ,          0          ,          0          ,       0      ,          0         ,         0         ,       0      ,          0         ,         0         ,  -R_2s-k_sa  ,   -delta_omega_s   ,         0         ,     0    ;
                          0          ,         k_as        ,          0          ,       0      ,          0         ,         0         ,       0      ,          0         ,         0         , delta_omega_s,     -R_2s-k_sa     , omega_1(time_iter),     0    ;
                          0          ,          0          ,         k_as        ,       0      ,          0         ,         0         ,       0      ,          0         ,         0         ,       0      , -omega_1(time_iter),     -R_1s-k_sa    , R_1s*M_0s;
                          0          ,          0          ,          0          ,       0      ,          0         ,         0         ,       0      ,          0         ,         0         ,       0      ,          0         ,         0         ,     0    ;];

            if omega_1(time_iter) == 0
                crush_counter = crush_counter + 1;
            end
            if crush_counter == 1%RF_delay/step_size%omega_1(time_iter) == 0 % crushing
                M = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                     0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                     0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                     0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;]*M;
            elseif crush_counter == crush_counter_max
%                 crush_counter = 0;
                crush_counter = crush_counter_init;
            end
            M = expm(A*step_size)*(M);
        end
    else
        A = [-R_2a-k_ab-k_ac-k_as,    -delta_omega_a   ,          0          ,     k_ba     ,       0       ,      0     ,     k_ca     ,       0       ,      0     ,     k_sa     ,       0       ,      0     ,     0    ;
                delta_omega_a    , -R_2a-k_ab-k_ac-k_as,     gyr*B1*2*pi     ,       0      ,      k_ba     ,      0     ,       0      ,      k_ca     ,      0     ,       0      ,      k_sa     ,      0     ,     0    ;
                      0          ,     -gyr*B1*2*pi    , -R_1a-k_ab-k_ac-k_as,       0      ,       0       ,    k_ba    ,       0      ,       0       ,    k_ca    ,       0      ,       0       ,    k_sa    , R_1a*M_0a;
                     k_ab        ,          0          ,          0          ,  -R_2b-k_ba  , -delta_omega_b,      0     ,       0      ,       0       ,      0     ,       0      ,       0       ,      0     ,     0    ;
                      0          ,         k_ab        ,          0          , delta_omega_b,   -R_2b-k_ba  , gyr*B1*2*pi,       0      ,       0       ,      0     ,       0      ,       0       ,      0     ,     0    ;
                      0          ,          0          ,         k_ab        ,       0      ,  -gyr*B1*2*pi , -R_1b-k_ba ,       0      ,       0       ,      0     ,       0      ,       0       ,      0     , R_1b*M_0b;
                     k_ac        ,          0          ,          0          ,       0      ,       0       ,      0     ,  -R_2c-k_ca  , -delta_omega_c,      0     ,       0      ,       0       ,      0     ,     0    ;
                      0          ,         k_ac        ,          0          ,       0      ,       0       ,      0     , delta_omega_c,   -R_2c-k_ca  , gyr*B1*2*pi,       0      ,       0       ,      0     ,     0    ;
                      0          ,          0          ,         k_ac        ,       0      ,       0       ,      0     ,       0      ,  -gyr*B1*2*pi , -R_1c-k_ca ,       0      ,       0       ,      0     , R_1c*M_0c;
                     k_as        ,          0          ,          0          ,       0      ,       0       ,      0     ,       0      ,       0       ,      0     ,  -R_2s-k_sa  , -delta_omega_s,      0     ,     0    ;
                      0          ,         k_as        ,          0          ,       0      ,       0       ,      0     ,       0      ,       0       ,      0     , delta_omega_s,   -R_2s-k_sa  , gyr*B1*2*pi,     0    ;
                      0          ,          0          ,         k_as        ,       0      ,       0       ,      0     ,       0      ,       0       ,      0     ,       0      ,  -gyr*B1*2*pi , -R_1s-k_sa , R_1s*M_0s;
                      0          ,          0          ,          0          ,       0      ,       0       ,      0     ,       0      ,       0       ,      0     ,       0      ,       0       ,      0     ,     0    ;];

       
        % during crusher, no B1
        A1 = [-R_2a-k_ab-k_ac-k_as,    -delta_omega_a   ,          0          ,     k_ba     ,       0       ,     0     ,     k_ca     ,       0       ,     0     ,     k_sa     ,       0       ,     0     ,     0    ;
                 delta_omega_a    , -R_2a-k_ab-k_ac-k_as,          0          ,       0      ,      k_ba     ,     0     ,       0      ,      k_ca     ,     0     ,       0      ,      k_sa     ,     0     ,     0    ;
                       0          ,          0          , -R_1a-k_ab-k_ac-k_as,       0      ,       0       ,    k_ba   ,       0      ,       0       ,    k_ca   ,       0      ,       0       ,    k_sa   , R_1a*M_0a;
                      k_ab        ,          0          ,          0          ,  -R_2b-k_ba  , -delta_omega_b,     0     ,       0      ,       0       ,     0     ,       0      ,       0       ,     0     ,     0    ;
                       0          ,         k_ab        ,          0          , delta_omega_b,   -R_2b-k_ba  ,     0     ,       0      ,       0       ,     0     ,       0      ,       0       ,     0     ,     0    ;
                       0          ,          0          ,         k_ab        ,       0      ,       0       , -R_1b-k_ba,       0      ,       0       ,     0     ,       0      ,       0       ,     0     , R_1b*M_0b;
                      k_ac        ,          0          ,          0          ,       0      ,       0       ,     0     ,  -R_2c-k_ca  , -delta_omega_c,     0     ,       0      ,       0       ,     0     ,     0    ;
                       0          ,         k_ac        ,          0          ,       0      ,       0       ,     0     , delta_omega_c,   -R_2c-k_ca  ,     0     ,       0      ,       0       ,     0     ,     0    ;
                       0          ,          0          ,         k_ac        ,       0      ,       0       ,     0     ,       0      ,       0       , -R_1c-k_ca,       0      ,       0       ,     0     , R_1c*M_0c;
                      k_as        ,          0          ,          0          ,       0      ,       0       ,     0     ,       0      ,       0       ,     0     ,  -R_2s-k_sa  , -delta_omega_s,     0     ,     0    ;
                       0          ,         k_as        ,          0          ,       0      ,       0       ,     0     ,       0      ,       0       ,     0     , delta_omega_s,   -R_2s-k_sa  ,     0     ,     0    ;
                       0          ,          0          ,         k_as        ,       0      ,       0       ,     0     ,       0      ,       0       ,     0     ,       0      ,       0       , -R_1s-k_sa, R_1s*M_0s;
                       0          ,          0          ,          0          ,       0      ,       0       ,     0     ,       0      ,       0       ,     0     ,       0      ,       0       ,     0     ,     0    ;];

          
        for k = 1:numPulse
            M = expm(A*tp)*M;
            M = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;]*M;
            M = expm(A1*ceil(td/step_size)*step_size)*M; % relaxation
        end
    end
    M_z_w(iter) = M(3);
end

zspectrum_fit =  M_z_w;
end