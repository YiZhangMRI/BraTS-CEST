%% pseudo-steady state signal evolution of TSE
function T2decay_filter = TSE_pss(T2, mask_size)
tesp = 3.6/1000; % echo space time 3.6ms
ETL = mask_size(2);
echo_times = tesp.*(1:ETL);
decay_amp = exp(-echo_times./T2);
T2decay_filter = zeros(size(decay_amp));
T2decay_filter(ceil((mask_size(2)+1)/2):mask_size(2)) = decay_amp(1:2:mask_size(2));
T2decay_filter(ceil((mask_size(2)+1)/2)-1:-1:1) = decay_amp(2:2:mask_size(2));
T2decay_filter = repmat(T2decay_filter,[mask_size(1),1]);
end