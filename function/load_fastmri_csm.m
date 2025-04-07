%% load CSM from fastMRI
function csm = load_fastmri_csm(csm_size, fastmri_folder_list)
csm_dir = [];
while isempty(csm_dir)
    patient_dir = char(fastmri_folder_list(randi([1,length(fastmri_folder_list)])));
    csm_dir = dir([patient_dir,'\','espirit*.mat']);
end
patient = fix(rand(1)*(length(csm_dir)-1))+1;
csm_dir = [patient_dir,'\',csm_dir(patient).name];
csm = load(csm_dir);
csm = imresize(csm.sensitivities, csm_size, 'nearest');
end