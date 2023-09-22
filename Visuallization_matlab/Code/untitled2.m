root = '/Users/ranwang/Documents/writen_paper/NER2020/717_prearti/';
load([root,'occ_value.mat']);
amplitudes_occ = amplitudes; freq_harm_diff_occ = freq_harm_diff;
load([root,'select_all_value.mat']);
amplitudes_ = repmat(reshape(amplitudes,[size(amplitudes,1),size(amplitudes,2),1,1,1]),[1,1,size(amplitudes_occ,3),size(amplitudes_occ,4),size(amplitudes_occ,5)]) - amplitudes_occ;
freq_harm_diff_ = repmat(reshape(freq_harm_diff,[size(freq_harm_diff,1),size(freq_harm_diff,2),1,1,1]),[1,1,size(freq_harm_diff_occ,3),size(freq_harm_diff_occ,4),size(freq_harm_diff_occ,5)]) - freq_harm_diff_occ;
load([root,'occ_cc_value.mat']);
save([root,'occ_cc_value.mat'],'spec','loudness','f0_hz','amplitudes','amplitude_formants_hamon','freq_formants_hamon','amplitude_formants_noise','freq_formants_noise','bandwidth_formants_noise_hz','freq_harm_diff','amplitudes_','freq_harm_diff_');