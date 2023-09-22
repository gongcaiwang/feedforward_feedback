attrs = {'spec','spec_','loudness','loudness_','f0_hz','f0_hz_','amplitudes','amplitudes_','amplitude_formants_hamon','amplitude_formants_hamon_','freq_formants_hamon','freq_formants_hamon_','amplitude_formants_noise','amplitude_formants_noise_','freq_formants_noise','freq_formants_noise_','bandwidth_formants_noise_hz','bandwidth_formants_noise_hz_','freq_harm_diff','freq_harm_diff_'};

% onregion = {'cSTG','mSTG','parstriangularis','parsopercularis','precentral','postcentral'};
% mask = isregion(regions,onregion);
for a = 1:length(attrs)
    data = eval(attrs{a});
    data = data(:,1:8:end,:,:,:,:);
    eval([attrs{a},'=','data;']);
end