% attrs = {'spec','loudness','f0_hz','amplitudes','amplitudes_','amplitude_formants_hamon','freq_formants_hamon','amplitude_formants_noise','freq_formants_noise','bandwidth_formants_noise_hz','freq_harm_diff','freq_harm_diff_'};
attrs = {'spec','spec_','loudness','loudness_','f0_hz','f0_hz_','amplitudes','amplitudes_','amplitude_formants_hamon','amplitude_formants_hamon_','freq_formants_hamon','freq_formants_hamon_','amplitude_formants_noise','amplitude_formants_noise_','freq_formants_noise','freq_formants_noise_','bandwidth_formants_noise_hz','bandwidth_formants_noise_hz_','freq_harm_diff','freq_harm_diff_'};

onregion = {'cSTG','mSTG','parstriangularis','parsopercularis','precentral','postcentral'};
% mask = isregion(regions,onregion);
mask = ones(15,15);
for a = 1:length(attrs)
    data = correctaud(eval(attrs{a}),mask,regions);
    eval([attrs{a},'=','data;']);
end

function [mask] = isregion(regions,target_region)
    sz = size(regions);
    mask = zeros(sz(1:length(sz)-1));
    for i =1:size(regions,1)
        for j=1:size(regions,2)
            flag = false;
            for s=1:length(target_region)
                if contains(reshape(regions(i,j,:),1,size(regions,3)),squeeze(target_region{s}))
                    flag=true;
                    break;
                end
            end
            if flag%contains(reshape(regions(i,j,:),1,size(regions,3)),'mSTG') || contains(reshape(regions(i,j,:),1,size(regions,3)),'cSTG')
                mask(i,j) = true;
            end
        end
    end
end

function [data] = correctaud(data,mask,regions)
%     onset = 48;
%     onset = 1;
%     onset = 128;
    onset = 6;
%     onset = 48;
    sz = size(regions);
    aud = zeros(sz(1:length(sz)-1));
    for i =1:size(regions,1)
        for j=1:size(regions,2)
            if contains(reshape(regions(i,j,:),1,size(regions,3)),'mSTG') || contains(reshape(regions(i,j,:),1,size(regions,3)),'cSTG') || contains(reshape(regions(i,j,:),1,size(regions,3)),'rSTG')
                aud(i,j) = true;
            end
        end
    end
    if length(size(data))==5
        aud_rep = repmat(reshape(aud,[1,1,1,size(aud,1),size(aud,2)]),size(data,1),onset,size(data,3),1,1);
        aud_rep = cat(2,aud_rep,zeros(size(data,1),size(data,2)-onset,size(data,3),size(data,4),size(data,5)));
    else
        aud_rep = repmat(reshape(aud,[1,1,1,size(aud,1),size(aud,2),1]),size(data,1),onset,size(data,3),1,1,size(data,6));
        aud_rep = cat(2,aud_rep,zeros(size(data,1),size(data,2)-onset,size(data,3),size(data,4),size(data,5),size(data,6)));
    end
    auds = find(aud_rep==1);
%     data(auds) = data(auds)*0.3;
    data(auds) = data(auds)*0.3;
%     if pre_max>0
%         pre = pre-pre_max;
%     end
end
