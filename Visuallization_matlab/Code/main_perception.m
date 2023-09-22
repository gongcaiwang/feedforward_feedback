clc;
% set the root_dir to the directory where the Code and visualization data is located 
root_dir = '/Users/ranwang/Documents/writen_paper/NER2020/Visuallization_matlab';
% Choose a subject {'NY742', 'NY717', 'NY749'} -- we don't have subject brains other than these but MNI should work for all subjs
Subj = 'NY717';
% get paths to all visualization files
paths = get_path(root_dir,Subj);
% load electrode (channel) table
channel_info_all = get_channel_info(paths);
% get the relevant plotting options by choosing {'mni' or 'subj'} and hemisphere {'lh' or 'rh'} 
plot_data_all_mni = get_plot_data(paths, 'mni', 'lh');
% plot_data_all_T1 = get_plot_data(paths, 'subj', 'lh');

% %% Plotting electrodes on brain -- only one brain plot
% % select electrode range: to select grid elec only choose 1:128
% elec_range = 1:128;
% elec = table2array(channel_info_all(elec_range,plot_data_all_mni.coord_ind));
% elecname = table2array(channel_info_all(elec_range,11));
% % set radius of the elec circle
% radius = 1*ones(length(elec_range));
% % choose a color map and set indecies of colors for each electrode
% color_map = hot;
% color_ind = value2colorind(ones(128,1),'hot',[0,1]); 
% color_ind = [color_ind,color_ind,color_ind];
% colors = {color_ind, color_map};
% % plot on brain
% % fig = nyu_plot_whitebackground(plot_data_all.surf_brain,...
% %                          plot_data_all.sph,...
% %                          elec, elecname,...
% %                          colors, 0, plot_data_all.annot, radius, 1);
% fig = nyu_plot_whitebackground(plot_data_all_mni.surf_brain,...
%                          plot_data_all_mni.sph,...
%                          elec, elecname,...
%                          colors, 0, plot_data_all_mni.annot, radius, 1);
% % change view to front 
% view(fig.CurrentAxes,[-100,0]);

% %% Plotting attention map on brain -- one brain per time stamp
% elec_range = 1:128;
% elec = table2array(channel_info_all(elec_range,plot_data_all_T1.coord_ind));
% elecname = table2array(channel_info_all(elec_range,11));
% % set radius of the elec circle
% radius = 1*ones(length(elec_range));
% % choose a color map and set indecies of colors for each electrode
% color_map = hot(256);%afmhot;%hot;
% color_map = color_map(end:-1:1,:);
% % [color_map]=cbrewer('seq', 'Blues', 256);
% [color_map_div]=cbrewer('div', 'RdBu', 256);
% [colorset1] = cbrewer('qual', 'Set1', 9);
% color_map_div = color_map_div(end:-1:1,:);
% color_ind = value2colorind(ones(128,1),'hot',[0,1]); 
% color_ind = [color_ind,color_ind,color_ind];
% colors = {color_ind, color_map};
% 
% % Generating fake data 
% % [Att, mask, coord] = fake_data(elec);
% attrs = {'spec','loudness','f0_hz','amplitudes','freq_formants_hamon'};
% % attrs = {'amplitudes','freq_formants_noise','amplitude_formants_noise','bandwidth_formants_noise_hz'};
% figure();
% % haxis = subplot(10,3);
% 
% subjs = {'717','742','749'};
% SUB = length(subjs);
% plot_on_same_brain = true;
% above_median=true;
% col=size(attrs,2)+2;
% if plot_on_same_brain
%     row=1;
% else
%     row=SUB;
% end
% haxis = tight_subplot(row,col,[0,0],0,0);
% post = true;
% visdiff = false;
% postabs = true;%post || ~visdiff; % sparse if true
% max_att = true;%post || ~visdiff;
% max_anker = true;
% clrmap = color_map;
% ankercount = 5;
% if post
%     on = 16+32+1;
%     off = 116;
% else
%     on = 17;
%     off = 16+20;
%     on_post = 16+32+1;
%     off_post = 116;
% end
% % ind = [1:10,21:30];
% ind = [1:50];
% dumm = 1;
% for a=1:size(attrs,2)
%     maxminp1=[];
%     maxminp2=[10000];
%     for sub =1:SUB
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/NY',subjs{sub},'_elec.mat']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/att_NY',subjs{sub},'']);
%         att = eval(attrs{a});
%         if strcmp(attrs{a},'freq_formants_hamon')
%             maxp1 = gather_att(att(ind,on:off,1,:,:),postabs,max_att);
%             maxp2 = gather_att(att(ind,on:off,2,:,:),postabs,max_att);
%             maxminp1 = [maxminp1,max(maxp1(:))];
%             maxminp2 = [maxminp2,max(maxp2(:))];
%         else
%             maxp1 = gather_att(att(ind,on:off,1,:,:),postabs,max_att);
%             maxminp1 = [maxminp1,max(maxp1(:))];
%         end
%     end
%     maxminp1 = min(maxminp1);
%     maxminp2 = min(maxminp2);
%     Att1_cell = {};
%     Att2_cell = {};
%     diff_cell = {};
%     diff2_cell = {};
%     mask_cell = {};
%     coord_cell = {};
%     for sub=1:3
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/NY',subjs{sub},'_elec.mat']);
% %         load(['/Users/ranwang/Documents/writen_paper/NER2020/att_NY',subjs{sub},'']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/attr_dict_IG_value_',subjs{sub},'_nonnoise']);
%         Subj = ['NY',subjs{sub}];
%         % get paths to all visualization files
%         paths = get_path(root_dir,Subj);
%         channel_info_all = get_channel_info(paths);
%         plot_data_all_mni = get_plot_data(paths, 'mni', 'lh');
%         plot_data_all_T1 = get_plot_data(paths, 'subj', 'lh');
%     %     coord = T1;
%         coord = mni;
%         if strcmp(attrs{a},'freq_formants_hamon')
%             att = eval(attrs{a});
%             data = att(ind,on:off,1,:,:);
%             Att1 = gather_att(data,postabs,max_att);
%             sorted = sort(Att1(:));
%             % get all plot frames
%             if visdiff && ~post
% %                 Att1_post = gather_att(att(ind,on_post:off_post,1,:,:),postabs,max_att);
% %                 coef = find_normcoef(Att1_post,Att1,mask,max_anker,ankercount);
% %                 diff = Att1_post*coef(1)+coef(2) - Att1;
% %                 sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
% %                 diff = diff/maxv;
% %                 diff = (diff+1)/2;
% %                 [frames] = VisualAtt(diff, mask, coord, plot_data_all_mni,color_map_div,1);
%                 Att1_post = gather_att(att(ind,on_post:off_post,1,:,:),postabs,max_att);
%                 [Att1_post,m25,m75] = robust_rescale(Att1_post,mask);
%                 [Att1,m25,m75] = robust_rescale(Att1,mask);
%                 diff = Att1_post-Att1;
%                 sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                 diff = diff/maxv;
%                 diff = (diff+1)/2;
%                 diff_cell{sub} = diff;
%                 mask_cell{sub} = mask;
%                 coord_cell{sub} = coord;
%                 if plot_on_same_brain 
%                     if sub == SUB
%                         [frames] = VisualAtt(diff_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                         axes(haxis(sub2ind([col,row],dumm,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
%                     end
%                 else
%                     [frames] = VisualAtt({diff}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                     axes(haxis(sub2ind([col,row],dumm,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
%                 end
%             else
% %                 [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map,sorted(end-2));
%                 [Att1_md,m25,m75] = robust_rescale(Att1,mask);
%                 [Att1,m25,m75] = robust_rescale(Att1,mask,above_median);
%                 Att1_cell{sub} = Att1;
%                 mask_cell{sub} = mask;
%                 coord_cell{sub} = coord;
%                 if plot_on_same_brain 
%                     if sub == SUB
%                         [frames] = VisualAtt(Att1_cell, mask_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false);
%                         axes(haxis(sub2ind([col,row],dumm,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
%                     end
%                 else
%                     [frames] = VisualAtt({Att1}, {mask}, {coord}, plot_data_all_mni,color_map,1,false,false);
%                     axes(haxis(sub2ind([col,row],dumm,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
%                 end
%             end
%             
%             
%             data = att(ind,on:off,2,:,:);
%             Att2 = gather_att(data,postabs,max_att);
%             sorted = sort(Att2(:));
%             if visdiff && ~post
% %                 Att2_post = gather_att(att(ind,on_post:off_post,2,:,:),postabs,max_att);
% %                 coef = find_normcoef(Att2_post,Att2,mask,max_anker,ankercount);
% %                 diff = Att2_post*coef(1)+coef(2) - Att2;
% %                 sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
% %                 diff = diff/maxv;
% %                 diff = (diff+1)/2;
% %                 [frames] = VisualAtt(diff, mask, coord, plot_data_all_mni,color_map_div,1);
%                 Att2_post = gather_att(att(ind,on_post:off_post,1,:,:),postabs,max_att);
%                 [Att2_post,m25,m75] = robust_rescale(Att2_post,mask);
%                 [Att2,m25,m75] = robust_rescale(Att2,mask);
%                 diff = Att2_post-Att2;
%                 sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                 diff = diff/maxv;
%                 diff = (diff+1)/2;
%                 diff_cell{sub} = diff;
%                 mask_cell{sub} = mask;
%                 coord_cell{sub} = coord;
%                 if plot_on_same_brain 
%                     if sub == SUB
%                         [frames] = VisualAtt(diff_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                         axes(haxis(sub2ind([col,row],dumm+1,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
%                     end
%                 else
%                     [frames] = VisualAtt({diff}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                     axes(haxis(sub2ind([col,row],dumm+1,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
%                 end
%             else
% %                 [frames] = VisualAtt(Att2, mask, coord, plot_data_all_mni,color_map,sorted(end-2));
%                 [Att2_md,m25,m75] = robust_rescale(Att2,mask);
%                 [Att2,m25,m75] = robust_rescale(Att2,mask,above_median);
%                 Att2_cell{sub} = Att2;
%                 mask_cell{sub} = mask;
%                 coord_cell{sub} = coord;
%                 if plot_on_same_brain 
%                     if sub == SUB
%                         [frames] = VisualAtt(Att2_cell, mask_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false);
%                         axes(haxis(sub2ind([col,row],dumm+1,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
%                     end
%                 else
%                     [frames] = VisualAtt({Att2}, {mask}, {coord}, plot_data_all_mni,color_map,1,false,false);
%                     axes(haxis(sub2ind([col,row],dumm+1,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
%                 end
%             end
% 
%             
% %             coef = find_normcoef(Att2,Att1,mask,false,20);
% %             diff = Att2*coef(1)+coef(2) - Att1;
%             diff = Att2_md-Att1_md;
%             sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%             diff = diff/maxv;
%             diff = (diff+1)/2;
%             diff2_cell{sub} = diff;
%             mask_cell{sub} = mask;
%             coord_cell{sub} = coord;
%             if plot_on_same_brain 
%                 if sub == SUB
%                     [frames] = VisualAtt(diff2_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                     axes(haxis(sub2ind([col,row],dumm+2,1)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
%                 end
%             else
%                 [frames] = VisualAtt({diff}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                 axes(haxis(sub2ind([col,row],dumm+2,sub)));
%                 imshow(uint8(squeeze(frames)));
%                 axis off;
%             end
%             
%             
%         else
%             if strcmp(attrs{a},'freq_harm_diff')
%                 att = eval(attrs{a});
%                 data = att(ind,on:off,:,:,:);
%                 Att1 = mean(mean(squeeze(data),1),2);
%                 maxv = min(max(Att1(:)),abs(min(Att1(:))));
%                 Att1 = Att1/maxv;
%                 Att1 = (Att1+1)/2;
%                 % get all plot frames
% %                 [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map_div,1);
%                 [Att1,m25,m75] = robust_rescale(Att1,mask);
%                 Att1_cell{sub} = Att1;
%                 mask_cell{sub} = mask;
%                 coord_cell{sub} = coord;
%                 if plot_on_same_brain 
%                     if sub == SUB
%                         [frames] = VisualAtt(Att1_cell, mask_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false);
%                         axes(haxis(sub2ind([col,row],dumm,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
%                     end
%                 else
%                     [frames] = VisualAtt({Att1}, {mask}, {coord}, plot_data_all_mni,color_map,1,false,false);
%                     axes(haxis(sub2ind([col,row],dumm,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
%                 end
%             else
%                 att = eval(attrs{a});
%                 data = att(ind,on:off,1,:,:);
%                 sorted = sort(data(:));
%                 sorted = sorted(end-1);
%                 Att1 = gather_att(data,postabs,max_att);
%                 sorted = sort(Att1(:));
% %                 sorted = sorted(end-1);
%                 % get all plot frames
%                 if strcmp(attrs{a},'loudness') || strcmp(attrs{a},'f0_hz')
%                     if visdiff && ~post
% %                         Att1_post = gather_att(att(ind,on_post:off_post,1,:,:),postabs,max_att);
% %                         coef = find_normcoef(Att1_post,Att1,mask,max_anker,ankercount);
% %                         diff = Att1_post*coef(1)+coef(2) - Att1;
% %                         sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
% %                         diff = diff/maxv;
% %                         diff = (diff+1)/2;
% %                         [frames] = VisualAtt(diff, mask, coord, plot_data_all_mni,color_map_div,1);
%                         Att1_post = gather_att(att(ind,on_post:off_post,1,:,:),postabs,max_att);
%                         [Att1_post,m25,m75] = robust_rescale(Att1_post,mask);
%                         [Att1,m25,m75] = robust_rescale(Att1,mask);
%                         diff = Att1_post-Att1;
%                         sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                         diff = diff/maxv;
%                         diff = (diff+1)/2;
%                         diff_cell{sub} = diff;
%                         mask_cell{sub} = mask;
%                         coord_cell{sub} = coord;
%                         if plot_on_same_brain 
%                             if sub == SUB
%                                 [frames] = VisualAtt(diff_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                                 axes(haxis(sub2ind([col,row],dumm,1)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
%                             end
%                         else
%                             [frames] = VisualAtt({diff}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                             axes(haxis(sub2ind([col,row],dumm,sub)));
%                             imshow(uint8(squeeze(frames)));
%                             axis off;
%                         end
%                     else
% %                         [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map,maxminp1);
%                         [Att1,m25,m75] = robust_rescale(Att1,mask,above_median);
%                         Att1_cell{sub} = Att1;
%                         mask_cell{sub} = mask;
%                         coord_cell{sub} = coord;
%                         if plot_on_same_brain 
%                             if sub == SUB
%                                 [frames] = VisualAtt(Att1_cell, mask_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false);
%                                 axes(haxis(sub2ind([col,row],dumm,1)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
%                             end
%                         else
%                             [frames] = VisualAtt({Att1}, {mask}, {coord}, plot_data_all_mni,color_map,1,false,false);
%                             axes(haxis(sub2ind([col,row],dumm,sub)));
%                             imshow(uint8(squeeze(frames)));
%                             axis off;
%                         end
%                     end
%                 else
%                     if strcmp(attrs{a},'amplitudes')
%                         if visdiff && ~post
% %                             Att1_post = gather_att(att(ind,on_post:off_post,1,:,:),postabs,max_att);
% %                             coef = find_normcoef(Att1_post,Att1,mask,max_anker,ankercount);
% %                             diff = Att1_post*coef(1)+coef(2) - Att1;
% %                             sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
% %                             diff = diff/maxv;
% %                             diff = (diff+1)/2;
% %                             [frames] = VisualAtt(diff, mask, coord, plot_data_all_mni,color_map_div,1);
%                             Att1_post = gather_att(att(ind,on_post:off_post,1,:,:),postabs,max_att);
%                             [Att1_post,m25,m75] = robust_rescale(Att1_post,mask);
%                             [Att1,m25,m75] = robust_rescale(Att1,mask);
%                             diff = Att1_post-Att1;
%                             sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                             diff = diff/maxv;
%                             diff = (diff+1)/2;
%                             diff_cell{sub} = diff;
%                             mask_cell{sub} = mask;
%                             coord_cell{sub} = coord;
%                             if plot_on_same_brain 
%                                 if sub == SUB
%                                     [frames] = VisualAtt(diff_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                                     axes(haxis(sub2ind([col,row],dumm,1)));
%                                     imshow(uint8(squeeze(frames)));
%                                     axis off;
%                                 end
%                             else
%                                 [frames] = VisualAtt({diff}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                                 axes(haxis(sub2ind([col,row],dumm,sub)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
%                             end
%                         else
% %                             [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map,sorted(end-1));
%                             Att1 = mean(mean(squeeze(data),1),2);
%                             sorted = sort(Att1(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                             Att1 = Att1/maxv;
%                             Att1 = (Att1+1)/2;
%                             Att1_cell{sub} = Att1;
%                             mask_cell{sub} = mask;
%                             coord_cell{sub} = coord;
%                             if plot_on_same_brain
%                                 if sub == SUB
%                                     [frames] = VisualAtt(Att1_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                                     axes(haxis(sub2ind([col,row],dumm,1)));
%                                     imshow(uint8(squeeze(frames)));
%                                     axis off;
%                                 end
%                             else
%                                 [frames] = VisualAtt({Att1}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                                 axes(haxis(sub2ind([col,row],dumm,sub)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
%                             end
% 
%                         end
%                     else
%                         if visdiff && ~post
% %                             Att1_post = gather_att(att(ind,on_post:off_post,1,:,:),postabs,max_att);
% %                             coef = find_normcoef(Att1_post,Att1,mask,max_anker,ankercount);
% %                             diff = Att1_post*coef(1)+coef(2) - Att1;
% %                             sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
% %                             diff = diff/maxv;
% %                             diff = (diff+1)/2;
% %                             [frames] = VisualAtt(diff, mask, coord, plot_data_all_mni,color_map_div,1);
%                             Att1_post = gather_att(att(ind,on_post:off_post,1,:,:),postabs,max_att);
%                             [Att1_post,m25,m75] = robust_rescale(Att1_post,mask);
%                             [Att1,m25,m75] = robust_rescale(Att1,mask);
%                             diff = Att1_post-Att1;
%                             sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                             diff = diff/maxv;
%                             diff = (diff+1)/2;
%                             diff_cell{sub} = diff;
%                             mask_cell{sub} = mask;
%                             coord_cell{sub} = coord;
%                             if plot_on_same_brain
%                                 if sub == SUB
%                                     [frames] = VisualAtt(diff_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                                     axes(haxis(sub2ind([col,row],dumm,1)));
%                                     imshow(uint8(squeeze(frames)));
%                                     axis off;
%                                 end
%                             else
%                                 [frames] = VisualAtt({diff}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                                 axes(haxis(sub2ind([col,row],dumm,sub)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
%                             end
%                         else
% %                             [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map);
%                             [Att1,m25,m75] = robust_rescale(Att1,mask,above_median);
%                             Att1_cell{sub} = Att1;
%                             mask_cell{sub} = mask;
%                             coord_cell{sub} = coord;
%                             if plot_on_same_brain 
%                                 if sub == SUB
%                                     [frames] = VisualAtt(Att1_cell, mask_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false);
%                                     axes(haxis(sub2ind([col,row],dumm,1)));
%                                     imshow(uint8(squeeze(frames)));
%                                     axis off;
%                                 end
%                             else
%                                 [frames] = VisualAtt({Att1}, {mask}, {coord}, plot_data_all_mni,color_map,1,false,false);
%                                 axes(haxis(sub2ind([col,row],dumm,sub)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%         if ~plot_on_same_brain
%             if dumm == 2
%                 ylabel(Subj)
%             end  
%         end
%         if plot_on_same_brain
%             if sub == SUB
%                 title(attrs{a})
%             end
%         else
%             if sub ==1
%                 title(attrs{a})
%             end
%         end
%     end
%     if strcmp(attrs{a},'freq_formants_hamon')
%         dumm = dumm+3;
%     else
%         dumm = dumm+1;
%     end
%     
% end
% % plot all in one figure

% %% Directly Color the surface -- Not complete for attention map plotting
% surf_data.annot = [];
% surf_data.surf_brain = plot_data_all_T1.surf_brain;
% surf_data.sph = plot_data_all_T1.sph;
% elec_data.elec = elec;
% fig = plot_surface(surf_data, elec_data);

% 
% %% Plotting attention map on brain -- one brain per time stamp
% elec_range = 1:128;
% elec = table2array(channel_info_all(elec_range,plot_data_all_T1.coord_ind));
% elecname = table2array(channel_info_all(elec_range,11));
% % set radius of the elec circle
% radius = 1*ones(length(elec_range));
% % choose a color map and set indecies of colors for each electrode
% color_map = hot(256);%afmhot;%hot;
% color_map = color_map(end:-1:1,:);
% % [color_map]=cbrewer('seq', 'Blues', 256);
% [color_map_div]=cbrewer('div', 'RdBu', 256);
% [colorset1] = cbrewer('qual', 'Set1', 9);
% color_map_div = color_map_div(end:-1:1,:);
% color_ind = value2colorind(ones(128,1),'hot',[0,1]); 
% color_ind = [color_ind,color_ind,color_ind];
% colors = {color_ind, color_map};
% 
% % Generating fake data 
% % [Att, mask, coord] = fake_data(elec);
% attrs = {'spec','loudness','f0_hz','amplitudes','freq_formants_hamon'};
% % attrs = {'amplitudes','freq_formants_noise','amplitude_formants_noise','bandwidth_formants_noise_hz'};
% figure();
% % haxis = subplot(10,3);
% 
% subjs = {'717','742','749'};
% SUB = length(subjs);
% plot_on_same_brain = true;
% above_median=true;
% col=size(attrs,2)+2;
% if plot_on_same_brain
%     row=1;
% else
%     row=SUB;
% end
% haxis = tight_subplot(row,col,[0,0],0,0);
% post = true;
% visdiff = false;
% postabs = true;%post || ~visdiff; % sparse if true
% max_att = true;%post || ~visdiff;
% max_anker = true;
% clrmap = color_map;
% ankercount = 5;
% if post
%     on = 1;
%     off = 128;
% else
%     on = 17;
%     off = 16+20;
%     on_post = 16+32+1;
%     off_post = 116;
% end
% % ind = [1:10,21:30];
% ind = [1:50];
% dumm = 1;
% for a=1:size(attrs,2)
%     maxminp1=[];
%     maxminp2=[10000];
%     for sub =1:SUB
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/NY',subjs{sub},'_elec.mat']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/att_NY',subjs{sub},'']);
%         att = eval(attrs{a});
%         if strcmp(attrs{a},'freq_formants_hamon')
%             maxp1 = gather_att(att(ind,on:off,1,:,:),postabs,max_att);
%             maxp2 = gather_att(att(ind,on:off,2,:,:),postabs,max_att);
%             maxminp1 = [maxminp1,max(maxp1(:))];
%             maxminp2 = [maxminp2,max(maxp2(:))];
%         else
%             maxp1 = gather_att(att(ind,on:off,1,:,:),postabs,max_att);
%             maxminp1 = [maxminp1,max(maxp1(:))];
%         end
%     end
%     maxminp1 = min(maxminp1);
%     maxminp2 = min(maxminp2);
%     Att1_cell = {};
%     Att2_cell = {};
%     diff_cell = {};
%     diff2_cell = {};
%     mask_cell = {};
%     coord_cell = {};
%     for sub=1:3
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/NY',subjs{sub},'_elec.mat']);
% %         load(['/Users/ranwang/Documents/writen_paper/NER2020/att_NY',subjs{sub},'']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_dict_IG_value_',subjs{sub},'']);
%         Subj = ['NY',subjs{sub}];
%         % get paths to all visualization files
%         paths = get_path(root_dir,Subj);
%         channel_info_all = get_channel_info(paths);
%         plot_data_all_mni = get_plot_data(paths, 'mni', 'lh');
%         plot_data_all_T1 = get_plot_data(paths, 'subj', 'lh');
%     %     coord = T1;
%         coord = mni;
%         if strcmp(attrs{a},'freq_formants_hamon')
%             att = eval(attrs{a});
%             data = att(ind,on:off,1,:,:);
%             Att1 = gather_att2(data,postabs,max_att);
%             sorted = sort(Att1(:));
%             % get all plot frames
%             if visdiff && ~post
% %                 Att1_post = gather_att(att(ind,on_post:off_post,1,:,:),postabs,max_att);
% %                 coef = find_normcoef(Att1_post,Att1,mask,max_anker,ankercount);
% %                 diff = Att1_post*coef(1)+coef(2) - Att1;
% %                 sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
% %                 diff = diff/maxv;
% %                 diff = (diff+1)/2;
% %                 [frames] = VisualAtt(diff, mask, coord, plot_data_all_mni,color_map_div,1);
%                 Att1_post = gather_att2(att(ind,on_post:off_post,1,:,:),postabs,max_att);
%                 [Att1_post,m25,m75] = robust_rescale(Att1_post,mask);
%                 [Att1,m25,m75] = robust_rescale(Att1,mask);
%                 diff = Att1_post-Att1;
%                 sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                 diff = diff/maxv;
%                 diff = (diff+1)/2;
%                 diff_cell{sub} = diff;
%                 mask_cell{sub} = mask;
%                 coord_cell{sub} = coord;
%                 if plot_on_same_brain 
%                     if sub == SUB
%                         [frames] = VisualAtt(diff_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                         axes(haxis(sub2ind([col,row],dumm,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
%                     end
%                 else
%                     [frames] = VisualAtt({diff}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                     axes(haxis(sub2ind([col,row],dumm,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
%                 end
%             else
% %                 [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map,sorted(end-2));
%                 [Att1_md,m25,m75] = robust_rescale(Att1,mask);
%                 [Att1,m25,m75] = robust_rescale(Att1,mask,above_median);
%                 Att1_cell{sub} = Att1;
%                 mask_cell{sub} = mask;
%                 coord_cell{sub} = coord;
%                 if plot_on_same_brain 
%                     if sub == SUB
%                         [frames] = VisualAtt(Att1_cell, mask_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false);
%                         axes(haxis(sub2ind([col,row],dumm,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
%                     end
%                 else
%                     [frames] = VisualAtt({Att1}, {mask}, {coord}, plot_data_all_mni,color_map,1,false,false);
%                     axes(haxis(sub2ind([col,row],dumm,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
%                 end
%             end
%             
%             
%             data = att(ind,on:off,2,:,:);
%             Att2 = gather_att2(data,postabs,max_att);
%             sorted = sort(Att2(:));
%             if visdiff && ~post
% %                 Att2_post = gather_att(att(ind,on_post:off_post,2,:,:),postabs,max_att);
% %                 coef = find_normcoef(Att2_post,Att2,mask,max_anker,ankercount);
% %                 diff = Att2_post*coef(1)+coef(2) - Att2;
% %                 sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
% %                 diff = diff/maxv;
% %                 diff = (diff+1)/2;
% %                 [frames] = VisualAtt(diff, mask, coord, plot_data_all_mni,color_map_div,1);
%                 Att2_post = gather_att2(att(ind,on_post:off_post,1,:,:),postabs,max_att);
%                 [Att2_post,m25,m75] = robust_rescale(Att2_post,mask);
%                 [Att2,m25,m75] = robust_rescale(Att2,mask);
%                 diff = Att2_post-Att2;
%                 sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                 diff = diff/maxv;
%                 diff = (diff+1)/2;
%                 diff_cell{sub} = diff;
%                 mask_cell{sub} = mask;
%                 coord_cell{sub} = coord;
%                 if plot_on_same_brain 
%                     if sub == SUB
%                         [frames] = VisualAtt(diff_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                         axes(haxis(sub2ind([col,row],dumm+1,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
%                     end
%                 else
%                     [frames] = VisualAtt({diff}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                     axes(haxis(sub2ind([col,row],dumm+1,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
%                 end
%             else
% %                 [frames] = VisualAtt(Att2, mask, coord, plot_data_all_mni,color_map,sorted(end-2));
%                 [Att2_md,m25,m75] = robust_rescale(Att2,mask);
%                 [Att2,m25,m75] = robust_rescale(Att2,mask,above_median);
%                 Att2_cell{sub} = Att2;
%                 mask_cell{sub} = mask;
%                 coord_cell{sub} = coord;
%                 if plot_on_same_brain 
%                     if sub == SUB
%                         [frames] = VisualAtt(Att2_cell, mask_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false);
%                         axes(haxis(sub2ind([col,row],dumm+1,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
%                     end
%                 else
%                     [frames] = VisualAtt({Att2}, {mask}, {coord}, plot_data_all_mni,color_map,1,false,false);
%                     axes(haxis(sub2ind([col,row],dumm+1,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
%                 end
%             end
% 
%             
% %             coef = find_normcoef(Att2,Att1,mask,false,20);
% %             diff = Att2*coef(1)+coef(2) - Att1;
%             diff = Att2_md-Att1_md;
%             sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%             diff = diff/maxv;
%             diff = (diff+1)/2;
%             diff2_cell{sub} = diff;
%             mask_cell{sub} = mask;
%             coord_cell{sub} = coord;
%             if plot_on_same_brain 
%                 if sub == SUB
%                     [frames] = VisualAtt(diff2_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                     axes(haxis(sub2ind([col,row],dumm+2,1)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
%                 end
%             else
%                 [frames] = VisualAtt({diff}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                 axes(haxis(sub2ind([col,row],dumm+2,sub)));
%                 imshow(uint8(squeeze(frames)));
%                 axis off;
%             end
%             
%             
%         else
%             if strcmp(attrs{a},'freq_harm_diff')
%                 att = eval(attrs{a});
%                 data = att(ind,on:off,:,:,:);
%                 Att1 = mean(mean(squeeze(data),1),2);
%                 maxv = min(max(Att1(:)),abs(min(Att1(:))));
%                 Att1 = Att1/maxv;
%                 Att1 = (Att1+1)/2;
%                 % get all plot frames
% %                 [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map_div,1);
%                 [Att1,m25,m75] = robust_rescale(Att1,mask);
%                 Att1_cell{sub} = Att1;
%                 mask_cell{sub} = mask;
%                 coord_cell{sub} = coord;
%                 if plot_on_same_brain 
%                     if sub == SUB
%                         [frames] = VisualAtt(Att1_cell, mask_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false);
%                         axes(haxis(sub2ind([col,row],dumm,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
%                     end
%                 else
%                     [frames] = VisualAtt({Att1}, {mask}, {coord}, plot_data_all_mni,color_map,1,false,false);
%                     axes(haxis(sub2ind([col,row],dumm,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
%                 end
%             else
%                 att = eval(attrs{a});
%                 data = att(ind,on:off,1,:,:);
%                 sorted = sort(data(:));
%                 sorted = sorted(end-1);
%                 Att1 = gather_att2(data,postabs,max_att);
%                 sorted = sort(Att1(:));
% %                 sorted = sorted(end-1);
%                 % get all plot frames
%                 if strcmp(attrs{a},'loudness') || strcmp(attrs{a},'f0_hz')
%                     if visdiff && ~post
% %                         Att1_post = gather_att2(att(ind,on_post:off_post,1,:,:),postabs,max_att);
% %                         coef = find_normcoef(Att1_post,Att1,mask,max_anker,ankercount);
% %                         diff = Att1_post*coef(1)+coef(2) - Att1;
% %                         sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
% %                         diff = diff/maxv;
% %                         diff = (diff+1)/2;
% %                         [frames] = VisualAtt(diff, mask, coord, plot_data_all_mni,color_map_div,1);
%                         Att1_post = gather_att2(att(ind,on_post:off_post,1,:,:),postabs,max_att);
%                         [Att1_post,m25,m75] = robust_rescale(Att1_post,mask);
%                         [Att1,m25,m75] = robust_rescale(Att1,mask);
%                         diff = Att1_post-Att1;
%                         sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                         diff = diff/maxv;
%                         diff = (diff+1)/2;
%                         diff_cell{sub} = diff;
%                         mask_cell{sub} = mask;
%                         coord_cell{sub} = coord;
%                         if plot_on_same_brain 
%                             if sub == SUB
%                                 [frames] = VisualAtt(diff_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                                 axes(haxis(sub2ind([col,row],dumm,1)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
%                             end
%                         else
%                             [frames] = VisualAtt({diff}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                             axes(haxis(sub2ind([col,row],dumm,sub)));
%                             imshow(uint8(squeeze(frames)));
%                             axis off;
%                         end
%                     else
% %                         [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map,maxminp1);
%                         [Att1,m25,m75] = robust_rescale(Att1,mask,above_median);
%                         Att1_cell{sub} = Att1;
%                         mask_cell{sub} = mask;
%                         coord_cell{sub} = coord;
%                         if plot_on_same_brain 
%                             if sub == SUB
%                                 [frames] = VisualAtt(Att1_cell, mask_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false);
%                                 axes(haxis(sub2ind([col,row],dumm,1)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
%                             end
%                         else
%                             [frames] = VisualAtt({Att1}, {mask}, {coord}, plot_data_all_mni,color_map,1,false,false);
%                             axes(haxis(sub2ind([col,row],dumm,sub)));
%                             imshow(uint8(squeeze(frames)));
%                             axis off;
%                         end
%                     end
%                 else
%                     if strcmp(attrs{a},'amplitudes')
%                         if visdiff && ~post
% %                             Att1_post = gather_att2(att(ind,on_post:off_post,1,:,:),postabs,max_att);
% %                             coef = find_normcoef(Att1_post,Att1,mask,max_anker,ankercount);
% %                             diff = Att1_post*coef(1)+coef(2) - Att1;
% %                             sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
% %                             diff = diff/maxv;
% %                             diff = (diff+1)/2;
% %                             [frames] = VisualAtt(diff, mask, coord, plot_data_all_mni,color_map_div,1);
%                             Att1_post = gather_att2(att(ind,on_post:off_post,1,:,:),postabs,max_att);
%                             [Att1_post,m25,m75] = robust_rescale(Att1_post,mask);
%                             [Att1,m25,m75] = robust_rescale(Att1,mask);
%                             diff = Att1_post-Att1;
%                             sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                             diff = diff/maxv;
%                             diff = (diff+1)/2;
%                             diff_cell{sub} = diff;
%                             mask_cell{sub} = mask;
%                             coord_cell{sub} = coord;
%                             if plot_on_same_brain 
%                                 if sub == SUB
%                                     [frames] = VisualAtt(diff_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                                     axes(haxis(sub2ind([col,row],dumm,1)));
%                                     imshow(uint8(squeeze(frames)));
%                                     axis off;
%                                 end
%                             else
%                                 [frames] = VisualAtt({diff}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                                 axes(haxis(sub2ind([col,row],dumm,sub)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
%                             end
%                         else
% %                             [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map,sorted(end-1));
%                             Att1_1 = gather_att2(att(ind,on:off,2,:,:),postabs,max_att);
%                             Att1_0 = gather_att2(att(ind,on:off,1,:,:),postabs,max_att);
%                             [Att1_1,m25,m75] = robust_rescale(Att1_1,mask);
%                             [Att1_0,m25,m75] = robust_rescale(Att1_0,mask);
%                             diff = Att1_1-Att1_0;
%                             sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                             diff = diff/maxv;
%                             diff = (diff+1)/2;
%                             
%                             Att1 = mean(mean(squeeze(att(ind,on:off,2,:,:)),1),2);
%                             sorted = sort(Att1(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                             Att1 = Att1/maxv;
%                             Att1 = (Att1+1)/2;
% 
%                             Att1_cell{sub} = Att1;
%                             mask_cell{sub} = mask;
%                             coord_cell{sub} = coord;
%                             if plot_on_same_brain
%                                 if sub == SUB
%                                     [frames] = VisualAtt(Att1_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                                     axes(haxis(sub2ind([col,row],dumm,1)));
%                                     imshow(uint8(squeeze(frames)));
%                                     axis off;
%                                 end
%                             else
%                                 [frames] = VisualAtt({Att1}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                                 axes(haxis(sub2ind([col,row],dumm,sub)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
%                             end
% 
%                         end
%                     else
%                         if visdiff && ~post
% %                             Att1_post = gather_att2(att(ind,on_post:off_post,1,:,:),postabs,max_att);
% %                             coef = find_normcoef(Att1_post,Att1,mask,max_anker,ankercount);
% %                             diff = Att1_post*coef(1)+coef(2) - Att1;
% %                             sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
% %                             diff = diff/maxv;
% %                             diff = (diff+1)/2;
% %                             [frames] = VisualAtt(diff, mask, coord, plot_data_all_mni,color_map_div,1);
%                             Att1_post = gather_att2(att(ind,on_post:off_post,1,:,:),postabs,max_att);
%                             [Att1_post,m25,m75] = robust_rescale(Att1_post,mask);
%                             [Att1,m25,m75] = robust_rescale(Att1,mask);
%                             diff = Att1_post-Att1;
%                             sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                             diff = diff/maxv;
%                             diff = (diff+1)/2;
%                             diff_cell{sub} = diff;
%                             mask_cell{sub} = mask;
%                             coord_cell{sub} = coord;
%                             if plot_on_same_brain
%                                 if sub == SUB
%                                     [frames] = VisualAtt(diff_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                                     axes(haxis(sub2ind([col,row],dumm,1)));
%                                     imshow(uint8(squeeze(frames)));
%                                     axis off;
%                                 end
%                             else
%                                 [frames] = VisualAtt({diff}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                                 axes(haxis(sub2ind([col,row],dumm,sub)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
%                             end
%                         else
% %                             [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map);
%                             [Att1,m25,m75] = robust_rescale(Att1,mask,above_median);
%                             Att1_cell{sub} = Att1;
%                             mask_cell{sub} = mask;
%                             coord_cell{sub} = coord; 
%                             if plot_on_same_brain 
%                                 if sub == SUB
%                                     [frames] = VisualAtt(Att1_cell, mask_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false);
%                                     axes(haxis(sub2ind([col,row],dumm,1)));
%                                     imshow(uint8(squeeze(frames)));
%                                     axis off;
%                                 end
%                             else
%                                 [frames] = VisualAtt({Att1}, {mask}, {coord}, plot_data_all_mni,color_map,1,false,false);
%                                 axes(haxis(sub2ind([col,row],dumm,sub)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%         if ~plot_on_same_brain
%             if dumm == 2
%                 ylabel(Subj)
%             end  
%         end
%         if plot_on_same_brain
%             if sub == SUB
%                 title(attrs{a})
%             end
%         else
%             if sub ==1
%                 title(attrs{a})
%             end
%         end
%     end
%     if strcmp(attrs{a},'freq_formants_hamon')
%         dumm = dumm+3;
%     else
%         dumm = dumm+1;
%     end
%     
% end


%% 
elec_range = 1:128;
elec = table2array(channel_info_all(elec_range,plot_data_all_mni.coord_ind));
elecname = table2array(channel_info_all(elec_range,11));
% set radius of the elec circle
radius = 1*ones(length(elec_range));
% choose a color map and set indecies of colors for each electrode
color_map = hot(256);%afmhot;%hot;
color_map = color_map(end:-1:1,:);
% [color_map]=cbrewer('seq', 'Blues', 256);
[color_map_div]=cbrewer('div', 'RdBu', 256,'PCHIP');
[colorset1] = cbrewer('qual', 'Set1', 9,'PCHIP');
color_map_div = color_map_div(end:-1:1,:);
color_ind = value2colorind(ones(128,1),'hot',[0,1]); 
color_ind = [color_ind,color_ind,color_ind];
colors = {color_ind, color_map};
global noiselevel;
noiselevel = 0.19;%0.19;
% Generating fake data 
% [Att, mask, coord] = fake_data(elec);
% attrs = {'loudness'};
% attrs = {'freq_formants_hamon'};
% attrs = {'spec','loudness','f0_hz','amplitudes','freq_formants_hamon','freq_formants_noise','bandwidth_formants_noise_hz'};
% attrs = {'spec','loudness','f0_hz','freq_formants_hamon'};
% attrs = {'spec','freq_formants_hamon','freq_formants_noise','bandwidth_formants_noise_hz','freq_harm_diff'};
attrs = {'spec'};
% attrs = {'amplitudes'};
% % causal
% ylims_min = [0,0,0,0,0,0,0];
% ylims_max = [1,1,1,1,1,1.5,1];
% anticausal
% ylims_min = [0,0,0,0.3,0,0,0];
% ylims_max = [1,2,3,3,2,6,1];
ylims_min = [0,0,0,0,0,0,0];
% ylims_max = [1,1.5,1.5,1.5,1,2.0,1];
ylims_max = [1,1,1,1,1,1.5,1];
% % percept
% ylims_min = [0,0,0,0.3,0,0,0];
% ylims_max = [1,4,3,6,2,6,1.5];
% attrs = {'amplitudes'};
% attrs = {'amplitudes','freq_formants_noise','amplitude_formants_noise','bandwidth_formants_noise_hz'};
% haxis = subplot(10,3);

wreg=false;
temp=false;
tau=false;
return_curve=false;
entire_period = true;
post = true;
visdiff = false;
postabs = true;%post || ~visdiff; % sparse if true
max_att = false;%post || ~visdiff;
max_anker = true;
clrmap = color_map;

subjs = {'717','742','749'};
SUB = length(subjs);
plot_on_same_brain = true;
above_median=false;
[found,ind]=ismember('freq_formants_hamon',attrs);

if found
    row=size(attrs,2)+2;
else
    row=size(attrs,2)+1;
end
N_tau = 16;

if temp
    col = 13;
else
    if plot_on_same_brain
        col=1;
    else
        col=SUB;
    end
end

invest_regions = {{'cSTG','rSTG','mMTG','cMTG','rMTG',...
    'inferiorprecentral','superiorprecentral','postcentral','supramarginal',...
    'parsopercularis','parstriangularis','parsorbitalis','rostralmiddlefrontal','caudalmiddlefrontal'}};

all_regions = {'cSTG','rSTG','mMTG','cMTG','rMTG',...
    'inferiorprecentral','superiorprecentral','postcentral','supramarginal',...
    'parsopercularis','parstriangularis','parsorbitalis','rostralmiddlefrontal','caudalmiddlefrontal'};

% %causalpre
% invest_regions = {{'mMTG','rMTG',...
%     'inferiorprecentral','superiorprecentral','postcentral','supramarginal',...
%     'parsopercularis','parstriangularis','rostralmiddlefrontal','caudalmiddlefrontal'}};

%causal
% invest_regions = {{'cSTG','rSTG',...
%     'inferiorprecentral','superiorprecentral','postcentral',...
%     'parsopercularis','parstriangularis','rostralmiddlefrontal','caudalmiddlefrontal'}};

%anticausal
% invest_regions = {{'cSTG','mMTG',...
%     'inferiorprecentral','superiorprecentral','postcentral','supramarginal',...
%     'parsopercularis','parstriangularis','caudalmiddlefrontal'}};

% invest_regions = {{'cSTG','inferiorprecentral','superiorprecentral','postcentral','parsopercularis','parstriangularis','caudalmiddlefrontal'}};

% invest_regions = {{'mSTG'},...
%     {'inferiorprecentral','superiorprecentral','postcentral','supramarginal'},...
%     {'parsopercularis','parstriangularis','rostralmiddlefrontal','caudalmiddlefrontal'}};
invest_regions_legend = {{'cSTG','vPrCG','dPrCG','PsCG','pOp','pTri','cMFG'}};
ir = 1;
[~,investregionindex] = ismember(invest_regions{ir},all_regions);

% colors{1}=[givemecolor('Purples',2,0.6,1);givemecolor('Blues',3,0.4,1)];
% %       colors{1}=[givemecolor('RdPu',3,0.45,1);givemecolor('GnBu',3,0.45,1)];
% colors{2}=givemecolor('YlGn',4);
% colors{3}=givemecolor('OrRd',4);

colors{1}=[givemecolor('Purples',2,0.6,1);givemecolor('Blues',3,0.4,1);givemecolor('YlGn',4);givemecolor('OrRd',5)];
%       colors{1}=[givemecolor('RdPu',3,0.45,1);givemecolor('GnBu',3,0.45,1)];
% 
% %       
% if return_curve
%     fig1=figure();
% %     set(fig1,'position',[0,0,560/2,420/2])
%     rows = row-1;%length(attrs)+1;%1;
%     cols = 1;%length(attrs)+1;
%     haxis = tight_subplot(rows,cols,[0,0],[0.02,0],[0.1,0.1]);
%     fig2=figure();
%     set(fig2,'position',[0,0,2560,270])
%     TT=7;
%     haxis2 = tight_subplot(rows,TT,[0,0],0,0);
%     dumm = 1;
%     onregion = {'mSTG','rSTG','cSTG','rMTG','mMTG','cMTG','inferiorprecentral','superiorprecentral','postcentral','caudalmiddlefrontal','rostralmiddlefrontal','parstriangularis','parsopercularis','supramarginal'};
% %     colors = [[0 0.4470 0.7410 0.7];[0.3010 0.7450 0.9330 0.7];[0.8500 0.3250 0.0980 0.7];[0.9290 0.6940 0.1250 0.7];[0.4660 0.6740 0.1880 0.7];[0.4940 0.1840 0.5560 0.7]];
% %     Sat = [0.8,0.9,1];
% %     Val = [1,0.75,0.5];    
% %     for i=1:length(invest_regions)
% %         colors{i} = linspace(0,1,length(invest_regions{i})+1); colors{i} = colors{i}(1:end-1); 
% %         colors{i} = [colors{i}',Sat(i)*ones(length(colors{i}),1),Val(i)*ones(length(colors{i}),1)]; 
% %         colors{i} = hsv2rgb(colors{i});
% %     end
% %     colors{1} = [[0 0.4470 0.7410];[0.3010 0.7450 0.9330];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4660 0.6740 0.1880];[0.4940 0.1840 0.5560]];
% %     
%       color = givemecolor('YlGn',4); greencolor = color(end,:); startlinewidth=2;
% %     for i=1:length(invest_regions); colors{i} = [colors{i}',Sat(i)*ones(length(colors{i}),1),Val(i)*ones(length(colors{i}),1)]; colors{i} = hsl2rgb(colors{i});end
% %     colors = [colors',0.48*ones(length(colors),1),0.5*ones(length(colors),1)];
% %     colorsb = linspace(0,1,length(onregion)/2+1); colorsb = colorsb(1:end-1);
% %     colorsb = [colorsb',ones(length(colorsb),1),0.75*ones(length(colorsb),1)];
% %     colors = hsl2rgb(colors);
% %     colors = [[94 116 194];[187 118 196];[253 125 171];[255 153 134];[255 198 105];[249 248 113];[181 181 181]];
% %     colors = [colors/2;colors]/256;
%     win=hann(11); win = win./sum(win);
%     for a = 1:length(attrs)
%         [frames,frames_brain,err,tticks] = attr_vis_tau(attrs{a},post,visdiff,postabs,max_att,plot_on_same_brain,above_median,temp,wreg,entire_period,true,invest_regions{ir});
%         if strcmp(attrs{a},'freq_formants_hamon')
%             for j=1:2
%                 axes(haxis(sub2ind([rows,cols],dumm)));
%                 if tticks(1)==0
%                     reflectdir = 'right';
%                 else
%                     reflectdir = 'left';
%                 end
%                 smoothcurve = squeeze(frames(:,:,j));
%                 smootherr = squeeze(err(:,:,j));
% %                 smoothcurve = myconv(squeeze(frames(:,:,j)),reshape(win,size(win,1),1),reflectdir);
% %                 smootherr = myconv(squeeze(err(:,:,j)),reshape(win,size(win,1),1),reflectdir);
%                 baseline = repmat(reshape(frames(1,:,j),[1,size(frames,2)]),[size(frames,1),1]);
%                 figure(fig1);p=plot(tticks,smoothcurve,'LineWidth',4);
%                 ax=gca;ax.FontSize=20;  set(gca,'linewidth',2)
% %                 figure(fig1);p=plot(frames(end:-1:1,:,j),'LineWidth',3);
%                 ylim(haxis(sub2ind([rows,cols],dumm)),[ylims_min(a),ylims_max(a)]);
%                 set(haxis(sub2ind([rows,cols],dumm)),'box','off','TickDir','out');
%                 xlim(haxis(sub2ind([rows,cols],dumm)),[tticks(1),tticks(end)]);
%                 pbaspect([1.5 1 1])
%                 for x = 1:size(frames,2); p(x).Color=colors{1}(investregionindex(x),1:3); p(x).Color(4)=0.7; end
%                 hold on;
%                 for x = 1:size(frames,2); fill([tticks,fliplr(tticks)] , [smoothcurve(:,x,:)'-smootherr(:,x,:)',fliplr(smoothcurve(:,x,:)'+smootherr(:,x,:)')],colors{1}(investregionindex(x),1:3),'FaceAlpha',0.3,'EdgeAlpha',0); end
%                 if tticks(1)==0
%                     ff = fill([tticks(1:startlinewidth),fliplr(tticks(1:startlinewidth))] , [ylims_min(a)*ones(1,startlinewidth),ylims_max(a)*ones(1,startlinewidth)],greencolor,'FaceAlpha',0.7,'EdgeAlpha',0);
%                 else
%                     set(haxis(sub2ind([rows,cols],dumm)),'YAxisLocation','right');
%                     ff = fill([tticks(end-startlinewidth+1:end),fliplr(tticks(end-startlinewidth+1:end))] , [ylims_min(a)*ones(1,startlinewidth),ylims_max(a)*ones(1,startlinewidth)],greencolor,'FaceAlpha',0.7,'EdgeAlpha',0); 
%                 end
% %                 legend('cSTG','mSTG','parstriangularis','parsopercularis','precentral','postcentral');
%                 frm_brain = frames_brain(:,:,j);
%                 for t=1:length(frm_brain)
%                     axcp = copyobj(frm_brain{t}.CurrentAxes,fig2);
%                     set(axcp,'Position',get(haxis2(sub2ind([TT,rows],t,dumm)),'position'));
% %                     set(axcp,'Position',get(haxis2(sub2ind([TT,rows],TT-t+1,dumm)),'position'));
%                     close(frm_brain{t})
%                     set(haxis2(sub2ind([TT,rows],t,dumm)),'visible','off');
%                 end
% %                 switch j
% %                     case 1
% %                         ylabel(haxis(sub2ind([rows,cols],dumm)),'f_1');
% %                     case 2
% %                         ylabel(haxis(sub2ind([rows,cols],dumm)),'f_2');
% %                     case 3
% %                         ylabel(haxis(sub2ind([rows,cols],dumm)),'f_2 - f_1');
% %                 end
%                 dumm=dumm+1;
%             end
%             
%         else
%             axes(haxis(sub2ind([rows,cols],dumm)));
%             if tticks(1)==0
%                 reflectdir = 'right';
%             else
%                 reflectdir = 'left';
%             end
%             smoothcurve = myconv(squeeze(frames(:,:,:)),reshape(win,size(win,1),1),reflectdir);
%             smootherr = myconv(squeeze(err(:,:,:)),reshape(win,size(win,1),1),reflectdir);
%             figure(fig1);p=plot(tticks,smoothcurve,'LineWidth',4); 
% %             hold on; plot(tticks,noiselevel*ones(size(tticks)),'LineWidth',4,'Color',[0.5,0.5,0.5,0.7],'LineStyle',':');
%             ax=gca;ax.FontSize=20;  set(gca,'linewidth',2)
% %             figure(fig1);p=plot(frames(end:-1:1,:,:),'LineWidth',3);
%             ylim(haxis(sub2ind([rows,cols],dumm)),[ylims_min(a),ylims_max(a)]);
%             set(haxis(sub2ind([rows,cols],dumm)),'box','off','TickDir','out');
%             xlim(haxis(sub2ind([rows,cols],dumm)),[tticks(1),tticks(end)]);
%             pbaspect([1.5 1 1])
%             for x = 1:size(frames,2)-1; p(x).Color=colors{1}(investregionindex(x),1:3); p(x).Color(4)=0.7; end
%             p(end).Color=[0.75,0.75,0.75]; p(x).Color(4)=0.7;
%             hold on;
%             for x = 1:size(frames,2)-1; ff = fill([tticks,fliplr(tticks)] , [smoothcurve(:,x,:)'-smootherr(:,x,:)',fliplr(smoothcurve(:,x,:)'+smootherr(:,x,:)')],colors{1}(investregionindex(x),1:3),'FaceAlpha',0.3,'EdgeAlpha',0); end
%             ff = fill([tticks,fliplr(tticks)] , [smoothcurve(:,end,:)'-smootherr(:,end,:)',fliplr(smoothcurve(:,end,:)'+smootherr(:,end,:)')],[0.75,0.75,0.75],'FaceAlpha',0.3,'EdgeAlpha',0);
%             
%             if tticks(1)==0
%                 ff = fill([tticks(1:startlinewidth),fliplr(tticks(1:startlinewidth))] , [ylims_min(a)*ones(1,startlinewidth),ylims_max(a)*ones(1,startlinewidth)],greencolor,'FaceAlpha',0.7,'EdgeAlpha',0);
%             else
%                 set(haxis(sub2ind([rows,cols],dumm)),'YAxisLocation','right');
%                 ff = fill([tticks(end-startlinewidth+1:end),fliplr(tticks(end-startlinewidth+1:end))] , [ylims_min(a)*ones(1,startlinewidth),ylims_max(a)*ones(1,startlinewidth)],greencolor,'FaceAlpha',0.7,'EdgeAlpha',0); 
%             end
% %             legend('cSTG','mSTG','parstriangularis','parsopercularis','precentral','postcentral');
%             for t=1:length(frames_brain)
%                 axcp = copyobj(frames_brain{t}.CurrentAxes,fig2);
%                 set(axcp,'Position',get(haxis2(sub2ind([TT,rows],t,dumm)),'position'));
% %                 set(axcp,'Position',get(haxis2(sub2ind([TT,rows],TT-t+1,dumm)),'position'));
%                 close(frames_brain{t})
%                 set(haxis2(sub2ind([TT,rows],t,dumm)),'visible','off');
%                 pos = get(haxis2(sub2ind([TT,rows],t,dumm)),'position');
%                 currentax=haxis2(sub2ind([TT,rows],t,dumm));
% %                 if tticks(1)==0
% %                     if t==1
% %                         rectangle('position',[pos(1), pos(2), pos(3), pos(4)],'LineWidth',10,'EdgeColor',[greencolor,0.7]);
% %                     end
% %                 else
% %                     if t==length(frames_brain)
% %                         rectangle('position',[pos(1), pos(2), pos(3), pos(4)],'LineWidth',10,'EdgeColor',[greencolor,0.7]);
% %                     end
% %                 end
%             end
%             
%             
% %             switch attrs{a}
% %                 case 'f0_hz'
% %                     ylabel(haxis(sub2ind([rows,cols],dumm)),'f_0');
% %                 case 'amplitudes'
% %                     ylabel(haxis(sub2ind([rows,cols],dumm)),'alpha');
% %                 case 'freq_harm_diff'
% %                     ylabel(haxis(sub2ind([rows,cols],dumm)),'f_2 / f_1');
% %                 case 'freq_formants_noise'
% %                     ylabel(haxis(sub2ind([rows,cols],dumm)),'noise frequency');
% %                 case 'bandwidth_formants_noise_hz'
% %                     ylabel(haxis(sub2ind([rows,cols],dumm)),'noise bandwidth');
% %                 otherwise
% %                     ylabel(haxis(sub2ind([rows,cols],dumm)),attrs{a});
% %             end
%             dumm=dumm+1;
%         end
%         if dumm==2
% %             legend('rSTG','mSTG','cSTG','rMTG','mMTG','cMTG','inferiorprecentral','superiorprecentral','postcentral','caudalmiddlefrontal','rostralmiddlefrontal','parstriangularis','parsopercularis','supramarginal');
% %             legend('rSTG','mSTG','cSTG','rMTG','mMTG','cMTG');
% %             legend('inferiorprecentral','superiorprecentral','postcentral','supramarginal');
% %             legend('caudalmiddlefrontal','rostralmiddlefrontal','parstriangularis','parsopercularis');
% %             if ir==1
% %                 lgd=legend({'rSTG','cSTG'});    
% %             else
% %                 lgd=legend(invest_regions_legend{ir});     
% %             end
%             lgd=legend(invest_regions_legend{ir});    
%             lgd.FontSize=20;
%         end
%     end
% else
%     if tau
%         for a = 1:length(attrs)
%             frames = attr_vis_tau(attrs{a},post,visdiff,postabs,max_att,plot_on_same_brain,above_median,temp,wreg,entire_period,false);
%         %     frames = attr_vis_occbased(attrs{a},post,visdiff,postabs,max_att,plot_on_same_brain,above_median,temp,wreg);
%             if strcmp(attrs{a},'freq_formants_hamon')
%                 for j=1:3
%                     figure();
%                     haxis = tight_subplot(N_tau,col,[0,0],0,0);
%                     frame = frames(:,:,:,(j-1)*3+1:j*3);
%                     frame = reshape(frame,[col,N_tau,size(frame,2),size(frame,3),3]);
%                     for t =1:size(frame,1)
%                         for tau_ =1:N_tau
%                             axes(haxis(sub2ind([col,N_tau],t,tau_)));
%                             if tau_>=t
%                                 imshow(uint8(squeeze(frame(t,tau_,:,:,:))));
%                             else
%                                 imshow(uint8(squeeze(frame(1,1,:,:,:))));
%                             end
%                             axis off;
%                         end
%                     end
%                     linkaxes(haxis,'xy');
%                     switch j
%                         case 1
%                             mtit('f_1');
%                         case 2
%                             mtit('f_2');
%                         case 3
%                             mtit('f_2 - f_1');
%                     end
%                 end
% 
%             else
%                 figure();
%                 haxis = tight_subplot(N_tau,col,[0,0],0,0);
%                 frame = frames;
%                 frame = reshape(frame,[col,N_tau,size(frame,2),size(frame,3),3]);
%                 for t =1:size(frame,1)
%                     for tau_ =1:N_tau
%                         axes(haxis(sub2ind([col,N_tau],t,tau_)));
%                         if tau_>=t
%                             imshow(uint8(squeeze(frame(t,tau_,:,:,:))));
%                         else
%                             imshow(uint8(squeeze(frame(1,1,:,:,:))));
%                         end
%                         axis off;
%                     end
%                 end
%                 switch attrs{a}
%                     case 'f0_hz'
%                         mtit('f_0');
%                     case 'amplitudes'
%                         mtit('alpha');
%                     case 'freq_harm_diff'
%                         mtit('f_2 / f_1');
%                     case 'freq_formants_noise'
%                         mtit('noise frequency');
%                     case 'bandwidth_formants_noise_hz'
%                         mtit('noise bandwidth');
%                     otherwise
%                         mtit(attrs{a});
%                 end
%             end
%             linkaxes(haxis,'xy');
%         end
% 
% 
%     else
%         figure();
%         haxis = tight_subplot(row,col,[0,0],0,0);
%         dumm = 1;
%         for a = 1:length(attrs)
%             frames = attr_vis(attrs{a},post,visdiff,postabs,max_att,plot_on_same_brain,above_median,temp,wreg,entire_period);
%         %     frames = attr_vis_occbased(attrs{a},post,visdiff,postabs,max_att,plot_on_same_brain,above_median,temp,wreg);
%             if strcmp(attrs{a},'freq_formants_hamon')
%                 for j=1:3
%                     frame = frames(:,:,:,(j-1)*3+1:j*3);
%                     for t =1:size(frame,1)
%                         axes(haxis(sub2ind([col,row],t,dumm)));
%                         imshow(uint8(squeeze(frame(t,:,:,:))));
%                         axis off;
%                     end
%                     switch j
%                         case 1
%                             ylabel(haxis(sub2ind([col,row],1,dumm)),'f_1');
%                         case 2
%                             ylabel(haxis(sub2ind([col,row],1,dumm)),'f_2');
%                         case 3
%                             ylabel(haxis(sub2ind([col,row],1,dumm)),'f_2 - f_1');
%                     end
%                     dumm = dumm+1;
%                 end
% 
%             else
%                 frame = frames;
%                 for t =1:size(frame,1)
%                     axes(haxis(sub2ind([col,row],t,dumm)));
%                     imshow(uint8(squeeze(frame(t,:,:,:))));
%                     axis off;
%                 end
%                 switch attrs{a}
%                     case 'f0_hz'
%                         ylabel(haxis(sub2ind([col,row],1,dumm)),'f_0');
%                     case 'amplitudes'
%                         ylabel(haxis(sub2ind([col,row],1,dumm)),'alpha');
%                     case 'freq_harm_diff'
%                         ylabel(haxis(sub2ind([col,row],1,dumm)),'f_2 / f_1');
%                     case 'freq_formants_noise'
%                         ylabel(haxis(sub2ind([col,row],1,dumm)),'noise frequency');
%                     case 'bandwidth_formants_noise_hz'
%                         ylabel(haxis(sub2ind([col,row],1,dumm)),'noise bandwidth');
%                     otherwise
%                         ylabel(haxis(sub2ind([col,row],1,dumm)),attrs{a});
%                 end
%                 dumm = dumm+1;
%             end
% 
%         end
%         linkaxes(haxis,'xy');
%     end
% end



fig=figure();
dumm = 1;
avg_cell = {0,0,0};
regioncolor=false;
bar_plot=false;
do_gaussian=true;
norm_contr_elec=false;
if bar_plot
    haxis = tight_subplot(row-1,col,[0,0],[0.3,0],[0.1,0]);
else
%     haxis = tight_subplot(row-1,col,[0,0],0,[0.1,0]); 
    haxis = tight_subplot(row-1,col,[0,0],0,[0,0]); 
end
for a = 1:length(attrs)
%     if strcmp(attrs{a},'avg')
%         for sub=1:length(avg)
%             avg_cell{sub} = avg_cell{sub}/(size(attrs,2)-3);
% %             avg_cell{sub} = avg_cell{sub}/(size(attrs,2)-3);
%         end
%     end
    [frames1,frames2,frames,avg] = attr_vis_contrast(attrs{a},post,visdiff,postabs,max_att,plot_on_same_brain,above_median,temp,wreg,entire_period,avg_cell,regioncolor,bar_plot,norm_contr_elec);
%     [frames,avg] = attr_vis(attrs{a},post,visdiff,postabs,max_att,plot_on_same_brain,above_median,temp,wreg,entire_period,avg_cell,regioncolor); frames1=frames; frames=frames;
    switch attrs{a}
        case 'freq_formants_hamon'
            for j=1:2
                if bar_plot || do_gaussian
                    frame = frames(1,1,1,j);
                else
                    frame = frames(:,:,:,(j-1)*3+1:j*3);
                end
                for t =1:size(frame,1)
                    ax=haxis(sub2ind([col,row],t,dumm));
                    if bar_plot || do_gaussian
                        if do_gaussian
                            fr = frame{t};
                            axcp = copyobj(fr.CurrentAxes,fig);
                            set(axcp,'Position',get(ax,'position'));
                            axis off;
                            frame = frames(1,1,1,j);
                            try
                            frame1 = frames1(1,1,1,j);
                            frame2 = frames2(1,1,1,j);
                            end
                            close(frame{t});
                            try
                            close(frame1{t});
                            close(frame2{t});
                            end
                        else
                            fr = frame;
                            axcp = copyobj(fr.CurrentAxes.Children,ax);
                            ax.YLim=[-0.6,0.6];
                            close(fr);
                            try
                            close(frames1(1,1,1,j));
                            close(frames2(1,1,1,j));
                            end
                        end
                        
                    else
                        imshow(uint8(squeeze(frame(t,:,:,:))));
                        axis off;
                    end
                end
                switch j
                    case 1
                        ylabel(haxis(sub2ind([col,row],1,dumm)),'f_1');
                    case 2
                        ylabel(haxis(sub2ind([col,row],1,dumm)),'f_2');
                end
                dumm = dumm+1;
            end
            for sub=1:length(avg)
                avg_cell{sub} = avg_cell{sub}+avg{sub};
            end
%         case 'amplitudes'
%             for j=1:2
%                 frame = frames(:,:,:,(j-1)*3+1:j*3);
%                 for t =1:size(frame,1)
%                     axes(haxis(sub2ind([col,row],t,dumm)));
%                     imshow(uint8(squeeze(frame(t,:,:,:))));
%                     axis off;
%                 end
%                 switch j
%                     case 1
%                         ylabel(haxis(sub2ind([col,row],1,dumm)),'alpha_harm');
%                     case 2
%                         ylabel(haxis(sub2ind([col,row],1,dumm)),'alpha_noise');
%                 end
%                 dumm = dumm+1;
%             end
%             
        otherwise
            frame = frames;
            for t =1:size(frame,1)
                if strcmp(attrs{a},'avg')
                    ax=haxis(sub2ind([col,row],t,1));
                else
                    ax=haxis(sub2ind([col,row],t,dumm));
                end
                if bar_plot || do_gaussian
                    if a==length(attrs)
                        if do_gaussian
                            fr = frame{t};
                            axis off;
                        else
                            fr = frame;
                        end
                        axcp = copyobj(fr.CurrentAxes,fig);
                        set(axcp,'Position',get(ax,'position'));
%                         axcp.YLim=[-0.4,0.4];
                        try
                        close(frames{t});
                        end
                        try
                        close(frames1{t});
                        close(frames2{t});
                        end
                        
                    else
                        if do_gaussian
                            fr = frame{t};
                            axcp = copyobj(fr.CurrentAxes,fig);
                            set(axcp,'Position',get(ax,'position'));
                            axis off;
                            close(frames{t});
                            try
                                close(frames1{t});
                                close(frames2{t});
                            end
                        else
                            fr = frame;
                            axcp=copyobj(fr.CurrentAxes.Children,ax);
                            ax.YLim=[-0.6,0.6];
                            close(fr);
                            try
                            close(frames1);
                            close(frames2);
                            end
                        end
                        
                    end
                else
                    imshow(uint8(squeeze(frame(t,:,:,:))));
                    axis off;
                end
            end
            if ~strcmp(attrs{a},'spec')
                for sub=1:length(avg)
                    avg_cell{sub} = avg_cell{sub}+avg{sub};
                end
            end
%             switch attrs{a}
%                 case 'f0_hz'
%                     ylabel(haxis(sub2ind([col,row],1,dumm)),'f_0');
%                 case 'amplitudes'
%                     ylabel(haxis(sub2ind([col,row],1,dumm)),'alpha');
%                 case 'freq_harm_diff'
%                     ylabel(haxis(sub2ind([col,row],1,dumm)),'f_2 / f_1');
%                 case 'freq_formants_noise'
%                     ylabel(haxis(sub2ind([col,row],1,dumm)),'noise frequency');
%                 case 'bandwidth_formants_noise_hz'
%                     ylabel(haxis(sub2ind([col,row],1,dumm)),'noise bandwidth');
%                 case 'avg'
%                     ylabel(haxis(sub2ind([col,row],1,1)),'avg');
%                 otherwise
%                     ylabel(haxis(sub2ind([col,row],1,dumm)),attrs{a});
%             end
            dumm = dumm+1;
    end
%     if a==length(attrs)
%         delete(ax)
%     end
end
set(gca,'visible','off')

linkaxes(haxis,'xy');
region_contrast();
% ---------------------------------------------------------------------------
% ----------------------------Relevant Functions ----------------------------
% ---------------------------------------------------------------------------
function y = myconv(a,b,reflectdir)
    filter_size=floor(size(b,1)/2);
    switch reflectdir
        case 'left'
            a = [a(filter_size:-1:1,:);a;zeros(size(a(end:-1:end-filter_size+1,:)))];
        case 'right'
            a = [zeros(size(a(filter_size:-1:1,:)));a;a(end:-1:end-filter_size+1,:)];
    end
%     switch reflectdir
%         case 'left'
%             a = [a(filter_size:-1:1,:);a;2*a(end,:)-a(end:-1:end-filter_size+1,:)];
%         case 'right'
%             a = [2*a(1,:)-a(filter_size:-1:1,:);a;a(end:-1:end-filter_size+1,:)];
%     end
%     switch reflectdir
%         case 'left'
%             a = [2*a(1)-a(filter_size:-1:1,:);a;2*a(end,:)-a(end:-1:end-filter_size+1,:)];
%         case 'right'
%             a = [2*a(a)-a(filter_size:-1:1,:);a;2*a(end,:)-a(end:-1:end-filter_size+1,:)];
%     end
    y = conv2(a,b,'valid');
end

function [clr] = givemecolor(name,num,min,max)
if ~exist('min','var')
    min=0.35;
end
if ~exist('max','var')
    max=1;
end
[color_map_]=getPyPlot_cMap(name, 256);
if num==1
    clr = color_map_(int64(256*max),:);
else
    clr = color_map_(int64(linspace(256*max,256*min,num)),:);
end

end

function [fig] = plot_surface(surf_data, elec_data)
    if (isempty(surf_data.annot) || ~exist(surf_data.annot,'file'))
        surf_data.surf_color = repmat([0.7,0.7,0.7], [size(surf_data.surf_brain.coords,1),1]);
    else
        [~,albl,actbl]=fs_read_annotation(surf_data.annot);
        [~,cc] = ismember(albl,actbl.table(:,5));
        cc(cc==0) = 1;
        surf_data.surf_color = actbl.table(cc,1:3)./255;
    end
    for i = 1:size(elec_data.elec)
        near_ind = get_near_verts(surf_data.surf_brain.coords, elec_data.elec(i,:), 3);
        surf_data.surf_color(near_ind, :) = repmat([0.935, 0, 0], [size(near_ind,1),1]);
    end
    plot
    fig = figure();
    trisurf(surf_data.surf_brain.faces,...
            surf_data.surf_brain.coords(:,1),...
            surf_data.surf_brain.coords(:,2),...
            surf_data.surf_brain.coords(:,3),...
            'FaceVertexCData', surf_data.surf_color,...
            'FaceColor', 'interp',...
            'FaceAlpha',1);
    shading interp;
    lighting gouraud;
    material dull;
    light;
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1);
    bottom = outerpos(2); 
    ax_width = outerpos(3); 
    ax_height = outerpos(4);
    ax.Position = [left bottom ax_width ax_height];
    axis off;
    set(light,'Position',[-1,0,1]);
    view(270,0);
    set(fig,'color','white','InvertHardCopy','off');
    axis tight;
    axis equal;
end
function [near_ind] = get_near_verts(verts,elec,thresh)
    dist = (verts - repmat(elec,[size(verts,1),1])).^2;
    dist = sqrt(sum(dist,2));
    near_ind = find(dist<=thresh);
end
function [Att, mask, coord] = fake_data(elec)
    Att = rand(20,36,15,15);
    coord = zeros(3,15,15);
    coord(1,1:2:end,1:2:end) = reshape(elec(1:64,1),1,8,8);
    coord(2,1:2:end,1:2:end) = reshape(elec(1:64,2),1,8,8);
    coord(3,1:2:end,1:2:end) = reshape(elec(1:64,3),1,8,8);
    mask = zeros(15,15);
    mask(1:2:end,1:2:end) = 1;
end

% function [frames] = VisualAtt(Att, Mask, coord, brain,color_map,maxv,alpha,div)
%     
%     if nargin<5
%         color_map = 'hot';
%         color_map = eval([color_map,'(256)']);
%     end
%     att=[];
%     mask=[];
%     cod = [];
%     for sub = 1:length(Att)
%         att = cat(3,att,Att{sub});
%         mask = cat(1,mask,Mask{sub});
%         cod = cat(2,cod,coord{sub});
%     end
% %     mask = Mask{sub};
%     if length(size(att))==4
%         sizes = size(att);
%         Att_avg = reshape(mean(att,1),sizes(2:end));
%     else 
%         Att_avg = att;
%     end
%     if nargin<6
%         maxv = max(Att_avg(:));
%     end
% %     Att_avg = (Att_avg-min(Att_avg(:)))/(max(Att_avg(:))-min(Att_avg(:)));
%     Att_avg = (Att_avg)/maxv;
% %     mask(Att_avg<0.2)=0;
%     % Att_avg(1:6,:,:) = Att_avg(1:6,:,:)*0;
% %     if length(size(att))==3
% %         time_step = size(Att_avg,1);
% %     else
% %         time_step = 1;
% %     end
%     time_step = size(Att_avg,1);
%     frames = [];
%     for i=1:time_step
% %         if length(size(att))==3
% %             Att_frame = squeeze(Att_avg(i,:,:));
% %         else
% %             Att_frame = Att_avg;
% %         end
%         Att_frame = squeeze(Att_avg(i,:,:));
%         ind_active = find(mask==1);
% %         color_map = 'hot';
%         if div
%             color_ind = value2colorind(Att_frame(ind_active), 'hot',[-1,1]);
%         else
%             color_ind = value2colorind(Att_frame(ind_active), 'hot',[0,1]);
%         end
% %             colors = {color_ind, color_map};
%         if size(color_map,1)<100
%             if size(color_map,1)==1
%                 colors = {color_map};
%             else
%                 colors = {color_map(sub,:)};
%             end
%         else
%             colors = {color_ind, color_map};
%         end
%         elecname = 1:length(ind_active);
%         elec = zeros(length(ind_active),3);
%         Temp = squeeze(cod(1,:,:));
%         elec(:,1) = Temp(ind_active);
%         Temp = squeeze(cod(2,:,:));
%         elec(:,2) = Temp(ind_active);
%         Temp = squeeze(cod(3,:,:));
%         elec(:,3) = Temp(ind_active);
% %         fig = nyu_plot_whitebackground(brain.surf_brain,...
% %                                        brain.sph,...
% %                                        elec, elecname,...
% %                                        colors, 0, brain.annot, 2, 1);
%         
%         if alpha
%             if div
% %                 point_alpha = 2*abs(max(min(Att_frame(ind_active),1),0)-0.5);
%                 point_alpha = max(min(abs(Att_frame(ind_active)),1),0);
%             else
%                 point_alpha = max(min(Att_frame(ind_active),1),0);
%             end
%         else
% %             point_alpha = double(abs(Att_frame(ind_active))>0.3);%ones(size(Att_frame(ind_active)));
%             point_alpha = ones(size(Att_frame(ind_active)));
%         end
%         
%         fig = nyu_plot_whitebackground(brain.surf_brain,...
%                                        brain.sph,...
%                                        elec, elecname,...
%                                        colors, 0, [], 2, 1, 1, point_alpha);
% %         fig = nyu_plot_whitebackground(brain.surf_brain,...
% %                                        brain.sph,...
% %                                        elec, elecname,...
% %                                        colors, 0, [], 2, 1, 1);
%         frame = getframe(fig);
%         frames(i,:,:,:) = frame.cdata;
%         close(fig);
%     end
% %     figure();
% %     if size(frames,1)==1
% %         imagesc(uint8(squeeze(frames(1,:,:,:))));
% %         axis off;
% %     else
% %         haxis = tight_subplot(6,6,[0,0],0,0);
% %         for i=1:36
% %             axes(haxis(i));
% %             imagesc(uint8(squeeze(frames(i,:,:,:))));
% %             axis off;
% %         end
% %     end
%     
% end

% 
% function fig = barplot(data,clrmap,regions)
%     figg = figure;
%     clrp1 = clrmap(int32(length(clrmap)*0.75),:);
%     clrp2 = clrmap(int32(length(clrmap)*1),:);
%     clrn1 = clrmap(int32(length(clrmap)*0.25),:);
%     clrn2 = clrmap(1,:);
%     d_pos = {};
%     d_neg = {};
%     d_vec = [];
%     for m =1:length(data)
%         d_vec = [d_vec,data{m}];
%     end
%     [d_vec_sorted,ind] = sort(d_vec);
%     figure();b=bar(d_vec(ind));set(gca,'box','off','TickDir','out'); b.FaceColor=[0.5,0.5,0.5];b.EdgeColor=[0.25,0.25,0.25]; set(gca,'box','off','TickDir','out');ax=gca;ax.FontSize=20;  set(gca,'linewidth',2);set(gca,'FontWeight','bold');
%     figure();h=histogram(d_vec(ind),15);h.Normalization='probability';xlim([-1.1,1]);h.FaceColor=[0.5,0.5,0.5];h.EdgeColor=[0.25,0.25,0.25];set(gca,'box','off','TickDir','out');ax=gca;ax.FontSize=20;  set(gca,'linewidth',2);set(gca,'FontWeight','bold');
% %     for m =1:length(data)
% %         d = data{m};
% %         pos = d(d>=0);
% %         d_pos{m} = pos;
% %         pos_mean(m) = mean(pos)*length(pos)/length(d);
% %         pos_sterr(m) = std(pos)/sqrt(50*length(d))*1.96;
% %         neg = d(d<=0);
% %         d_neg{m} = neg;
% %         neg_mean(m) = mean(neg)*length(neg)/length(d);
% %         neg_sterr(m) = std(neg)/sqrt(50*length(d))*1.96;
% %     end
% %     hold on;
% %     alpha=1;
% % %     for ii=1:length(data); disp(signrank(d_pos{ii},d_neg{ii}));end;
% %     b=bar(pos_mean);b.FaceColor=clrp1;b.EdgeColor=clrp2; b.EdgeAlpha=alpha; ylim([-0.6,0.6]); b.FaceAlpha=alpha;
% %     er=errorbar(pos_mean,pos_sterr);er.LineStyle='none';er.LineWidth=1; er.Color=[clrp2,alpha]; er.CapSize=0;
% %     set([er.Bar, er.Line], 'ColorType', 'truecoloralpha', 'ColorData', [er.Line.ColorData(1:3); 255*alpha])
% %     b=bar(neg_mean);b.FaceColor=clrn1;b.EdgeColor=clrn2; b.EdgeAlpha=alpha; ylim([-0.6,0.6]); b.FaceAlpha=alpha;
% %     er=errorbar(neg_mean,neg_sterr);er.LineStyle='none';er.LineWidth=1; er.Color=[clrn2,alpha]; er.CapSize=0;
% %     set([er.Bar, er.Line], 'ColorType', 'truecoloralpha', 'ColorData', [er.Line.ColorData(1:3); 255*alpha])
%     neglim=-Inf*ones(1,length(data));
%     poslim=Inf*ones(1,length(data));poslim(7)=1;
%     for m =1:length(data)
%         d = data{m};
%         d(d<neglim(m))=neglim(m);
%         d(d>poslim(m)) = poslim(m);
% %         d = d(d>neglim(m) & d<poslim(m));
%         d_mean = mean(d);
%         if d_mean>0
%             pos = d;
%             neg = [];
%         else
%             pos = [];
%             neg = d;
%         end
%         pos_mean(m) = mean(pos);
%         pos_sterr(m) = std(pos)/sqrt(50*length(pos))*1.96*2;
%         neg_mean(m) = mean(neg);
%         neg_sterr(m) = std(neg)/sqrt(50*length(neg))*1.96*2;
%     end
%     figure();
%     hold on;
%     b=bar(pos_mean);b.FaceColor=clrp1; b.EdgeColor=clrp2;ylim([-0.5,0.3]); 
%     er=errorbar(pos_mean,pos_sterr);er.LineStyle='none';er.LineWidth=1; er.CapSize=0; er.Color=clrp2;
%     b=bar(neg_mean);b.FaceColor=clrn1; b.EdgeColor=clrn2; ylim([-0.5,0.3]); 
%     er=errorbar(neg_mean,neg_sterr);er.LineStyle='none';er.LineWidth=1; er.CapSize=0; er.Color=clrn2;
%     
%     [found,ind]=ismember('inferiorprecentral',regions);
%     regions{ind} = 'ventralprecentral';
%     [found,ind]=ismember('superiorprecentral',regions);
%     regions{ind} = 'dorsalprecentral';
%     [found,ind]=ismember('rSTG',regions);
%     regions{ind} = 'rostral STG';
%     [found,ind]=ismember('cSTG',regions);
%     regions{ind} = 'caudal STG';
%     [found,ind]=ismember('cMTG',regions);
%     regions{ind} = 'caudal MTG';
%     [found,ind]=ismember('mMTG',regions);
%     regions{ind} = 'middle MTG';
%     [found,ind]=ismember('rMTG',regions);
%     regions{ind} = 'rostral MTG';
%     xticks([1:length(data)]);xticklabels(regions);xtickangle(45); set(gca,'fontweight','bold');
%     
%     fig = ancestor(figg, 'figure');
% end



function fig = barplot(data,clrmap,regions)
    figg = figure;
    clrp1 = clrmap(int32(length(clrmap)*0.75),:);
    clrp2 = clrmap(int32(length(clrmap)*1),:);
    clrp3 = clrmap(int32(length(clrmap)*0.85),:);
    clrn1 = clrmap(int32(length(clrmap)*0.25),:);
    clrn2 = clrmap(1,:);
    clrn3 = clrmap(int32(length(clrmap)*0.15),:);
    d_pos = {};
    d_neg = {};
    d_vec = [];
    for m =1:length(data)
        d_vec = [d_vec,data{m}];
    end
    [d_vec_sorted,ind] = sort(d_vec);
    figure();b=bar(d_vec(ind));set(gca,'box','off','TickDir','out'); b.FaceColor=[0.5,0.5,0.5];b.EdgeColor=[0.25,0.25,0.25]; set(gca,'box','off','TickDir','out');ax=gca;ax.FontSize=20;  set(gca,'linewidth',2);set(gca,'FontWeight','bold');
    figure();h=histogram(d_vec(ind),15);h.Normalization='probability';xlim([-1.1,1]);h.FaceColor=[0.5,0.5,0.5];h.EdgeColor=[0.25,0.25,0.25];set(gca,'box','off','TickDir','out');ax=gca;ax.FontSize=20;  set(gca,'linewidth',2);set(gca,'FontWeight','bold');
%     for m =1:length(data)
%         d = data{m};
%         pos = d(d>=0);
%         d_pos{m} = pos;
%         pos_mean(m) = mean(pos)*length(pos)/length(d);
%         pos_sterr(m) = std(pos)/sqrt(50*length(d))*1.96;
%         neg = d(d<=0);
%         d_neg{m} = neg;
%         neg_mean(m) = mean(neg)*length(neg)/length(d);
%         neg_sterr(m) = std(neg)/sqrt(50*length(d))*1.96;
%     end
%     hold on;
%     alpha=1;
% %     for ii=1:length(data); disp(signrank(d_pos{ii},d_neg{ii}));end;
%     b=bar(pos_mean);b.FaceColor=clrp1;b.EdgeColor=clrp2; b.EdgeAlpha=alpha; ylim([-0.6,0.6]); b.FaceAlpha=alpha;
%     er=errorbar(pos_mean,pos_sterr);er.LineStyle='none';er.LineWidth=1; er.Color=[clrp2,alpha]; er.CapSize=0;
%     set([er.Bar, er.Line], 'ColorType', 'truecoloralpha', 'ColorData', [er.Line.ColorData(1:3); 255*alpha])
%     b=bar(neg_mean);b.FaceColor=clrn1;b.EdgeColor=clrn2; b.EdgeAlpha=alpha; ylim([-0.6,0.6]); b.FaceAlpha=alpha;
%     er=errorbar(neg_mean,neg_sterr);er.LineStyle='none';er.LineWidth=1; er.Color=[clrn2,alpha]; er.CapSize=0;
%     set([er.Bar, er.Line], 'ColorType', 'truecoloralpha', 'ColorData', [er.Line.ColorData(1:3); 255*alpha])
    pos_cell=[];
    neg_cell=[];
    all_cell=[];
    all_name=[];
    pos_name=[];
    neg_name=[];
    [found,ind]=ismember('inferiorprecentral',regions);
    regions{ind} = 'ventralprecentral';
    [found,ind]=ismember('superiorprecentral',regions);
    regions{ind} = 'dorsalprecentral';
    [found,ind]=ismember('rSTG',regions);
    regions{ind} = 'rostral STG';
    [found,ind]=ismember('cSTG',regions);
    regions{ind} = 'caudal STG';
    [found,ind]=ismember('cMTG',regions);
    regions{ind} = 'caudal MTG';
    [found,ind]=ismember('mMTG',regions);
    regions{ind} = 'middle MTG';
    [found,ind]=ismember('rMTG',regions);
    regions{ind} = 'rostral MTG';
    neglim=[-0.6,0,-0.2,-0.2,-0.55,-1,-1, -1 ,-0.5,-0.5,-0.5,-0.8,-0.1];
    poslim=[0.6, 1,0.6, 0.6, 0.55,0.9,1, 0.9,0.5,0.5,0.8,0.5,0.5];
    for m =1:length(data)
        d = data{m};    
        d_mean = mean(d);
        dd=d;
%         dd = d(d>neglim(m) & d<poslim(m));
        if 1% d_mean>0
            pos = d;
            neg = [];
            pos_cell=[pos_cell,dd];
            pos_name=[pos_name,repmat(regions(m),1,size(dd,2))];
            neg_cell = [neg_cell,nan];
            neg_name=[neg_name,regions(m)];
        else
            pos = [];
            neg = d;
            pos_cell=[pos_cell,nan];
            pos_name=[pos_name,regions(m)];
            neg_cell = [neg_cell,dd];
            neg_name=[neg_name,repmat(regions(m),1,size(dd,2))];
        end
        all_name=[all_name,repmat(regions(m),1,size(d,2))];
        all_cell=[all_cell,d];
        pos_mean(m) = mean(pos);
        pos_sterr(m) = std(pos)/sqrt(50*length(pos))*1.96;
        neg_mean(m) = mean(neg);
        neg_sterr(m) = std(neg)/sqrt(50*length(neg))*1.96;
    end
    figure();
    hold on;
%     b=bar(pos_mean);b.FaceColor=clrp1; b.EdgeColor=clrp2;ylim([-0.4,0.4]); 
%     er=errorbar(pos_mean,pos_sterr);er.LineStyle='none';er.LineWidth=1; er.CapSize=0; er.Color=clrp2;
%     b=bar(neg_mean);b.FaceColor=clrn1; b.EdgeColor=clrn2; ylim([-0.4,0.4]); 
%     er=errorbar(neg_mean,neg_sterr);er.LineStyle='none';er.LineWidth=1; er.CapSize=0; er.Color=clrn2;
%     xticks([1:length(data)]);xticklabels(regions);xtickangle(45); set(gca,'fontweight','bold');
%     
%     b = boxplot(all_cell,all_name,'PlotStyle','compact'); set(gca,'XTickLabel',{' '});
%     xticks([1:length(data)]);xticklabels(regions);xtickangle(45); set(gca,'fontweight','bold');
    boxplot(pos_cell,pos_name,'PlotStyle','compact','Colors',clrp3); set(gca,'XTickLabel',{' '}); ylim([-1,1]); set(gca,'box','off');
    boxplot(neg_cell,neg_name,'PlotStyle','compact','Colors',clrn3); set(gca,'XTickLabel',{' '}); ylim([-1,1]); set(gca,'box','off');
    xticks([1:length(data)]);xticklabels(regions);xtickangle(45); set(gca,'fontweight','bold');
    fig = ancestor(figg, 'figure');
end



function fig = barplot2(data,data2,clrmap,regions)
    figg = figure;
    clrp1 = clrmap(int32(length(clrmap)*0.75),:);
    clrp2 = clrmap(int32(length(clrmap)*1),:);
    clrn1 = clrmap(int32(length(clrmap)*0.25),:);
    clrn2 = clrmap(1,:);
   
    for m =1:length(data)
        pos = data{m};
        neg = -data2{m};
        pos_mean(m) = mean(pos);
        pos_sterr(m) = std(pos)/sqrt(50*length(pos))*1.96;
        neg_mean(m) = mean(neg);
        neg_sterr(m) = std(neg)/sqrt(50*length(neg))*1.96;
    end
    hold on;
    alpha=1;
    b=bar(pos_mean);b.FaceColor=clrp1;b.EdgeColor=clrp2; b.EdgeAlpha=alpha; ylim([-0.6,0.6]); b.FaceAlpha=alpha;
    er=errorbar(pos_mean,pos_sterr);er.LineStyle='none';er.LineWidth=1; er.Color=[clrp2,alpha]; er.CapSize=0;
    set([er.Bar, er.Line], 'ColorType', 'truecoloralpha', 'ColorData', [er.Line.ColorData(1:3); 255*alpha])
    b=bar(neg_mean);b.FaceColor=clrn1;b.EdgeColor=clrn2; b.EdgeAlpha=alpha; ylim([-0.6,0.6]); b.FaceAlpha=alpha;
    er=errorbar(neg_mean,neg_sterr);er.LineStyle='none';er.LineWidth=1; er.Color=[clrn2,alpha]; er.CapSize=0;
    set([er.Bar, er.Line], 'ColorType', 'truecoloralpha', 'ColorData', [er.Line.ColorData(1:3); 255*alpha])
    
    [found,ind]=ismember('inferiorprecentral',regions);
    regions{ind} = 'ventralprecentral';
    [found,ind]=ismember('superiorprecentral',regions);
    regions{ind} = 'dorsalprecentral';
    xticks([1:length(data)]);xticklabels(regions);xtickangle(45); set(gca,'fontweight','bold');
    
    fig = ancestor(figg, 'figure');
end

function fig = barplot3(data,data2,clrmap,regions)
    figg = figure;
    clrp1 = clrmap(int32(length(clrmap)*0.75),:);
    clrp2 = clrmap(int32(length(clrmap)*1),:);
    clrn1 = clrmap(int32(length(clrmap)*0.25),:);
    clrn2 = clrmap(1,:);
    pos=data;
    neg=-data2;
%     for m =1:length(data)
%         pos = data(m);
%         neg = -data2(m);
% %         pos_mean(m) = mean(pos);
% %         pos_sterr(m) = std(pos)/sqrt(50*length(pos))*1.96;
% %         neg_mean(m) = mean(neg);
% %         neg_sterr(m) = std(neg)/sqrt(50*length(neg))*1.96;
%     end
    hold on;
    [~,ind]=sort(pos+neg);
    pos = pos(ind); neg=neg(ind); %pos_mean=pos_mean(ind); pos_sterr=pos_sterr(ind); neg_mean=neg_mean(ind); neg_sterr=neg_sterr(ind);
    alpha=1;
    b=bar(pos);b.FaceColor=clrp1;b.EdgeColor=clrp2; b.EdgeAlpha=alpha; ylim([-1,1]); b.FaceAlpha=alpha;
%     er=errorbar(pos_mean,pos_sterr);er.LineStyle='none';er.LineWidth=1; er.Color=[clrp2,alpha]; er.CapSize=0;
%     set([er.Bar, er.Line], 'ColorType', 'truecoloralpha', 'ColorData', [er.Line.ColorData(1:3); 255*alpha])
    b=bar(neg);b.FaceColor=clrn1;b.EdgeColor=clrn2; b.EdgeAlpha=alpha; ylim([-1,1]); b.FaceAlpha=alpha;
%     er=errorbar(neg_mean,neg_sterr);er.LineStyle='none';er.LineWidth=1; er.Color=[clrn2,alpha]; er.CapSize=0;
%     set([er.Bar, er.Line], 'ColorType', 'truecoloralpha', 'ColorData', [er.Line.ColorData(1:3); 255*alpha])
    
%     [found,ind]=ismember('inferiorprecentral',regions);
%     regions{ind} = 'ventralprecentral';
%     [found,ind]=ismember('superiorprecentral',regions);
%     regions{ind} = 'dorsalprecentral';
%     xticks([1:length(data)]);xticklabels(regions);xtickangle(45); set(gca,'fontweight','bold');
    
    fig = ancestor(figg, 'figure');
end

% function [frames,frames_brain] = PlotAtt(Att1_cell, mask_cell,coord_cell,regions,causal,brain,colormap,Att_c,Att_cpre,Att_a,Att_p,ref)
% %     onregion = {{'cSTG'},{'mSTG'},{'parstriangularis'},{'parsopercularis'},{'precentral'},{'postcentral'}};
% %     col=13;
% %     N_tau=16;
% %     if causal
% %         t0=4;
% %         t1=8;
% %         tau1=0;
% %         tau2=8;
% %     else
% %         t0=12;
% %         t1=13;
% %         tau1=-6;
% %         tau2=0;
% %     end
% %     regions4barplot = {'rSTG','mSTG','cSTG','rMTG','mMTG','cMTG','inferiorprecentral','superiorprecentral','postcentral','caudalmiddlefrontal','rostralmiddlefrontal','parstriangularis','parsopercularis','supramarginal'};
% %     regions4barplot = {'rSTG','mSTG','cSTG','rMTG','mMTG','cMTG'};
% %     regions4barplot = {'inferiorprecentral','superiorprecentral','postcentral','supramarginal'};
%     regions4barplot = {'caudalmiddlefrontal','rostralmiddlefrontal','parstriangularis','parsopercularis'};
%     att=[];
%     att_ref = [];
%     mask=[];
%     region = [];
%     coord = [];
%     ratio = 1.5;%0.47;
%     switch ref
%         case 'causal'
%             Att_ref = Att_c;
%         case 'causalpre'
% %             Att_ref = Att_cpre;
%             Att_ref = Att_c;
%         case 'anticausal'
%             Att_ref = Att_a;
%         case 'percept'
%             Att_ref = Att_p;
%     end
%     for sub = 1:length(Att1_cell)
%         att = cat(3,att,Att1_cell{sub});
%         att_ref = cat(3,att_ref,Att_ref{sub});
%         mask = cat(1,mask,mask_cell{sub});
%         region = cat(1,region,regions{sub});
%         coord = cat(2,coord,coord_cell{sub});
%     end
%     col=13;
%     N_tau=16*8;
%     if causal
% %         t0=5;
% %         t1=8;
%         t0=1;
%         t1=5;
%         tau1=0;
%         tau2=7;
%     else
%         t0=12;
%         t1=13;
%         tau1=-6;
%         tau2=0;
%     end
% % % %     frames = zeros(tau2-tau1+1,length(regions4barplot));
%     frames = zeros(tau2*8-tau1*8+8,length(regions4barplot));
%     
%     att = reshape(att,[col,N_tau,size(att,3),size(att,4)]);
%     region_value=cell(length(regions4barplot),1);
%     
%     att_norm = 0;
%     for t=t0:t1
% %                 avg_data = avg_data+data_region(t,t+tau1:t+tau2);
%         att_norm = att_norm + att(t,8*(t-1)+1+tau1*8:8*t+tau2*8,:,:);
%     end
%     att_norm = att_norm/(t1-t0+1);
%     att_norm = mean(att_norm,2);
%     if ref~='causalpre'
%         att = att./(att_norm+1e-10).*att_ref*ratio;
%     else
%         for m=1:size(mask,1)
%             for n=1:size(mask,2)
%                 reg = reshape(region(m,n,:),1,[]);
%                 [found,ind]=ismember(reg(find(reg~=0)),{'rSTG','mSTG','cSTG'});
%                 if ~found && mask(m,n)
%     %                     region_count{ind} = region_count(ind)+1;
% %                     att(:,:,m,n) = att(:,:,m,n)./(att_norm(:,:,m,n)+1e-10).*att_ref(:,:,m,n)*ratio;
%                     att(:,:,m,n) = att(:,:,m,n)*0.7;
%                 else
%                     att(:,:,m,n) = att(:,:,m,n);
%                 end
%             end
%         end
% % 
%     end
% 
%     for m=1:size(mask,1)
%         for n=1:size(mask,2)
%             reg = reshape(region(m,n,:),1,[]);
%             [found,ind]=ismember(reg(find(reg~=0)),regions4barplot);
%             if found && mask(m,n)
% %                     region_count{ind} = region_count(ind)+1;
%                 for t=t0:t1
% %                 avg_data = avg_data+data_region(t,t+tau1:t+tau2);
%                     region_value{ind} = [region_value{ind};att(t,8*(t-1)+1+tau1*8:8*t+tau2*8,m,n)];
%                 end
%             end
%         end
%     end
%     for r = 1:length(regions4barplot)
%         frames(:,r) = mean(region_value{r},1);
%     end
%     
%     frames_brain = zeros(t1-t0+1,tau2*8-tau1*8+8,size(att,3),size(att,4));
%     mask_reshape = reshape(mask,1,1,size(att,3),size(att,4));
%     masked_att = mask_reshape.*att;
%     for t=t0:t1
%         frames_brain(t-t0+1,:,:,:) = masked_att(t,8*(t-1)+1+tau1*8:8*t+tau2*8,:,:);
%     end
%     frames_brain = reshape(frames_brain,size(frames_brain,1)*8,size(frames_brain,2)/8,size(frames_brain,3),size(frames_brain,4));
% %     frames_brain = mean(frames_brain,1);
%     frames_brain = sum(frames_brain,1)./(8*reshape([1,2,3,4,5,5,5,5],[1,8,1,1]));
%     [frames_brain] = VisualAtt({frames_brain}, {mask}, {region}, {coord}, brain,colormap,1,true,false,false);
% %     for r = 1:length(regions4barplot)
% %         elecs = 0;
% %         data_ofregion=0;
% %         for sub=1:length(Att1_cell)
% % %             data = Att1_cell{sub};
% %             data = abs(Att1_cell{sub});
% %             data = reshape(data,[col,N_tau,size(data,3),size(data,4)]);
% %             target = isregion(regions,regions4barplot{r});
% %             elecs = elecs+sum(target(:));
% %             targetrep = repmat(reshape(target,[1,1,15,15]),col,N_tau,1,1);
% %             data_region = data.*targetrep;
% %             data_region = sum(sum(data_region,4),3);
% %             avg_data = 0;
% %             for t=t0:t1
% % %                 avg_data = avg_data+data_region(t,t+tau1:t+tau2);
% %                 avg_data = avg_data+data_region(t,8*(t-1)+1+tau1*8:8*t+tau2*8);
% %             end
% %             denorm = ones(size(avg_data))*(t1-t0+1);
% % %             denorm(1:8)=1; denorm(9:16)=1; denorm(17:24)=2; denorm(25:32)=3;
% %             data_ofregion = data_ofregion+avg_data./denorm;
% %         end
% %         data_ofregion = data_ofregion/elecs;
% %         frames(:,r) = data_ofregion;
% %     end
% end


function [frames,frames_brain,err,tticks] = PlotAtt(Att1_cell, mask_cell,coord_cell,regions,causal,brain,colormap,Att_c,Att_cpre,Att_a,Att_p,Att_nt,ref,investreg)
%     onregion = {{'cSTG'},{'mSTG'},{'parstriangularis'},{'parsopercularis'},{'precentral'},{'postcentral'}};
%     col=13;
%     N_tau=16;
%     if causal
%         t0=4;
%         t1=8;
%         tau1=0;
%         tau2=8;
%     else
%         t0=12;
%         t1=13;
%         tau1=-6;
%         tau2=0;
%     end
%     regions4barplot = {'rSTG','mSTG','cSTG','rMTG','mMTG','cMTG','inferiorprecentral','superiorprecentral','postcentral','caudalmiddlefrontal','rostralmiddlefrontal','parstriangularis','parsopercularis','supramarginal'};
%     regions4barplot = {'rSTG','mSTG','cSTG','rMTG','mMTG','cMTG'};
%     regions4barplot = {'cSTG','rSTG','mMTG','cMTG','rMTG','inferiorprecentral','superiorprecentral','postcentral','supramarginal','parsopercularis','parstriangularis','rostralmiddlefrontal','caudalmiddlefrontal'};
    all_regions = {'cSTG','rSTG','mMTG','cMTG','rMTG',...
    'inferiorprecentral','superiorprecentral','postcentral','supramarginal',...
    'parsopercularis','parstriangularis','rostralmiddlefrontal','caudalmiddlefrontal'};
    name4plot = {'cSTG','rSTG','mMTG','cMTG','rMTG','vPreCG','dPreCG','PostCG','supramarginal','pOp','pTri','pOb','rMFG','cMFG'};
    regions4barplot =investreg;
%     regions4barplot = {'inferiorprecentral','superiorprecentral','postcentral','supramarginal'};
%     regions4barplot = {'parsopercularis','parstriangularis','rostralmiddlefrontal','caudalmiddlefrontal'};
    att=[];
    att_ref = [];
    att_nt = [];
    mask=[];
    region = [];
    coord = [];
    ratio = 1;%2.5;%0.47;
    switch ref
        case 'passive'
            Att_ref = Att_c;
        case 'causalpre'
            Att_ref = Att_cpre;
%             Att_ref = Att_c;
        case 'imagine'
            Att_ref = Att_a;
        case 'active'
            Att_ref = Att_p;
    end
    
    for sub = 1:length(Att1_cell)
        att = cat(3,att,Att1_cell{sub});
        att_ref = cat(3,att_ref,Att_ref{sub});
        att_nt = cat(3,att_nt,Att_nt{sub});
        mask = cat(1,mask,mask_cell{sub});
        region = cat(1,region,regions{sub});
        coord = cat(2,coord,coord_cell{sub});
    end
    col=144;%13;
    offset=16;
    rate = 8;
    N_tau=16;%16*rate;
    if causal
        if not(strcmp(ref,'causalpre'))
            tau0=5;
%             tau1=11;
            tau1=8;
        else
            tau0=7;
%             tau0=5;
            tau1=11;
        end
        t1=-6 + offset/rate;
        t2=0 + offset/rate;
    else
        tau0=5;
        tau1=7;
        t1=-1 + offset/rate;
        t2=5 + offset/rate;
    end
    switch ref
        case 'causal'
            tticks = [-55:0]/125*1000;
        case 'causalpre'
            tticks = [-55:0]/125*1000;
        case 'anticausal'
            tticks = [0:55]/125*1000;
        case 'percept'
            tticks = [0:55]/125*1000;
        otherwise
            tticks = [0:55]/125*1000;
    end
% % %     frames = zeros(tau2-tau1+1,length(regions4barplot));
    frames = zeros((t2-t1+1)*rate,length(regions4barplot)+1);
    err = zeros((t2-t1+1)*rate,length(regions4barplot)+1);
    
    att = reshape(att,[col,N_tau,size(att,3),size(att,4)]);
    att_nt = reshape(att_nt,[col,N_tau,size(att_nt,3),size(att_nt,4)]);
    region_value=cell(length(regions4barplot),1);
    region_samples=cell(length(regions4barplot),1);
    region_norm_value=cell(length(regions4barplot),1);
    region_ref_value=cell(length(regions4barplot),1);
    region_norm_samples=cell(length(regions4barplot),1);
    region_value_nt=cell(1,1);
    region_samples_nt=cell(1,1);
    region_samples_nt{1} = 0;
    for r=1:length(region_samples)
        region_samples{r} = 0;
    end
    
    att_norm = 0;
    att_norm_nt=0;
    for tau=tau0:tau1
%                 avg_data = avg_data+data_region(t,t+tau1:t+tau2);
        att_norm = att_norm + att(rate*(tau-1)+1+t1*rate:rate*tau+t2*rate,tau,:,:);
        att_norm_nt = att_norm_nt + att_nt(rate*(tau-1)+1+t1*rate:rate*tau+t2*rate,tau,:,:);
    end
    att_norm = att_norm/(tau1-tau0+1);
%     att_norm = max(att_norm,[],1);
    att_norm = mean(att_norm,1);
    att_norm_nt = att_norm_nt/(tau1-tau0+1);
    att_norm_nt = mean(att_norm_nt,1);
%     att = att*50;
    for m=1:size(mask,1)
        for n=1:size(mask,2)
            reg = reshape(region(m,n,:),1,[]);
            [found,ind]=ismember(reg(find(reg~=0)),all_regions);
            [containmstg,~] = ismember('mSTG',all_regions);
            if containmstg
                [ismstg,~]=ismember(reg(find(reg~=0)),{'mSTG'});
                if ismstg
                    if coord(2,m,n)>=(0.553486529*coord(3,m,n)-2.527049117)
                        [~,ind]=ismember('rSTG',all_regions);
                    else
                        [~,ind]=ismember('cSTG',all_regions);
                    end
                end
            end
            if found && mask(m,n)
                region_norm_value{ind} = [region_norm_value{ind},att_norm(:,:,m,n)];
                region_ref_value{ind} = [region_ref_value{ind},att_ref(:,:,m,n)];
            end
        end
    end
    for ind=1:length(region_value)
        region_norm_value{ind} = mean(region_norm_value{ind});
        region_ref_value{ind} = mean(region_ref_value{ind});
    end
    
    global noiselevel
    att_nt = att_nt./(att_norm_nt+1e-10).*noiselevel*ratio;
    if not(strcmp(ref,'causalpre'))
        att = att./(att_norm+1e-10).*att_ref*ratio;
%         att = att;
    else
        for m=1:size(mask,1)
            for n=1:size(mask,2)
                reg = reshape(region(m,n,:),1,[]);
                [found,ind]=ismember(reg(find(reg~=0)),{'rSTG','mSTG','cSTG'});
                if ~found && mask(m,n)
    %                     region_count{ind} = region_count(ind)+1;
                    att(:,:,m,n) = att(:,:,m,n)./(att_norm(:,:,m,n)+1e-10).*att_ref(:,:,m,n)*ratio;
%                     att(:,:,m,n) = att(:,:,m,n);
                else
                    att(:,:,m,n) = 0.6*att(:,:,m,n)./(att_norm(:,:,m,n)+1e-10).*att_ref(:,:,m,n)*ratio;
%                     att(:,:,m,n) = att(:,:,m,n)*0.6;
                end
            end
        end
    end

%     for m=1:size(mask,1)
%         for n=1:size(mask,2)
%             if not(strcmp(ref,'causalpre'))
%                 reg = reshape(region(m,n,:),1,[]);
%                 [found,indd]=ismember(reg(find(reg~=0)),all_regions);
%                 if found && mask(m,n)
%                     att(:,:,m,n) = att(:,:,m,n)./(region_norm_value{indd}+1e-10).*region_ref_value{indd}*ratio;
%                 end
%             else
%                 reg = reshape(region(m,n,:),1,[]);
%                 [found,indd]=ismember(reg(find(reg~=0)),all_regions);
%                 [foundd,ind]=ismember(reg(find(reg~=0)),{'rSTG','mSTG','cSTG'});
%                 if found && mask(m,n)
%                     if ~foundd && mask(m,n)
%         %                     region_count{ind} = region_count(ind)+1;
%         %                     att(:,:,m,n) = att(:,:,m,n)./(att_norm(:,:,m,n)+1e-10).*att_ref(:,:,m,n)*ratio;
%                         att(:,:,m,n) = att(:,:,m,n)./(region_norm_value{indd}+1e-10).*region_ref_value{indd}*ratio;
%                     else
%         %                     att(:,:,m,n) = 0.6*att(:,:,m,n)./(att_norm(:,:,m,n)+1e-10).*att_ref(:,:,m,n)*ratio;
%                         att(:,:,m,n) = 0.6*att(:,:,m,n)./(region_norm_value{indd}+1e-10).*region_ref_value{indd}*ratio;
%                     end
%                 end
%             end
%         end
%     end
    
        
        
%     if strcmp(ref,'anticausal')
%         att = att/1.5;
%     end
    win=hann(11); win = win./sum(win);
    for m=1:size(mask,1)
        for n=1:size(mask,2)
            if mask(m,n)
                for tau=tau0:tau1
                    switch ref
                        case 'causalpre'
                            region_value_nt{1} = [region_value_nt{1},[att_nt(rate*(tau-1)+1+t1*rate:rate*tau0+t2*rate,tau,m,n);zeros(rate*tau+t2*rate-(rate*tau0+t2*rate),1,1,1)]];
                            region_samples_nt{1} = region_samples_nt{1}+[ones(rate*tau0+t2*rate-(rate*(tau-1)+1+t1*rate)+1,1,1,1);zeros(rate*tau+t2*rate-(rate*tau0+t2*rate),1,1,1)];
                        otherwise
                            region_value_nt{1} = [region_value_nt{1},att_nt(rate*(tau-1)+1+t1*rate:rate*tau+t2*rate,tau,m,n)];
                    end
                end
            end
        end
    end
    if tticks(1)==0
        reflectdir = 'right';
    else
        reflectdir = 'left';
    end
    region_value_nt{1} = myconv(region_value_nt{1},reshape(win,size(win,1),1),reflectdir);
    switch  ref
        case 'causalpre'
            frames(:,end) = sum(region_value_nt{1},2)./(region_samples_nt{1}+1e-10);
    %                 frames(:,r) = frames(:,r)-min(frames(:,r));
            err(:,end) = sqrt(sum((region_value_nt{1}-frames(:,end)).^2,2)./(sum(region_samples_nt{1},2)+1e-10))./(sqrt(sum(region_samples_nt{1},2)*50/15)+1e-10)*1.96;
    %             frames(:,r) = mean(region_value{r},2);
    %             err(:,r) = std(region_value{r},1,2)/sqrt(size(region_value{r},2))/2;
    %             case 'causal'
    %                 frames(:,r) = sum(region_value{r},2)./(region_samples{r}+1e-10);
    %     %             frames(:,r) = frames(:,r)-min(frames(:,r));
    %                 err(:,r) = sqrt(sum((region_value{r}-frames(:,r)).^2,2)./(sum(region_samples{r},2)+1e-10))./(sqrt(sum(region_samples{r},2)*50)+1e-10);
    %     %             frames(:,r) = mean(region_value{r},2);
    %     %             err(:,r) = std(region_value{r},1,2)/sqrt(size(region_value{r},2))/2;
        otherwise
            frames(:,end) = mean(region_value_nt{1},2);
    %             frames(:,r) = frames(:,r)-min(frames(:,r));
            err(:,end) = std(region_value_nt{1},1,2)/sqrt(size(region_value_nt{1},2)*50/15)*1.96;
    end
    switch  ref
        case 'causalpre'
            frames(:,end) = frames(:,end).*reshape([1:size(frames(:,end),1)]*(-0.09396)+10,size(frames(:,end),1),1);
            frames(:,end) = frames(:,end)/max(frames(:,end))*noiselevel+0.03;
        case 'causal'
            frames(:,end) = frames(:,end).*reshape([1:size(frames(:,end),1)]*(-0.06396)+10,size(frames(:,end),1),1);
            frames(:,end) = frames(:,end)/max(frames(:,end))*noiselevel;
        otherwise
            frames(:,end) = frames(:,end).*reshape([1:size(frames(:,end),1)]*0.1396+10,size(frames(:,end),1),1);
            frames(:,end) = frames(:,end)/max(frames(:,end))*noiselevel+0.01;
    end
    
    
    for m=1:size(mask,1)
        for n=1:size(mask,2)
            reg = reshape(region(m,n,:),1,[]);
            [found,ind]=ismember(reg(find(reg~=0)),regions4barplot);
            [containmstg,~] = ismember('mSTG',regions4barplot);
            if containmstg
                [ismstg,~]=ismember(reg(find(reg~=0)),{'mSTG'});
                if ismstg
                    if coord(2,m,n)>=(0.553486529*coord(3,m,n)-2.527049117)
                        [~,ind]=ismember('rSTG',regions4barplot);
                    else
                        [~,ind]=ismember('cSTG',regions4barplot);
                    end
                end
            end
            if found && mask(m,n)
%                     region_count{ind} = region_count(ind)+1;
                for tau=tau0:tau1
%                 avg_data = avg_data+data_region(t,t+tau1:t+tau2);
                    switch ref
                        case 'causalpre'
                            region_value{ind} = [region_value{ind},[att(rate*(tau-1)+1+t1*rate:rate*tau0+t2*rate,tau,m,n);zeros(rate*tau+t2*rate-(rate*tau0+t2*rate),1,1,1)]];
                            region_samples{ind} = region_samples{ind}+[ones(rate*tau0+t2*rate-(rate*(tau-1)+1+t1*rate)+1,1,1,1);zeros(rate*tau+t2*rate-(rate*tau0+t2*rate),1,1,1)];
%                         case 'causal'
%                             start = max(rate*(tau-1)+1+t1*rate,4*rate+offset+1); endd = rate*tau+t2*rate;
%                             region_value{ind} = [region_value{ind},[zeros(start-(rate*(tau-1)+1+t1*rate),1,1,1);att(start:endd,tau,m,n)]];
%                             region_samples{ind} = region_samples{ind}+[zeros(start-(rate*(tau-1)+1+t1*rate),1,1,1);ones(endd-start+1,1,1,1)];
% %                             region_value{ind} = [region_value{ind},att(rate*(tau-1)+1+t1*rate:rate*tau+t2*rate,tau,m,n)];
                        otherwise
                            region_value{ind} = [region_value{ind},att(rate*(tau-1)+1+t1*rate:rate*tau+t2*rate,tau,m,n)];
                    end
%                     if not(strcmp(ref,'causalpre'))
%                         region_value{ind} = [region_value{ind},att(rate*(tau-1)+1+t1*rate:rate*tau+t2*rate,tau,m,n)];
%                     else
%                         region_value{ind} = [region_value{ind},[att(rate*(tau-1)+1+t1*rate:rate*tau0+t2*rate,tau,m,n);zeros(rate*tau+t2*rate-(rate*tau0+t2*rate),1,1,1)]];
%                         region_samples{ind} = region_samples{ind}+[ones(rate*tau0+t2*rate-(rate*(tau-1)+1+t1*rate)+1,1,1,1);zeros(rate*tau+t2*rate-(rate*tau0+t2*rate),1,1,1)];
%                     end
                end
            end
        end
    end
    disp('---------');
    global noiselevel;
    for ir=1:length(regions4barplot)
        xx = zeros(size(region_value{ir},2),1); 
        for sa=1:size(xx,1)
            xx(sa,1)=mean(region_value{ir}(region_value{ir}(:,sa)>0,sa)); 
        end; 
        [p,~,z]=signrank(xx-noiselevel);
        disp(p);disp(z);
    end
    
    if tticks(1)==0
        reflectdir = 'right';
    else
        reflectdir = 'left';
    end
    for r = 1:length(regions4barplot)
        region_value{r} = myconv(region_value{r},reshape(win,size(win,1),1),reflectdir);
    end
    for r = 1:length(regions4barplot)
        switch  ref
            case 'causalpre'
                frames(:,r) = sum(region_value{r},2)./(region_samples{r}+1e-10);
%                 frames(:,r) = frames(:,r)-min(frames(:,r));
                err(:,r) = sqrt(sum((region_value{r}-frames(:,r)).^2,2)./(sum(region_samples{r},2)+1e-10))./(sqrt(sum(region_samples{r},2)*50)+1e-10)*1.96;
    %             frames(:,r) = mean(region_value{r},2);
    %             err(:,r) = std(region_value{r},1,2)/sqrt(size(region_value{r},2))/2;
%             case 'causal'
%                 frames(:,r) = sum(region_value{r},2)./(region_samples{r}+1e-10);
%     %             frames(:,r) = frames(:,r)-min(frames(:,r));
%                 err(:,r) = sqrt(sum((region_value{r}-frames(:,r)).^2,2)./(sum(region_samples{r},2)+1e-10))./(sqrt(sum(region_samples{r},2)*50)+1e-10);
%     %             frames(:,r) = mean(region_value{r},2);
%     %             err(:,r) = std(region_value{r},1,2)/sqrt(size(region_value{r},2))/2;
            otherwise
                frames(:,r) = mean(region_value{r},2);
    %             frames(:,r) = frames(:,r)-min(frames(:,r));
                err(:,r) = std(region_value{r},1,2)/sqrt(size(region_value{r},2)*50)*1.96;
        end
    end
    clrmap =cbrewer('div', 'RdBu', 256,'PCHIP');
    clrp1 = clrmap(int32(length(clrmap)*0.75),:);
    clrp2 = clrmap(int32(length(clrmap)*1),:);
    clrn1 = clrmap(int32(length(clrmap)*0.25),:);
    clrn2 = clrmap(int32(1),:);
    if  tticks(1)==0
        thres = 1.5;% 1.2;%1.5;
    else
        thres= 1.5;%0.8;%1.5;
    end
    x = [];
    g = {};
    x_mean=[];
    g_mean=[];
    max_cell={};
    maxvalue_cell={};
    name_cell={};
    max_mean=[];
    max_meanmean=[];
    rr=1;
%     maxplotregion_causal={'inferiorprecentral','superiorprecentral','postcentral','parsopercularis'};
%     maxplotregion_anticausal={'cSTG','inferiorprecentral','superiorprecentral','postcentral'};
%     maxplotregion_causal={'inferiorprecentral','superiorprecentral','postcentral','parsopercularis','parstriangularis'};
%     maxplotregion_anticausal={'cSTG','mMTG','cMTG','rMTG','inferiorprecentral','superiorprecentral','postcentral','supramarginal','parsopercularis','parstriangularis','rostralmiddlefrontal','caudalmiddlefrontal'};
    maxplotregion_causal=regions4barplot;
    maxplotregion_anticausal=regions4barplot;
    for r=1:length(regions4barplot)
        if  tticks(1)==0
            [found,~]=ismember(regions4barplot{r},maxplotregion_anticausal);
        else
            [found,~]=ismember(regions4barplot{r},maxplotregion_causal);
        end
        if ~found
            continue;
        end
%         [maxvalue,argmax]=max(region_value{r}(:,max(region_value{r},[],1)>thres));
%         [maxvalue,argmax]=max(frames(:,r));
        win=hann(7); win = win./sum(win); smoothcurve = myconv(squeeze(frames(:,r)),reshape(win,size(win,1),1),reflectdir);
        [maxvalue,argmax]=max(smoothcurve);
        
        argmax = (argmax-1)*1000/125;
        if tticks(1)~=0
            argmax=argmax-size(region_value{r},1)*1000/125;
        end
%         if size(argmax,2)<3
%             continue;
%         end
        max_cell{rr}=argmax;
        maxvalue_cell{rr}=maxvalue;
        name_cell{rr}=repmat({name4plot{r}},1,size(argmax,2));
        max_mean=[max_mean,median(argmax)];
        max_meanmean=[max_meanmean,sum(argmax.*maxvalue)/sum(maxvalue(:))];
        rr=rr+1;
    end
    [sorted_max,index] = sort(max_meanmean);
%     [sorted_max,index] = sort(max_mean);
    for rr=1:length(max_cell)
        x = [x,max_cell{index(rr)}];
        g = [g,name_cell{index(rr)}];
        x_mean = [x_mean,sum(max_cell{index(rr)}.*maxvalue_cell{index(rr)})/sum(maxvalue_cell{index(rr)})];
        g_mean = [g_mean,{name_cell{index(rr)}{1}}];
    end
    figure();boxplot(x,g);xtickangle(45)
%     if tticks(1)==0
%         figure();b=violin(max_cell,'facecolor',clrn1,'edgecolor',clrn2,'mc',[],'medc',[],'facealpha',1);set(gca, 'XTickLabel', g_mean);xtickangle(45)
%     else
%         figure();b=violin(max_cell,'facecolor',clrp1,'edgecolor',clrp2,'mc',[],'medc',[],'facealpha',1);set(gca, 'XTickLabel', g_mean);xtickangle(45)
%     end
    
    figure();b=bar(x_mean);set(gca, 'XTickLabel', g_mean);xtickangle(45)
    if tticks(1)==0
        b.FaceColor=clrn1;b.EdgeColor=clrn2; b.EdgeAlpha=1; b.FaceAlpha=1;
    else
        b.FaceColor=clrp1;b.EdgeColor=clrp2; b.EdgeAlpha=1; b.FaceAlpha=1;
    end
    ax=gca;ax.FontSize=20;  set(gca,'linewidth',2)
    
    
    frames_brain = zeros(t2*rate-t1*rate+rate,tau1-tau0+1,size(att,3),size(att,4));
    mask_reshape = reshape(mask,1,1,size(att,3),size(att,4));
    masked_att = mask_reshape.*att;
    summ=0;
    for tau=tau0:tau1
        switch ref
            case 'causalpre'
                frames_brain(:,tau-tau0+1,:,:) = [masked_att(rate*(tau-1)+1+t1*rate:rate*tau0+t2*rate,tau,:,:);zeros(rate*tau+t2*rate-(rate*tau0+t2*rate),1,size(masked_att,3),size(masked_att,4))];
                summ = summ+[ones(rate*tau0+t2*rate-(rate*(tau-1)+1+t1*rate)+1,1,size(masked_att,3),size(masked_att,4));zeros(rate*tau+t2*rate-(rate*tau0+t2*rate),1,size(masked_att,3),size(masked_att,4))];
%             case 'causal'
%                 start = max(rate*(tau-1)+1+t1*rate,4*rate+offset+1); endd = rate*tau+t2*rate;
%                 frames_brain(:,tau-tau0+1,:,:) = [zeros(start-(rate*(tau-1)+1+t1*rate),1,size(masked_att,3),size(masked_att,4));masked_att(start:endd,tau,:,:)];
%                 summ = summ+[zeros(start-(rate*(tau-1)+1+t1*rate),1,size(masked_att,3),size(masked_att,4));ones(endd-start+1,1,size(masked_att,3),size(masked_att,4))];
            otherwise
                frames_brain(:,tau-tau0+1,:,:) = masked_att(rate*(tau-1)+1+t1*rate:rate*tau+t2*rate,tau,:,:);
                summ = summ+ones(rate*tau+t2*rate-(rate*(tau-1)+t1*rate),1,size(masked_att,3),size(masked_att,4));
        end
%         frames_brain(:,tau-tau0+1,:,:) = masked_att(rate*(tau-1)+1+t1*rate:rate*tau+t2*rate,tau,:,:);
    end
    
%     [color_map_]=getPyPlot_cMap('turbo', 256);
    [color_map_] = jet(256);
%     color_map_ = magma(256); color_map_=color_map_(end:-1:1,:);
    color_map_(1:2,:) = ones(size(color_map_(1:2,:)));
    frames_brain_time = mean(frames_brain,2);
    if  tticks(1)==0
        thres = 1.5;
    else
        thres=1;
    end
    maxv = 40;
    mask_peakvalue = max(frames_brain_time,[],1)>thres;
    mask_peakvalue = squeeze(mask_peakvalue).*mask;
    [~,peak_time] = max(frames_brain_time);
    if tticks(1)~=0
        peak_time = size(frames_brain,1)-peak_time;
    end
%     maxv = max(peak_time(find(reshape(mask_peakvalue,1,1,size(mask_peakvalue,1),size(mask_peakvalue,2)))));
    
    [color_map_] = jet(256);
%     color_map_ = magma(256); color_map_=color_map_(end:-1:1,:);
    color_map_(1:int64(15/maxv*256),:) = ones(size(color_map_(1:int64(15/maxv*256),:)));
%     VisualAtt({peak_time}, {mask_peakvalue}, {region}, {coord}, brain,color_map_,maxv,true,false,false,false,true,[0,1],7);
    
    frames_brain = permute(frames_brain,[2,1,3,4]);
    frames_brain = reshape(frames_brain,size(frames_brain,1)*rate,size(frames_brain,2)/rate,size(frames_brain,3),size(frames_brain,4));
%     frames_brain = mean(frames_brain,1)/1.5;
    summ = permute(summ,[2,1,3,4]);
    summ = reshape(summ,size(summ,1)*rate,size(summ,2)/rate,size(summ,3),size(summ,4));
    frames_brain = sum(frames_brain,1)./sum(summ,1)/1.5;
    
%     frames_brain = max(frames_brain-min(frames_brain(:)),0);
%     frames_brain = permute(frames_brain,[2,1,3,4]);
%     frames_brain = sum(frames_brain,1)./(8*reshape([1,2,3,4,5,5,5,5],[1,8,1,1]));
    [frames_brain] = VisualAtt({frames_brain}, {mask}, {region}, {coord}, brain,colormap,1,true,false,false,false,true,[0,0.7],4);
%     for r = 1:length(regions4barplot)
%         elecs = 0;
%         data_ofregion=0;
%         for sub=1:length(Att1_cell)
% %             data = Att1_cell{sub};
%             data = abs(Att1_cell{sub});
%             data = reshape(data,[col,N_tau,size(data,3),size(data,4)]);
%             target = isregion(regions,regions4barplot{r});
%             elecs = elecs+sum(target(:));
%             targetrep = repmat(reshape(target,[1,1,15,15]),col,N_tau,1,1);
%             data_region = data.*targetrep;
%             data_region = sum(sum(data_region,4),3);
%             avg_data = 0;
%             for t=t0:t1
% %                 avg_data = avg_data+data_region(t,t+tau1:t+tau2);
%                 avg_data = avg_data+data_region(t,8*(t-1)+1+tau1*8:8*t+tau2*8);
%             end
%             denorm = ones(size(avg_data))*(t1-t0+1);
% %             denorm(1:8)=1; denorm(9:16)=1; denorm(17:24)=2; denorm(25:32)=3;
%             data_ofregion = data_ofregion+avg_data./denorm;
%         end
%         data_ofregion = data_ofregion/elecs;
%         frames(:,r) = data_ofregion;
%     end
end

function [frames,prob,prob0] = VisualAtt(Att, Mask, Region, coord, brain,color_map,maxv,alpha,div,annotation,bar_plot,gaussianplot,scale,gausskernel)
    if ~exist('annotation','var')
        annotation = [];
    end
    if ~exist('gausskernel','var')
        gausskernel = 10;
    end
    if ~exist('gaussianplot','var')
        gaussianplot = true;
    end
    if ~exist('bar_plot','var')
        bar_plot = false;
    end
    
    if nargin<5
        color_map = 'hot';
        color_map = eval([color_map,'(256)']);
    end
    annot = brain.annot;
    regions4barplot = {'rSTG','cSTG','rMTG','mMTG','cMTG','inferiorprecentral','superiorprecentral','postcentral','caudalmiddlefrontal','rostralmiddlefrontal','parstriangularis','parsopercularis','parsorbitalis','supramarginal'};
    att=[];
    mask=[];
    cod = [];
    region = [];
    for sub = 1:length(Att)
        att = cat(3,att,Att{sub});
        mask = cat(1,mask,Mask{sub});
        cod = cat(2,cod,coord{sub});
        region = cat(1,region,Region{sub});
    end
    load('annot_regions.mat');
    if bar_plot
        annot_regions = regions4barplot;
    end
%     mask = Mask{sub};
    if length(size(att))==4
        sizes = size(att);
        Att_avg = reshape(mean(att,1),sizes(2:end));
    else 
        Att_avg = att;
    end
    if nargin<6
        maxv = max(Att_avg(:));
    end
%     Att_avg = (Att_avg-min(Att_avg(:)))/(max(Att_avg(:))-min(Att_avg(:)));
    Att_avg = (Att_avg)/maxv;
%     mask(Att_avg<0.2)=0;
    % Att_avg(1:6,:,:) = Att_avg(1:6,:,:)*0;
%     if length(size(att))==3
%         time_step = size(Att_avg,1);
%     else
%         time_step = 1;
%     end
    time_step = size(Att_avg,1);
    if gaussianplot
        frames={};
    else
        frames = [];
    end
    for i=1:time_step
%         if length(size(att))==3
%             Att_frame = squeeze(Att_avg(i,:,:));
%         else
%             Att_frame = Att_avg;
%         end
        Att_frame = squeeze(Att_avg(i,:,:));
        region_value = cell(length(annot_regions),1);
%         region_count = cell(length(annot_regions),1);
%         region_value = {};
        for m=1:size(mask,1)
            for n=1:size(mask,2)
                reg = reshape(region(m,n,:),1,[]);
                [found,ind]=ismember(reg(find(reg~=0)),annot_regions);
                [ismstg,~]=ismember(reg(find(reg~=0)),{'mSTG'});
                if ismstg
                    if cod(2,m,n)>=0.553486529*cod(3,m,n)-2.527049117
                        [~,ind]=ismember('rSTG',annot_regions);
                    else
                        [~,ind]=ismember('cSTG',annot_regions);
                    end
                end
                if found && mask(m,n)
%                     region_count{ind} = region_count(ind)+1;
                    region_value{ind} = [region_value{ind},Att_frame(m,n)];
                end
            end
        end
        for m=1:length(region_value)
            region_value_mean(m) = mean(region_value{m});
            region_value_stderr(m) = std(region_value{m})/sqrt(length(region_value{m}));
        end
        disp('------');
        if bar_plot
            for ii=1:length(region_value)
                p=zeros(1000,1);z=zeros(1000,1);
                n=15;
%                 n=1;
                for iii=1:1000;[p(iii),~,zval]=signrank(randn(size(region_value{ii},2)*n,1)*std(region_value{ii})+mean(region_value{ii}));z(iii)=zval.zval;end;disp(mean(p));disp(mean(z));
%                 [p,~,z]=signrank(region_value{ii});
%                 disp(p);disp(z);
            end;
            fig = barplot(region_value,color_map,regions4barplot); ax=gca;ax.FontSize=13;  set(gca,'linewidth',2)
            frames = fig;
            prob = frames;
            prob0 = frames;
%             close(fig);
        else
            region_value = region_value_mean;
%             if div
%                 region_value = (region_value-0.5)*2+0.5;
%             else
%                 region_value = region_value*2;
%             end
%             ind_active = find(mask==1& double(Att_frame<0.16));
            ind_active = find(mask==1);
    %         color_map = 'hot';

            if div
                color_ind = value2colorind(Att_frame(ind_active), 'hot',[0,1]);
                color_ind_region = value2colorind(region_value, 'hot',[0,1]);
            else
                color_ind = value2colorind(Att_frame(ind_active), 'hot',[0,1]);
                color_ind_region = value2colorind(region_value, 'hot',[0,1]);
            end
            if annotation
                regioncolor = color_map(color_ind_region,:);
                regioncolor(isnan(region_value),:) = repmat([0.92 0.92 0.92],length(find(isnan(region_value))),1);
            end
    %             colors = {color_ind, color_map};
            if size(color_map,1)<100
                if size(color_map,1)==1
                    colors = {color_map};
                else
                    colors = {color_map(sub,:)};
                end
            else
                colors = {color_ind, color_map};
            end
            elecname = 1:length(ind_active);
            elec = zeros(length(ind_active),3);
            Temp = squeeze(cod(1,:,:));
            elec(:,1) = Temp(ind_active);
            Temp = squeeze(cod(2,:,:));
            elec(:,2) = Temp(ind_active);
            Temp = squeeze(cod(3,:,:));
            elec(:,3) = Temp(ind_active);
    %         fig = nyu_plot_whitebackground(brain.surf_brain,...
    %                                        brain.sph,...
    %                                        elec, elecname,...
    %                                        colors, 0, brain.annot, 2, 1);

            if alpha
                if div
                    point_alpha = 2*abs(max(min(Att_frame(ind_active),1),0)-0.5);
    %                 point_alpha = max(min(abs(Att_frame(ind_active)),1),0);
                else
                    point_alpha = max(min(Att_frame(ind_active),1),0);
                end
                point_alpha = double(point_alpha>0.1);
            else
                point_alpha = double(abs(Att_frame(ind_active)-0.5)>0.07);%ones(size(Att_frame(ind_active)));
    %             point_alpha = ones(size(Att_frame(ind_active)));
            end
            if annotation
                braincolors = regioncolor;
                point_alpha = point_alpha*0;
            else
                braincolors = [];
                
                brain.annot = [];
            end
            regionreshape = reshape(region,size(region,1)*size(region,2),size(region,3));
            validregion = regionreshape(ind_active,:);
            if gaussianplot
                load('S.mat')
                Lcrtx = load('ch2_template_lh_pial.mat');
                Lcrtx.tri = Lcrtx.faces;Lcrtx.vert = Lcrtx.coords;
                if div
                    data = (Att_frame(ind_active)-0.5)*2;
                else
                    data = Att_frame(ind_active);
                end
                [fig,prob_all]=ctmr_gauss_plot_edited(Lcrtx,elec,ones(size(elec,1),1),S.cax,0,S.cm,S.gsp*gausskernel,[],annot,validregion); 
                prob0 = prob_all;
                close(fig)
                global noiselevel;
                if div
                    S.cax=[-0.7,0.7];
                else
                    S.cax = [0,0.7];
                end
                if exist('scale','var')
                    S.cax = scale;
%                     S.cax(1) = noiselevel;
                end
%                 S.cax = [0 0.5];
%                 global noiselevel;
                color_map(1:int64(size(color_map,1)*noiselevel/S.cax(2)),:)=1;
                [fig,prob]=ctmr_gauss_plot_edited(Lcrtx,elec,data,S.cax,0,color_map,S.gsp*gausskernel,prob_all,annot,validregion); 
                viewangle = 'l';   
                pltshowplanes = 0;
                cameratoolbar('setmode',''); 
                litebrain(viewangle,.8); 
                hold on; colormap(gca,color_map);  c = colorbar;
                axis off;
                frames{i}=fig;
            else
                fig = nyu_plot_whitebackground(brain.surf_brain,...
                                               brain.sph,...
                                               elec, elecname,...
                                               colors, 0, brain.annot, 2, 1, 1, point_alpha,braincolors);
                frame = getframe(fig);
                frames(i,:,:,:) = frame.cdata;
                prob = frames;
                prob0 = frames;
%                 close(fig);
            end
    %         fig = nyu_plot_whitebackground(brain.surf_brain,...
    %                                        brain.sph,...
    %                                        elec, elecname,...
    %                                        colors, 0, [], 2, 1, 1);

        end
    end
%     figure();
%     if size(frames,1)==1
%         imagesc(uint8(squeeze(frames(1,:,:,:))));
%         axis off;
%     else
%         haxis = tight_subplot(6,6,[0,0],0,0);
%         for i=1:36
%             axes(haxis(i));
%             imagesc(uint8(squeeze(frames(i,:,:,:))));
%             axis off;
%         end
%     end
    
end

function [frames,prob,prob0] = VisualAtt4barplot(Att, Mask, Region, coord, brain,color_map,maxv,alpha,div,annotation,bar_plot,gaussianplot)
    if ~exist('annotation','var')
        annotation = [];
    end
    if ~exist('gaussianplot','var')
        gaussianplot = true;
    end
    if ~exist('bar_plot','var')
        bar_plot = false;
    end
    if nargin<5
        color_map = 'hot';
        color_map = eval([color_map,'(256)']);
    end
    annot = brain.annot;
    regions4barplot = {'rSTG','cSTG','rMTG','mMTG','cMTG','inferiorprecentral','superiorprecentral','postcentral','caudalmiddlefrontal','rostralmiddlefrontal','parstriangularis','parsopercularis','supramarginal'};
    att=[];
    mask=[];
    cod = [];
    region = [];
    att2=[];
    for sub = 1:length(Mask)
        att = cat(3,att,Att{1}{sub});
        att2 = cat(3,att2,Att{2}{sub});
        mask = cat(1,mask,Mask{sub});
        cod = cat(2,cod,coord{sub});
        region = cat(1,region,Region{sub});
    end
    load('annot_regions.mat');
    if bar_plot
        annot_regions = regions4barplot;
    end
%     mask = Mask{sub};
    if length(size(att))==4
        sizes = size(att);
        Att_avg = reshape(mean(att,1),sizes(2:end));
        Att_avg2 = reshape(mean(att2,1),sizes(2:end));
    else 
        Att_avg2 = att2;
    end
    if nargin<6
        maxv = max(Att_avg(:));
        maxv2 = max(Att_avg2(:));
    end
%     Att_avg = (Att_avg-min(Att_avg(:)))/(max(Att_avg(:))-min(Att_avg(:)));
    Att_avg = (Att_avg)/maxv;
    Att_avg2 = (Att_avg2)/maxv;
%     mask(Att_avg<0.2)=0;
    % Att_avg(1:6,:,:) = Att_avg(1:6,:,:)*0;
%     if length(size(att))==3
%         time_step = size(Att_avg,1);
%     else
%         time_step = 1;
%     end
    time_step = size(Att_avg,1);
    if gaussianplot
        frames={};
    else
        frames = [];
    end
    for i=1:time_step
%         if length(size(att))==3
%             Att_frame = squeeze(Att_avg(i,:,:));
%         else
%             Att_frame = Att_avg;
%         end
        Att_frame = squeeze(Att_avg(i,:,:));
        Att_frame2 = squeeze(Att_avg2(i,:,:));
        region_value = cell(length(annot_regions),1);
        region_value2 = cell(length(annot_regions),1);
        elec_value = [];
        elec_value2=[];
%         region_count = cell(length(annot_regions),1);
%         region_value = {};
        for m=1:size(mask,1)
            for n=1:size(mask,2)
                reg = reshape(region(m,n,:),1,[]);
                [found,ind]=ismember(reg(find(reg~=0)),annot_regions);
                [ismstg,~]=ismember(reg(find(reg~=0)),{'mstg'});
                if ismstg
                    if cod(2,m,n)>=0.553486529*cod(3,m,n)-2.527049117
                        [~,ind]=ismember('rstg',annot_regions);
                    else
                        [~,ind]=ismember('cstg',annot_regions);
                    end
                end
                if found && mask(m,n)
%                     region_count{ind} = region_count(ind)+1;
                    region_value{ind} = [region_value{ind},Att_frame(m,n)];
                    region_value2{ind} = [region_value2{ind},Att_frame2(m,n)];
                    elec_value = [elec_value,Att_frame(m,n)];
                    elec_value2 = [elec_value2,Att_frame2(m,n)];
                end
            end
        end
        for ii=1:length(region_value); disp(signrank(region_value{ii},region_value2{ii}));end;
        fig = barplot2(region_value,region_value2,color_map,regions4barplot);
%         fig = barplot3(elec_value,elec_value2,color_map,regions4barplot);
        frames = fig;
        prob = frames;
        prob0 = frames;
%             close(fig);
    end
%     figure();
%     if size(frames,1)==1
%         imagesc(uint8(squeeze(frames(1,:,:,:))));
%         axis off;
%     else
%         haxis = tight_subplot(6,6,[0,0],0,0);
%         for i=1:36
%             axes(haxis(i));
%             imagesc(uint8(squeeze(frames(i,:,:,:))));
%             axis off;
%         end
%     end
    
end

function [frames,prob,prob0] = VisualAtt4normcontr(Att,Att_c,Att_a, Mask, Region, coord, brain,color_map,maxv,alpha,div,annotation)
   if ~exist('annotation','var')
        annotation = [];
    end
    if ~exist('gausskernel','var')
        gausskernel = 10;
    end
    if ~exist('gaussianplot','var')
        gaussianplot = true;
    end
    if ~exist('bar_plot','var')
        bar_plot = false;
    end
    
    if nargin<5
        color_map = 'hot';
        color_map = eval([color_map,'(256)']);
    end
    annot = brain.annot;
    regions4barplot = {'rSTG','cSTG','rMTG','mMTG','cMTG','inferiorprecentral','superiorprecentral','postcentral','caudalmiddlefrontal','rostralmiddlefrontal','parstriangularis','parsopercularis','supramarginal'};
    att=[];
    att_c=[];
    att_a=[];
    mask=[];
    cod = [];
    region = [];
    for sub = 1:length(Att)
        att = cat(3,att,Att{sub});
        att_c = cat(3,att_c,Att_c{sub});
        att_a = cat(3,att_a,Att_a{sub});
        mask = cat(1,mask,Mask{sub});
        cod = cat(2,cod,coord{sub});
        region = cat(1,region,Region{sub});
    end
    the_power=0.43;
    att_c = (max(att_c,0)).^the_power;
    att_a = (max(att_a,0)).^the_power;
    global noiselevel
    mask_largelength = squeeze(double((att_c>noiselevel) | (att_a>noiselevel)));
    att = att/2+0.5;
    maxv = 1;
    load('annot_regions.mat');
    if bar_plot
        annot_regions = regions4barplot;
    end
%     mask = Mask{sub};
    if length(size(att))==4
        sizes = size(att);
        Att_avg = reshape(mean(att,1),sizes(2:end));
    else 
        Att_avg = att;
    end
    if nargin<6
        maxv = max(Att_avg(:));
    end
%     Att_avg = (Att_avg-min(Att_avg(:)))/(max(Att_avg(:))-min(Att_avg(:)));
    Att_avg = (Att_avg)/maxv;
%     mask(Att_avg<0.2)=0;
    % Att_avg(1:6,:,:) = Att_avg(1:6,:,:)*0;
%     if length(size(att))==3
%         time_step = size(Att_avg,1);
%     else
%         time_step = 1;
%     end
    time_step = size(Att_avg,1);
    if gaussianplot
        frames={};
    else
        frames = [];
    end
    for i=1:time_step
%         if length(size(att))==3
%             Att_frame = squeeze(Att_avg(i,:,:));
%         else
%             Att_frame = Att_avg;
%         end
        Att_frame = squeeze(Att_avg(i,:,:));
        region_value = cell(length(annot_regions),1);
%         region_count = cell(length(annot_regions),1);
%         region_value = {};
        for m=1:size(mask,1)
            for n=1:size(mask,2)
                reg = reshape(region(m,n,:),1,[]);
                [found,ind]=ismember(reg(find(reg~=0)),annot_regions);
                [ismstg,~]=ismember(reg(find(reg~=0)),{'mSTG'});
                if ismstg
                    if cod(2,m,n)>=0.553486529*cod(3,m,n)-2.527049117
                        [~,ind]=ismember('rSTG',annot_regions);
                    else
                        [~,ind]=ismember('cSTG',annot_regions);
                    end
                end
                if found && mask(m,n)
%                     region_count{ind} = region_count(ind)+1;
                    region_value{ind} = [region_value{ind},Att_frame(m,n)];
                end
            end
        end
        for m=1:length(region_value)
            region_value_mean(m) = mean(region_value{m});
            region_value_stderr(m) = std(region_value{m})/sqrt(length(region_value{m}));
        end
        region_value = region_value_mean;
%             if div
%                 region_value = (region_value-0.5)*2+0.5;
%             else
%                 region_value = region_value*2;
%             end
%         ind_active = find(mask==1 & abs(Att_frame*2-1)>0.9);
        ind_active = find(mask==1);
%         color_map = 'hot';

        if div
            color_ind = value2colorind(Att_frame(ind_active), 'hot',[0,1]);
            color_ind_region = value2colorind(region_value, 'hot',[0,1]);
        else
            color_ind = value2colorind(Att_frame(ind_active), 'hot',[0,1]);
            color_ind_region = value2colorind(region_value, 'hot',[0,1]);
        end
        if annotation
            regioncolor = color_map(color_ind_region,:);
            regioncolor(isnan(region_value),:) = repmat([0.92 0.92 0.92],length(find(isnan(region_value))),1);
        end
%             colors = {color_ind, color_map};
        if size(color_map,1)<100
            if size(color_map,1)==1
                colors = {color_map};
            else
                colors = {color_map(sub,:)};
            end
        else
            color_ind(color_ind>128) = 217;
            color_ind(color_ind<=128) = 38;
            colors = {color_ind, color_map};
        end
        cmap = zeros(256,3);
        for cc=1:3;cmap(1:128,cc) = linspace(color_map(217,cc),1,128);end
        for cc=1:3;cmap(129:256,cc) = linspace(1,color_map(38,cc),128);end
        elecname = 1:length(ind_active);
        elec = zeros(length(ind_active),3);
        Temp = squeeze(cod(1,:,:));
        elec(:,1) = Temp(ind_active);
        Temp = squeeze(cod(2,:,:));
        elec(:,2) = Temp(ind_active);
        Temp = squeeze(cod(3,:,:));
        elec(:,3) = Temp(ind_active);
%         fig = nyu_plot_whitebackground(brain.surf_brain,...
%                                        brain.sph,...
%                                        elec, elecname,...
%                                        colors, 0, brain.annot, 2, 1);

        if alpha
            if div
                point_alpha = 2*abs(max(min(Att_frame(ind_active),1),0)-0.5);
%                 point_alpha = max(min(abs(Att_frame(ind_active)),1),0);
            else
                point_alpha = max(min(Att_frame(ind_active),1),0);
            end
            point_alpha = double(point_alpha>0.1);
        else
%             point_alpha = ones(size(Att_frame(ind_active))).*mask_largelength(ind_active);
%             point_alpha = double(abs(Att_frame(ind_active)*2-1)>0.9).*mask_largelength(ind_active);
            point_alpha = double(abs(Att_frame(ind_active)*2-1)).*mask_largelength(ind_active);
            

        end
        if annotation
            braincolors = regioncolor;
            point_alpha = point_alpha*0;
        else
            braincolors = [];

            brain.annot = [];
        end
        regionreshape = reshape(region,size(region,1)*size(region,2),size(region,3));
        validregion = regionreshape(ind_active,:);
        

        fig = nyu_plot_whitebackground(brain.surf_brain,...
                                       brain.sph,...
                                       elec, elecname,...
                                       colors, 0, brain.annot, 2, 1, 1, point_alpha,braincolors);
        frame = getframe(fig);
        frames(i,:,:,:) = frame.cdata;
        prob = frames;
        prob0 = frames;
        close(fig);

%         fig = nyu_plot_whitebackground(brain.surf_brain,...
%                                        brain.sph,...
%                                        elec, elecname,...
%                                        colors, 0, [], 2, 1, 1);
    end
%     figure();
%     if size(frames,1)==1
%         imagesc(uint8(squeeze(frames(1,:,:,:))));
%         axis off;
%     else
%         haxis = tight_subplot(6,6,[0,0],0,0);
%         for i=1:36
%             axes(haxis(i));
%             imagesc(uint8(squeeze(frames(i,:,:,:))));
%             axis off;
%         end
%     end
    
end

function [rgb_inds] = value2colorind(value, c_map, c_range)
    if ~exist('c_map','var') || isempty(c_map)
        c_map = 'hot';
    end
    if ~exist('c_range','var') || isempty(c_range)
        c_range = [mean(value(:)), max(value(:))];
    end
    cmap = eval([c_map,'(256)']);
    % clip the out of range values
    value(value<c_range(1)) = c_range(1);
    value(value>c_range(2)) = c_range(2);
    % normalize the values
    rgb_inds = round(((value-c_range(1))/(c_range(2)-c_range(1)))*255)+1;
    rgb_inds(isnan(rgb_inds)) = 1;
end
function [elec_all, elecname_all] = multi_subj_elecs()
    root_dir = pwd;
    Subjs = {'NY717','NY742','NY749'};
    elec_all = [];
    elecname_all = [];
    for i=1:length(Subjs)
        Subj = Subjs{i};
        paths = get_path(root_dir,Subj);
        channel_info_all = get_channel_info(paths);
        plot_data_all = get_plot_data(paths, 'mni', 'lh');
        elec_range = 1:128;
        elec = table2array(channel_info_all(elec_range,plot_data_all.coord_ind));
        elecname = table2array(channel_info_all(elec_range,11));
        elec_all = [elec_all;elec];
        elecname_all = [elecname_all; elecname];
    end
end
function [plot_data_all] = get_plot_data(paths, visual_mode, sph)
    if strcmp(lower(visual_mode), 'mni')
        if strcmp(lower(sph), 'lh')
            plot_data_all.surf_brain = load(paths.MNI_lh);
            plot_data_all.sph = 'lh';
            plot_data_all.annot = paths.MNI_lh_annot;
            plot_data_all.coord_ind = 2:4;
        elseif strcmp(lower(sph), 'rh')
            plot_data_all.surf_brain = load(paths.MNI_rh);
            plot_data_all.sph = 'rh';
            plot_data_all.annot = paths.MNI_rh_annot;
            plot_data_all.coord_ind = 2:4;
        else
            error('sph mode not supported! use lh or rh');
        end
    elseif strcmp(lower(visual_mode), 'subj')
        if strcmp(lower(sph), 'lh')
            plot_data_all.surf_brain = load(paths.Subj_lh);
            plot_data_all.sph = 'lh';
            plot_data_all.annot = paths.Subj_lh_annot;
            plot_data_all.coord_ind = 6:8;
        elseif strcmp(lower(sph), 'rh')
            plot_data_all.surf_brain = load(paths.Subj_rh);
            plot_data_all.sph = 'rh';
            plot_data_all.annot = paths.Subj_rh_annot;
            plot_data_all.coord_ind = 6:8;
        else
            error('sph mode not supported! use lh or rh');
        end
    else
        error('visual_mode not supported! use mni or subj');
    end
end
function [channel_info_all] = get_channel_info(paths)
    % load coordinate files
    coordinates = readtable(paths.coordinate_file);
    ind = [1:8,10,11];
    coordinates = coordinates(:,ind);
    coordinates.channel= [1:size(coordinates,1)]';
    Subj.sebject = repmat(paths.Subj,size(coordinates,1),1);
    Subj = struct2table(Subj);
%     Subj = table('Size',[size(coordinates,1),1],...
%                  'VariableTypes', {'string'},...
%                  'VariableNames', {'subject'});
%     Subj.subject = paths.Subj;
    channel_info_all = [coordinates, Subj];
end
function [paths] = get_path(root_dir, Subj)
    addpath(genpath(root_dir));
    paths.root_dir = root_dir;
    paths.Subj = Subj;
    paths.Subj_path = [root_dir, filesep, 'Data', filesep, Subj];
    paths.coordinate_file = [paths.Subj_path, filesep, 'coordinates.csv'];
    paths.MNI_lh = [paths.root_dir, filesep, 'Data', filesep, 'MNI', filesep, 'ch2_template_lh_pial_120519.mat'];
    paths.MNI_rh = [paths.root_dir, filesep, 'Data', filesep, 'MNI', filesep, 'ch2_template_rh_pial_120519.mat'];
    paths.MNI_lh_annot = [paths.root_dir, filesep, 'Data', filesep, 'MNI', filesep, 'ch2_template.lh.aparc.split_STG_MTG.annot'];
    paths.MNI_rh_annot = [paths.root_dir, filesep, 'Data', filesep, 'MNI', filesep, 'ch2_template.rh.aparc.split_STG_MTG.annot'];
    lh_pial_file = dir([paths.Subj_path,filesep,'*lh_pial_surf.mat']);
    lh_annot_file = dir([paths.Subj_path,filesep,'lh*_STG_MTG.annot']);
    rh_pial_file = dir([paths.Subj_path,filesep,'*rh_pial_surf.mat']);
    rh_annot_file = dir([paths.Subj_path,filesep,'rh*_STG_MTG.annot']);
    if (length(lh_pial_file)~=0)
        paths.Subj_lh = [paths.Subj_path, filesep, lh_pial_file(1).name];
    else
        paths.Subj_lh = [];
    end
    if (length(rh_pial_file)~=0)
        paths.Subj_rh = [paths.Subj_path, filesep, rh_pial_file(1).name];
    else
        paths.Subj_rh = [];
    end
    if (length(lh_annot_file)~=0)
        paths.Subj_lh_annot = [paths.Subj_path, filesep, lh_annot_file(1).name];
    else
        paths.Subj_lh_annot = [];
    end
    if (length(rh_annot_file)~=0)
        paths.Subj_rh_annot = [paths.Subj_path, filesep, rh_annot_file(1).name];
    else
        paths.Subj_rh_annot = [];
    end
end
function coef = find_normcoef(Att2,Att1,mask,max_anker,count)
sorted = sort(Att2(:));
num = length(sorted);
x0 = mean(sorted(1:count));
if ~max_anker
    sorted = sort(Att2(repmat(reshape(mask,[1,1,size(mask,1),size(mask,2)]),size(Att2,1),size(Att2,2),1,1)==1));
    x1 = mean(sorted(1:count));
%     x1 = mean(sorted(int32(num/2-count/2):int32(num/2+count/2)));
%     x1 = mean(sorted);
else
    x1 = mean(sorted(end));
end
sorted = sort(Att1(:));
y0 = mean(sorted(1:count));
if ~max_anker
    sorted = sort(Att1(repmat(reshape(mask,[1,1,size(mask,1),size(mask,2)]),size(Att1,1),size(Att1,2),1,1)==1));
    y1 = mean(sorted(1:count));
%     y1 = mean(sorted(int32(num/2-count/2):int32(num/2+count/2)));
%     y1 = mean(sorted);
else
    y1 = mean(sorted(end));
end
coef = [x0,1;x1,1]\[y0;y1];
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
function [mask] = isaud(regions)
    sz = size(regions);
    mask = zeros(sz(1:length(sz)-1));
    for i =1:size(regions,1)
        for j=1:size(regions,2)
            if contains(reshape(regions(i,j,:),1,size(regions,3)),'mSTG') || contains(reshape(regions(i,j,:),1,size(regions,3)),'cSTG')
                mask(i,j) = true;
            end
        end
    end
end
% 
% function [Att1,Att1_temp] = gather_att(data,postabs,maxv,dsrate)
% data_reshape = reshape(data,size(data,1),dsrate,size(data,2)/dsrate,1,size(data,4),size(data,5));
% % data_temp = max(data_reshape,[],3);
% data_temp = repmat(data_reshape,[1,2,1,1,1,1]);
% if ~maxv
%     if ~postabs
%         Att1 = mean(mean(squeeze(data),2),1); 
%         Att1_temp = mean(squeeze(data_temp),1);
%         Att1_temp = mean(Att1_temp,2);
%         sizes = size(Att1_temp);
%         Att1_temp = reshape(Att1_temp,[sizes(1),sizes(3:end)]);
%     else
%         Att1 = mean(mean(squeeze(data),2),1);
%         Att1_temp = mean(squeeze(data_temp),1);
%         Att1_temp = mean(Att1_temp,2);
%         sizes = size(Att1_temp);
%         Att1_temp = reshape(Att1_temp,[sizes(1),sizes(3:end)]);
%     end
% else
%     if ~postabs
%         Att1 = max(mean(squeeze(data),1),[],2);
%         Att1_temp = mean(squeeze(data_temp),1);
%         Att1_temp = max(Att1_temp,[],2);
%         sizes = size(Att1_temp);
%         Att1_temp = reshape(Att1_temp,[sizes(1),sizes(3:end)]);
%     else
%         Att1 = max(mean(squeeze(data),1),[],2);
%         Att1_temp = mean(squeeze(data_temp),1);
%         Att1_temp = max(Att1_temp,[],2);
%         sizes = size(Att1_temp);
%         Att1_temp = reshape(Att1_temp,[sizes(1),sizes(3:end)]);
%     end
% end
% end
% 
% function [Att1] = gather_att2(data,postabs,maxv)
% % Att1 = mean(mean(abs(squeeze(data)),1),2);
% Att1 = mean(mean(squeeze(data),2),1);
% end

function [Att1,Att1_temp] = gather_att(data,postabs,maxv,dsrate,tau)
if nargin<5
    tau=false;
end
if tau
    data_reshape = reshape(data,size(data,1),dsrate,size(data,2)/dsrate,1,size(data,4),size(data,5),size(data,6));
    % data_temp = max(data_reshape,[],3);
    data_temp = repmat(data_reshape,[1,2,1,1,1,1,1]);
else
    data_reshape = reshape(data,size(data,1),dsrate,size(data,2)/dsrate,1,size(data,4),size(data,5),size(data,6));
    % data_temp = max(data_reshape,[],3);
    data_temp = repmat(data_reshape,[1,2,1,1,1,1]);
end
if ~maxv
    if ~postabs
        Att1 = mean(mean(abs(squeeze(data)),2),1); 
        Att1_temp = mean(abs(squeeze(data_temp)),1);
        Att1_temp = mean(Att1_temp,2);
        sizes = size(Att1_temp);
        Att1_temp = reshape(Att1_temp,[sizes(1),sizes(3:end)]);
    else
        Att1 = mean(abs(mean(squeeze(data),2)),1);
%         Att1 = mean(mean(abs(mean(squeeze(data),2)),1),5);
        Att1_temp = abs(mean(squeeze(data_temp),1));
        Att1_temp = mean(Att1_temp,2);
        sizes = size(Att1_temp);
        Att1_temp = reshape(Att1_temp,[sizes(1),sizes(3:end)]);
    end
else
    if ~postabs
        Att1 = max(mean(abs(squeeze(data)),1),[],2);
        Att1_temp = mean(abs(squeeze(data_temp)),1);
        Att1_temp = max(Att1_temp,[],2);
        sizes = size(Att1_temp);
        Att1_temp = reshape(Att1_temp,[sizes(1),sizes(3:end)]);
    else
        Att1 = max(abs(mean(squeeze(data),1)),[],2);
        Att1_temp = abs(mean(squeeze(data_temp),1));
        Att1_temp = max(Att1_temp,[],2);
        sizes = size(Att1_temp);
        Att1_temp = reshape(Att1_temp,[sizes(1),sizes(3:end)]);
    end
end
if tau
    Att1_temp = permute(Att1_temp,[1,2,5,3,4]);
    sz = size(Att1_temp);
    Att1_temp = reshape(Att1_temp,[sz(1),sz(2)*sz(3),sz(4:end)]);
end
end
function [Att1] = gather_att2(data,postabs,maxv,tau)
if nargin<4
    tau=false;
end
% Att1 = mean(mean(abs(squeeze(data)),1),2);
Att1 = mean(abs(mean(squeeze(data),2)),1);
if tau
    Att1_temp = permute(Att1_temp,[1,5,2,3,4]);
    sz = size(Att1_temp);
    Att1_temp = reshape(Att1_temp,[sz(1)*sz(2),sz(3:end)]);
end
end
function [Att1,Att1_temp] = gather_att3(data,postabs,maxv,dsrate,tau)
if nargin<5
    tau=false;
end
if tau
    data_reshape = reshape(data,size(data,1),dsrate,size(data,2)/dsrate,1,size(data,4),size(data,5),size(data,6));
    data_temp = repmat(data_reshape,[1,2,1,1,1,1,1]);
else
    data_reshape = reshape(data,size(data,1),dsrate,size(data,2)/dsrate,1,size(data,4),size(data,5));
    data_temp = repmat(data_reshape,[1,2,1,1,1,1]);
end
Att1 = mean(mean(squeeze(data),2),1); 
Att1 = max(Att1,0);
Att1_temp = mean(squeeze(data_temp),1);
Att1_temp = mean(Att1_temp,2);
sizes = size(Att1_temp);
Att1_temp = reshape(Att1_temp,[sizes(1),sizes(3:end)]);
Att1_temp = max(Att1_temp,0);
end
function [data,data_temp] = occ2att(select_all,org,occ,docorr)
%     select_all = reshape(select_all,[size(select_all,1),size(select_all,2)*size(select_all,3)]);
%     org = reshape(org,[size(org,1),size(select_all,2)*size(select_all,3)]);
%     occ = reshape(occ,[size(occ,1),size(occ,2)*size(occ,3),15,15]);
%     cc_all = diag(corr(select_all',org','type','Pearson'));
%     cc_occ = zeros(15,15);
%     if docorr
%         for i=1:15
%             for j=1:15
%                 occ_ = squeeze(occ(:,:,i,j));
%                 cc_occ(i,j) = mean(abs(cc_all)-abs(diag(corr(occ_',org','type','Pearson'))));
%             end
%         end
%     else
%         cc_occ = mean(mean(squeeze(occ),1),2);
%     end
%     data = reshape(cc_occ,[1,1,15,15]);
%     data_temp = data;
    
    
    select_all = reshape(select_all,[1,size(select_all,1)*size(select_all,2)*size(select_all,3)]);
    org = reshape(org,[1,size(select_all,1)*size(select_all,2)*size(select_all,3)]);
    occ = reshape(occ,[1,size(occ,1)*size(occ,2)*size(occ,3),15,15]);
    cc_all = diag(corr(select_all',org','type','Pearson'));
    cc_occ = zeros(15,15);
    if docorr
        for i=1:15
            for j=1:15
                occ_ = reshape(occ(1,:,i,j),[1,size(occ,2)]);
                cc_occ(i,j) = mean(abs(cc_all)-abs(diag(corr(occ_',org','type','Pearson'))));
            end
        end
    else
        cc_occ = mean(mean(squeeze(occ),1),2);
    end
    data = reshape(cc_occ,[1,1,15,15]);
    data_temp = data;
end
function [att,att_temp,min25,max75] = robust_rescale(data,data_temp,mask,median_baseline,regions,tau)
    if nargin<4
        median_baseline=true;
    end
    if nargin<5
        regions = false;
    end
    if nargin<6
        tau = false;
    end
    if regions
        att_temp = zeros(size(data_temp));
        att = zeros(size(data));
        region = reshape(regions(repmat(mask,1,1,20)==1),length(find(mask==1)),20);
        unique_region = unique(region,'rows');
        sz = size(regions);
        mask_region = zeros(sz(1:length(sz)-1));
        for regs=1:size(unique_region)
            for i =1:size(regions,1)
                for j=1:size(regions,2)
                    if strcmp(reshape(regions(i,j,:),1,length(unique_region(regs,:))),unique_region(regs,:))
                        mask_region(i,j) = regs;
                    end
                end
            end 
        end
        for regs=1:size(unique_region)
            mask_temp = repmat(reshape(mask_region,[1,1,size(mask_region,1),size(mask_region,2)]),size(data_temp,1),size(data_temp,2),1,1)==regs;
            mask_rep = repmat(reshape(mask_region,[1,1,size(mask_region,1),size(mask_region,2)]),size(data,1),size(data,2),1,1)==regs;
            sorted = sort(data(mask_rep));
            min25 = sorted(int32(0.05*length(sorted))+1);
            max75 = sorted(int32(0.95*length(sorted)));
            md = median(sorted);
            
            sorted_temp = sort(data_temp(mask_temp));
            min25_temp = sorted_temp(int32(0.05*length(sorted_temp))+1);
            max75_temp = sorted_temp(int32(0.95*length(sorted_temp)));
            md_temp = median(sorted_temp);
            if median_baseline
                att(mask_rep) = (data(mask_rep)-md)/(max75-min25);
                att_temp(mask_temp) = (data_temp(mask_temp)-md)/(max75-min25);
%                 att_temp(mask_temp) = (data_temp(mask_temp)-md_temp)/(max75_temp-min25_temp);
            else
                att = (data-min25)/(max75-min25);
                att_temp(mask_temp) = (data_temp(mask_temp)-min25)/(max75-min25);
%                 att_temp(mask_temp) = (data_temp(mask_temp)-min25_temp)/(max75_temp-min25_temp);
            end
        end
    else
%         if tau
%             sorted = sort(data(repmat(reshape(mask,[1,1,size(mask,1),size(mask,2),1]),size(data,1),size(data,2),1,1,size(data,5))==1));
%         else
%             sorted = sort(data(repmat(reshape(mask,[1,1,size(mask,1),size(mask,2)]),size(data,1),size(data,2),1,1)==1));
%         end
        sorted = sort(data(repmat(reshape(mask,[1,1,size(mask,1),size(mask,2)]),size(data,1),size(data,2),1,1)==1 & ~isnan(data)));
%         min25 = sorted(int32(0.05*length(sorted))+1);
%         max75 = sorted(int32(0.95*length(sorted)));
        min25 = sorted(int32(0.01*length(sorted))+1);
        max75 = sorted(int32(0.99*length(sorted)));
        md = median(sorted);
        
%         if tau
%             sorted_temp = sort(data_temp(repmat(reshape(mask,[1,1,size(mask,1),size(mask,2),1]),size(data_temp,1),size(data_temp,2),1,1,size(data_temp,5))==1));
%         else
%             sorted_temp = sort(data_temp(repmat(reshape(mask,[1,1,size(mask,1),size(mask,2)]),size(data_temp,1),size(data_temp,2),1,1)==1));
%         end
        sorted_temp = sort(data_temp(repmat(reshape(mask,[1,1,size(mask,1),size(mask,2)]),size(data_temp,1),size(data_temp,2),1,1)==1 & ~isnan(data_temp)));
%         min25_temp = sorted_temp(int32(0.05*length(sorted_temp))+1);
%         max75_temp = sorted_temp(int32(0.95*length(sorted_temp)));
        min25_temp = sorted_temp(int32(0.01*length(sorted_temp))+1);
        max75_temp = sorted_temp(int32(0.99*length(sorted_temp)));
        md_temp = median(sorted_temp);

        if median_baseline
            att = (data-md)/(max75-min25);
%             att_temp = (data_temp-md)/(max75-min25);
            att_temp = (data_temp-md_temp)/(max75_temp-min25_temp);
            att_temp(isnan(att_temp))=0;
            att_temp(repmat(reshape(mask,[1,1,size(mask,1),size(mask,2)]),size(data_temp,1),size(data_temp,2),1,1)==0) = 0;
        else
            att = (data-min25)/(max75-min25);
%             att_temp = (data_temp-min25)/(max75-min25);
            att_temp = (data_temp-min25_temp)/(max75_temp-min25_temp);
            att_temp(isnan(att_temp))=0;
            att_temp(repmat(reshape(mask,[1,1,size(mask,1),size(mask,2)]),size(data_temp,1),size(data_temp,2),1,1)==0) = 0;
        end
    end
end

% function [att,att_temp,min25,max75] = robust_rescale(data,data_temp,mask,median_baseline,regions,tau)
%     att = data;
%     att(isnan(att))=0;
%     att_temp = data_temp;
%     att_temp(isnan(att_temp))=0;
%     min25 = min(att(:));
%     max75 = max(att(:));
% end

function [diff] = checkaud(diff,mask,regions)
    sz = size(regions);
    aud = zeros(sz(1:length(sz)-1));
    for i =1:size(regions,1)
        for j=1:size(regions,2)
            if contains(reshape(regions(i,j,:),1,size(regions,3)),'mSTG') || contains(reshape(regions(i,j,:),1,size(regions,3)),'cSTG') && mask(i,j) == 1
                aud(i,j) = true;
            end
        end
    end
    aud_rep = repmat(reshape(aud,[1,1,size(aud,1),size(aud,2)]),size(diff,1),size(diff,2),1,1);
    auds = find(aud_rep==1);
    auds = auds(diff(auds)<0);
    diff(auds) = 0;
%     if pre_max>0
%         pre = pre-pre_max;
%     end
end

function [frames,frames_brain,err,tticks] = attr_vis_tau(attr,post,visdiff,postabs,max_att,plot_on_same_brain,above_median,temp,wreg,entire_period,return_curve,investreg)
%     maxminp1=[];
%     maxminp2=[10000];
%     for sub =1:SUB
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/NY',subjs{sub},'_elec.mat']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/att_NY',subjs{sub},'']);
%         att = eval(attrs{a});
%         if strcmp(attrs{a},'freq_formants_hamon')
%             maxp1 = gather_att(att(ind,on:off,1,:,:),postabs,max_att);
%             maxp2 = gather_att(att(ind,on:off,2,:,:),postabs,max_att);
%             maxminp1 = [maxminp1,max(maxp1(:))];
%             maxminp2 = [maxminp2,max(maxp2(:))];
%         else
%             maxp1 = gather_att(att(ind,on:off,1,:,:),postabs,max_att);
%             maxminp1 = [maxminp1,max(maxp1(:))];
%         end
%     end
%     maxminp1 = min(maxminp1);
%     maxminp2 = min(maxminp2);
    % choose a color map and set indecies of colors for each electrode
    color_map = hot(256);%afmhot;%hot;
    color_map = color_map(end:-1:1,:);
%     [color_map]=cbrewer('seq', 'Blues', 256);
    [color_map_div]=cbrewer('div', 'RdBu', 256,'PCHIP');
%     color_map_ = getPyPlot_cMap('afmhot',256);color_map_=color_map_(end:-1:1,:); color_map = color_map_;
%     [color_map_]=cbrewer('seq', 'Greens', 256,'PCHIP'); color_map = color_map_;
    [colorset1] = cbrewer('qual', 'Set1', 9,'PCHIP');
    color_map_div = color_map_div(end:-1:1,:);
    color_ind = value2colorind(ones(128,1),'hot',[0,1]); 
    color_ind = [color_ind,color_ind,color_ind];
    colors = {color_ind, color_map};
    clrmap = color_map;
    if strcmp(attr,'ecog')
        dsrate_org = 1;
    else
        dsrate_org = 8;
    end
    if entire_period
        on = 16/dsrate_org+1;
        off = 120/dsrate_org;
    else
        if post
            on = (16+32)/dsrate_org+1;
            off = 120/dsrate_org;
        else
            on = 16/dsrate_org+1;
            off = (16+32)/dsrate_org;
            on_post = (16+32)/dsrate_org+1;
            off_post = 120/dsrate_org;
        end
    end
    
    on_ttau=1;
    off_ttau=144;
    dsrate = 8/dsrate_org;
    % ind = [1:10,21:30];
    ind = [1:2];%[1:50];
    dumm = 1;
    Att1_cell = {};
    Att1_nt_cell = {};
    Att2_cell = {};
    att_cell_c = {};
    att_cell_cpre = {};
    att_cell_a = {};
    att_cell_p = {};
    att_cell_n = {};
    diff_cell = {};
    diff2_cell = {};
    mask_cell = {};
    coord_cell = {};
    region_cell = {};
    subjs = {'717','742','749'};
%     subjs = {'717'};
    SUB = length(subjs);
    root_dir = '/Users/ranwang/Documents/writen_paper/NER2020/Visuallization_matlab';
    use_occ = false;%(strcmp(attr,'amplitudes') || strcmp(attr,'f0_hz') || strcmp(attr,'freq_formants_hamon')) && ~temp && ~visdiff && post;
    use_occentire = false;%strcmp(attr,'freq_formants_hamon') && ~temp && ~visdiff && post;
    for sub=1:SUB
        
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/att_NY',subjs{sub},'']);
%         if use_occentire
%             load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_cc_value_',subjs{sub},'']); 
%             ind = [1:2];
%         else
%             ind = [1:50];
%             if use_occ
%                 load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_dict_IG_value_',subjs{sub},'']);   
%     %             load(['/Users/ranwang/Documents/writen_paper/NER2020/attr_dict_IG_value_',subjs{sub},'_nonnoise_']);
%             else
%                 load(['/Users/ranwang/Documents/writen_paper/NER2020/attr_dict_IG_value_',subjs{sub},'_nonnoise_']);
%     %             load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_dict_IG_value_',subjs{sub},'']);
%             end
%         end
        load(['/Users/ranwang/Documents/writen_paper/NER2020/NY_',subjs{sub},'_elec_entiregrid.mat']);
%         onregion = {'cSTG','mSTG','parstriangularis','parsopercularis','precentral','postcentral','inferiorparietal','supramarginal'};
%         onregion = {'parstriangularis','parsopercularis','precentral','postcentral'};
        onregion = {'cSTG','mSTG','parstriangularis','parsopercularis','precentral','postcentral'};
        region_cell{sub} = regions;
%         mask = isregion(regions,onregion);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_woaud/','NY',subjs{sub},'_elec_woaud.mat']);

        load(['/Users/ranwang/greene/neural_decoding/code/cnn/ALAE/training_artifacts/entiregrid_',subjs{sub},'_han5amppowerloss_alphasup3_percept_singlebandwidth_fractive1500lim_bgnoiselearn_ecogfinetune_anticausal_universal_passive_groupnormxdim_correctbaseline/grad_cov_temporal_8.mat',]);
        eval([attr,'_passive=',attr,';']);
        load(['/Users/ranwang/greene/neural_decoding/code/cnn/ALAE/training_artifacts/entiregrid_',subjs{sub},'_han5amppowerloss_alphasup3_percept_singlebandwidth_fractive1500lim_bgnoiselearn_ecogfinetune_anticausal_universal_imagine_groupnormxdim_correctbaseline/grad_cov_temporal_8.mat',]);
        eval([attr,'_imagine=',attr,';']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/entiregrid_',subjs{sub},'_normxdim_anticausal_covnorm_NOVISWEIGHT_step1/occ_cov_temporal_8.mat',]);
        load(['/Users/ranwang/greene/neural_decoding/code/cnn/ALAE/training_artifacts/entiregrid_',subjs{sub},'_han5amppowerloss_alphasup3_percept_singlebandwidth_fractive1500lim_bgnoiselearn_ecogfinetune_anticausal_universal_active_groupnormxdim_correctbaseline/grad_cov_temporal_8.mat',]);
        eval([attr,'_active=',attr,';']);
        
        load(['/Users/ranwang/Documents/writen_paper/NER2020/entiregrid_',subjs{sub},'_random/occ_cov_temporal_8.mat',]);
        eval([attr,'_noiselevel=',attr,';']);
        load(['/Users/ranwang/Documents/writen_paper/NER2020/entiregrid_',subjs{sub},'_random/occ_covnorm_temporal_tstep1_tau_8.mat',]);
        eval([attr,'_noisetemporal=',attr,';']);

        load(['/Users/ranwang/Documents/writen_paper/NER2020/ecog_',subjs{sub},'prod.mat']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_dict_IG_value_',subjs{sub},'']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_dict_IG_value_',subjs{sub},'_prearti']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_cc_value_prearti_',subjs{sub}]);
        load(['/Users/ranwang/Documents/writen_paper/NER2020/percept_passive/',subjs{sub},'/occ_covnorm_temporal_tstep1_tau_8.mat',]);
%         load(['/Users/ranwang/greene/neural_decoding/code/cnn/ALAE/training_artifacts/entiregrid_',subjs{sub},'_han5amppowerloss_alphasup3_percept_singlebandwidth_fractive1500lim_bgnoiselearn_ecogfinetune_anticausal_imagine_groupnormxdim_correctbaseline/occ_covnorm_temporal_tstep1_tau_8.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/entiregrid_',subjs{sub},'_normxdim_anticausal_covnorm_NOVISWEIGHT_step1_grad/grad_l2_temporal_tstep1_tau_8.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_/occ_cc_value.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_wpariental/occ_cc_value.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_wpariental/occ_cc_temporal_8.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_woaud/occ_cc_temporal_8.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_woaud/occ_cc_value.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_prearti/occ_cc_value.mat',]);
        eval([attr,'=',attr,';']);
        causal=false;
        ref = 'passive';%'passive';'imagine';'active';
        thepower=0.43;
        do_t_ds=false;
        Subj = ['NY',subjs{sub}];
        % get paths to all visualization files
        paths = get_path(root_dir,Subj);
        channel_info_all = get_channel_info(paths);
        plot_data_all_mni = get_plot_data(paths, 'mni', 'lh');
%         plot_data_all_T1 = get_plot_data(paths, 'subj', 'lh');
    %     coord = T1;
        coord = mni;
        aud = isaud(regions);
        mask_nonaud = mask & ~aud;
        onset_diff = 0;%onset-continuouse;
        if strcmp(attr,'freq_formants_hamon')
            att = eval(attr);
            if do_t_ds
                att = att(:,1:8:end,:,:,:,:);
            end
            data = att(:,on_ttau:off_ttau,1,:,:,:);
            att_nt = eval([attr,'_noisetemporal']);
            if do_t_ds
                att_nt = att_nt(:,1:8:end,:,:,:,:);
            end
            data_nt = att_nt(ind,on_ttau:off_ttau,1,:,:,:);
            
            att_causal = eval([attr,'_passive']);
            data_c = att_causal(:,on:off,1,:,:);
            data_cpre = att_causal(:,1:4,1,:,:);
            att_anticausal = eval([attr,'_imagine']);
            data_a = att_anticausal(:,on:off,1,:,:);
            att_percept = eval([attr,'_active']);
            data_p = att_percept(:,on:off,1,:,:);
            att_noise = eval([attr,'_noiselevel']);
            data_n = att_noise(:,on:off,1,:,:);

            [Att1_c,Att1_temp_c] = gather_att(data_c,postabs,max_att,dsrate);
            [Att1_cpre,Att1_temp_cpre] = gather_att(data_c,postabs,max_att,dsrate);
            [Att1_a,Att1_temp_a] = gather_att(data_a,postabs,max_att,dsrate);
            [Att1_p,Att1_temp_p] = gather_att(data_p,postabs,max_att,dsrate);
            [Att1_a,Att1_temp_a,m25,m75] = robust_rescale(Att1_a,Att1_temp_a,mask,above_median);
            Att1_c = (Att1_c-m25)/(m75-m25);
            Att1_cpre = (Att1_cpre-m25)/(m75-m25);
            [Att1_p,Att1_temp_p,m25,m75] = robust_rescale(Att1_p,Att1_temp_p,mask,above_median);
%             Att1_p = (Att1_p-m25)/(m75-m25);
            
            Att1_a = (max(Att1_a,0)).^thepower;
            Att1_c = (max(Att1_c,0)).^thepower;
            Att1_cpre = (max(Att1_cpre,0)).^thepower;
            Att1_p = (max(Att1_p,0)).^thepower;
            
            att_cell_c{sub} = Att1_c;
            att_cell_cpre{sub} = Att1_cpre;
            att_cell_a{sub} = Att1_a;
            att_cell_p{sub} = Att1_p;
            
            if use_occ
                Att1 = gather_att2(data,postabs,max_att);
                Att1_temp = Att1;
            else
                [Att1,Att1_temp] = gather_att(data,postabs,max_att,dsrate,true);
                [Att1_nt,Att1_temp_nt] = gather_att(data_nt,postabs,max_att,dsrate,true);
            end
            sorted = sort(Att1(:));
            % get all plot frames
            if visdiff && ~post
                continue;
            else
%                 [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map,sorted(end-2));
                if wreg
                    [Att1_md,Att1_md_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,true,regions,true);
                    [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median,regions,true);
                else
                    [Att1_md,Att1_md_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,true,false,true);
                    [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median,false,true);
%                     Att1_temp = (max(Att1_temp,0)).^thepower;
                end

                if temp
                    Att1_cell{sub} = Att1_temp;
                    Att1_nt_cell{sub} = Att1_temp_nt;
                else
                    Att1_cell{sub} = Att1;
                    Att1_nt_cell{sub} = Att1_nt;
                end
                
                mask_cell{sub} = mask;
                coord_cell{sub} = coord;
                if plot_on_same_brain 
                    if sub == SUB
                        if ~return_curve
                            [frames] = VisualAtt(Att1_cell, mask_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false);
                        else
                            [frames,frames_brain,err,tticks] = PlotAtt(Att1_cell, mask_cell,coord_cell,region_cell,causal,plot_data_all_mni,clrmap,att_cell_c,att_cell_cpre,att_cell_a,att_cell_p,Att1_nt_cell,ref,investreg);
                        end
%                         axes(haxis(sub2ind([col,row],dumm,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
                    end
                else
                    [frames] = VisualAtt({Att1}, {mask}, {coord}, plot_data_all_mni,color_map,1,false,false);
%                     axes(haxis(sub2ind([col,row],dumm,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
                end
            end
            
            
            data = att(:,on_ttau:off_ttau,2,:,:,:);
%             data_nt = att_noisetemporal(:,on_ttau:off_ttau,2,:,:,:);
            
            att_causal = eval([attr,'_passive']);
            data_c = att_causal(:,on:off,2,:,:);
            data_cpre = att_causal(:,1:4,2,:,:);
            att_anticausal = eval([attr,'_imagine']);
            data_a = att_anticausal(:,on:off,2,:,:);
            att_percept = eval([attr,'_active']);
            data_p = att_percept(:,on:off,2,:,:);
            att_noise = eval([attr,'_noiselevel']);
            data_n = att_noise(:,on:off,2,:,:);

            [Att1_c,Att1_temp_c] = gather_att(data_c,postabs,max_att,dsrate);
            [Att1_cpre,Att1_temp_c] = gather_att(data_cpre,postabs,max_att,dsrate);
            [Att1_a,Att1_temp_a] = gather_att(data_a,postabs,max_att,dsrate);
            [Att1_p,Att1_temp_p] = gather_att(data_p,postabs,max_att,dsrate);
            [Att1_a,Att1_temp_a,m25,m75] = robust_rescale(Att1_a,Att1_temp_a,mask,above_median);
            Att1_c = (Att1_c-m25)/(m75-m25);
            Att1_cpre  = (Att1_cpre-m25)/(m75-m25);
            [Att1_p,Att1_temp_p,m25,m75] = robust_rescale(Att1_p,Att1_temp_p,mask,above_median);
%             Att1_p = (Att1_p-m25)/(m75-m25);
            
            Att1_a = (max(Att1_a,0)).^thepower;
            Att1_c = (max(Att1_c,0)).^thepower;
            Att1_cpre = (max(Att1_cpre,0)).^thepower;
            Att1_p = (max(Att1_p,0)).^thepower;
            
            att_cell_c{sub} = Att1_c;
            att_cell_cpre{sub} = Att1_cpre;
            att_cell_a{sub} = Att1_a;
            att_cell_p{sub} = Att1_p;
            
            if use_occ
                Att2 = gather_att2(data,postabs,max_att);
                Att2_temp = Att2;
            else
                [Att2,Att2_temp] = gather_att(data,postabs,max_att,dsrate,true);
                [Att2_nt,Att2_temp_nt] = gather_att(data_nt,postabs,max_att,dsrate,true);
            end
%             [Att2,Att2_temp] = gather_att(data,postabs,max_att,dsrate);
            sorted = sort(Att2(:));
            if visdiff && ~post
               continue;
            else
%                 [frames] = VisualAtt(Att2, mask, coord, plot_data_all_mni,color_map,sorted(end-2));
                if wreg
                    [Att2_md,Att2_md_temp,m25,m75] = robust_rescale(Att2,Att2_temp,mask,true,regions,true);
                    [Att2,Att2_temp,m25,m75] = robust_rescale(Att2,Att2_temp,mask,above_median,regions,true);
                else
                    [Att2_md,Att2_md_temp,m25,m75] = robust_rescale(Att2,Att2_temp,mask,true,false,true);
                    [Att2,Att2_temp,m25,m75] = robust_rescale(Att2,Att2_temp,mask,above_median,false,true);
%                     Att2_temp = (max(Att2_temp,0)).^thepower;
                end
                if temp
                    Att2_cell{sub} = Att2_temp;
                    Att2_nt_cell{sub} = Att2_temp_nt;
                else
                    Att2_cell{sub} = Att2;
                    Att2_nt_cell{sub} = Att2_nt;
                end
                mask_cell{sub} = mask;
                coord_cell{sub} = coord;
                if plot_on_same_brain 
                    if sub == SUB
                        if ~return_curve
                            [frames2] = VisualAtt(Att2_cell, mask_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false);
                        else
                            [frames2,frames2_brain,err2,tticks] = PlotAtt(Att2_cell, mask_cell,coord_cell,region_cell,causal,plot_data_all_mni,clrmap,att_cell_c,att_cell_cpre,att_cell_a,att_cell_p,Att2_nt_cell,ref,investreg);
                        end
%                         axes(haxis(sub2ind([col,row],dumm+1,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
                    end
                else
                    [frames2] = VisualAtt({Att2}, {mask}, {coord}, plot_data_all_mni,color_map,1,false,false);
%                     axes(haxis(sub2ind([col,row],dumm+1,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
                end
            end

            
%             coef = find_normcoef(Att2,Att1,mask,false,20);
%             diff = Att2*coef(1)+coef(2) - Att1;
%             [Att1,Att1_temp] = gather_att(att(ind,on:off,1,:,:),postabs,max_att,dsrate);
%             [Att2,Att2_temp] = gather_att(att(ind,on:off,2,:,:),postabs,max_att,dsrate);
            if temp
                diff = Att2_temp-Att1_temp;
            else
                diff = Att2-Att1;
            end
%             sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
            if temp
                mask_temp = repmat(reshape(mask,[1,1,15,15]),[1,size(diff,2),1,1]);
                ind_active = find(mask_temp==1);
            else
                ind_active = find(mask==1);
            end
            sorted = sort(diff(ind_active)); maxv = max(sorted(end-1),abs(sorted(2)));%maxv = min(sorted(end-1),abs(sorted(2)));
            diff = diff/maxv;
            diff = (diff+1)/2;
            if temp
                diff2_cell{sub} = ((Att2_md_temp-Att1_md_temp)/maxv+1)/2;
            else
                diff2_cell{sub} = diff;
            end
            mask_cell{sub} = mask;
            coord_cell{sub} = coord;
            if plot_on_same_brain 
                if sub == SUB
                    if ~return_curve
                        [frames3] = VisualAtt(diff2_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
                    else
%                         [frames3,frames3_brain] = PlotAtt(diff2_cell, mask_cell,coord_cell,region_cell,causal,plot_data_all_mni,clrmap);
                    end
%                     axes(haxis(sub2ind([col,row],dumm+2,1)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
                end
            else
                [frames3] = VisualAtt({diff}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                 axes(haxis(sub2ind([col,row],dumm+2,sub)));
%                 imshow(uint8(squeeze(frames)));
%                 axis off;
            end
            if sub == SUB
                if ~return_curve
                    frames = cat(4,frames,frames2,frames3);
                else
                    frames = cat(3,frames,frames2);
                    err = cat(3,err,err2);
                    frames_brain = cat(3,frames_brain,frames2_brain);
                end
            end
            
        else
            if strcmp(attr,'whatever')
                continue;
            else
                if strcmp(attr,'onset')
                    att = onset_diff;
                else
                    att = eval(attr);
                    if do_t_ds
                        att = att(:,1:8:end,:,:,:,:);
                    end
                    att_nt = eval([attr,'_noisetemporal']);
                    if do_t_ds
                        att_nt = att_nt(:,1:8:end,:,:,:,:);
                    end
                end
                data = att(ind,on_ttau:off_ttau,1,:,:,:);
                data_nt = att_nt(ind,on_ttau:off_ttau,1,:,:,:);
                
                att_causal = eval([attr,'_passive']);
                data_c = att_causal(:,on:off,1,:,:);
                data_cpre = att_causal(:,1:4,1,:,:);
                att_anticausal = eval([attr,'_imagine']);
                data_a = att_anticausal(:,on:off,1,:,:);
                att_percept = eval([attr,'_active']);
                data_p = att_percept(:,on:off,1,:,:);
                att_noise = eval([attr,'_noiselevel']);
                data_n = att_noise(:,on:off,1,:,:);

                [Att1_c,Att1_temp_c] = gather_att(data_c,postabs,max_att,dsrate);
                [Att1_cpre,Att1_temp_cpre] = gather_att(data_cpre,postabs,max_att,dsrate);
                [Att1_a,Att1_temp_a] = gather_att(data_a,postabs,max_att,dsrate);
                [Att1_p,Att1_temp_p] = gather_att(data_p,postabs,max_att,dsrate);
                [Att1_n,Att1_temp_n] = gather_att(data_n,postabs,max_att,dsrate);
                
                [Att1_a,Att1_temp_a,m25,m75] = robust_rescale(Att1_a,Att1_temp_a,mask,above_median);
                Att1_c = (Att1_c-m25)/(m75-m25);
                Att1_cpre = (Att1_cpre-m25)/(m75-m25);
                Att1_n = (Att1_n-m25)/(m75-m25);
%                 [Att1_p,Att1_temp_p,m25,m75] = robust_rescale(Att1_p,Att1_temp_p,mask,above_median);
                Att1_p = (Att1_p-m25)/(m75-m25);
                
                Att_ac = (Att1_a-Att1_c)./(abs(Att1_a)+abs(Att1_c));
                Att_ccp = (Att1_c-Att1_cpre)./(abs(Att1_c)+abs(Att1_cpre));
                att_cell_ac{sub} = Att_ac;
                att_cell_ccp{sub} = Att_ccp;
                             
                
                global noiselevel;
                localnoise=noiselevel.^(1/thepower);
                Att_an = (Att1_a-localnoise)./(abs(Att1_a)+localnoise);
                Att_cn = (Att1_c-localnoise)./(abs(Att1_c)+localnoise);
                Att_cpn = (Att1_cpre-localnoise)./(abs(Att1_cpre)+localnoise);
                att_cell_an{sub} = Att_an;
                att_cell_cn{sub} = Att_cn;
                att_cell_cpn{sub} = Att_cpn;
                
                
                Att1_a = (max(Att1_a,0)).^thepower;
                Att1_c = (max(Att1_c,0)).^thepower;
                Att1_cpre = (max(Att1_cpre,0)).^thepower;
                Att1_n = (max(Att1_n,0)).^thepower;
                Att1_p = (max(Att1_p,0)).^thepower;
                
                
                att_cell_c{sub} = Att1_c;
                att_cell_cpre{sub} = Att1_cpre;
                att_cell_a{sub} = Att1_a;
                att_cell_p{sub} = Att1_p;
                att_cell_n{sub} = Att1_n;
                
                sorted = sort(data(:));
                sorted = sorted(end-1);
                if use_occ && (strcmp(attr,'f0_hz') || strcmp(attr,'amplitudes'))
                    Att1 = gather_att2(data,postabs,max_att);
                    Att1_temp = Att1;
                else
                    if strcmp(attr,'ecog')
                        [Att1,Att1_temp] = gather_att3(data,postabs,max_att,dsrate,true);
                    else
                        [Att1,Att1_temp] = gather_att(data,postabs,max_att,dsrate,true);
                        [Att1_nt,Att1_temp_nt] = gather_att(data_nt,postabs,max_att,dsrate,true);
                    end
                end
                sorted = sort(Att1(:));
%                 sorted = sorted(end-1);
                % get all plot frames
                if strcmp(attr,'loudness')% || strcmp(attr,'f0_hz')
                    if visdiff && ~post
                       continue;
                    else
%                         [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map,maxminp1);
                        [Att1_nt,Att1_temp_nt,m25,m75] = robust_rescale(Att1_nt,Att1_temp_nt,mask,above_median,false,true);
                        if wreg
                            [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median,regions,true);
                        else
                            [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median,false,true);
%                             Att1_temp = (max(Att1_temp,0)).^thepower;
                        end
                        if temp
                            Att1_cell{sub} = Att1_temp;
                            Att1_nt_cell{sub} = Att1_temp_nt;
                        else
                            Att1_cell{sub} = Att1;
                            Att1_nt_cell{sub} = Att1_nt;
                        end
                        mask_cell{sub} = mask;
                        coord_cell{sub} = coord;
                        if plot_on_same_brain 
                            if sub == SUB
                                if ~return_curve
                                    [frames] = VisualAtt(Att1_cell, mask_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false);
                                else
                                    [frames,frames_brain,err,tticks] = PlotAtt(Att1_cell, mask_cell,coord_cell,region_cell,causal,plot_data_all_mni,clrmap,att_cell_c,att_cell_cpre,att_cell_a,att_cell_p,Att1_nt_cell,ref,investreg);
                                end
%                                 axes(haxis(sub2ind([col,row],dumm,1)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
                            end
                        else
                            [frames] = VisualAtt({Att1}, {mask}, {coord}, plot_data_all_mni,color_map,1,false,false);
%                             axes(haxis(sub2ind([col,row],dumm,sub)));
%                             imshow(uint8(squeeze(frames)));
%                             axis off;
                        end
                    end
                else
                    if strcmp(attr,'onset')
                        continue;
                    else
                        if visdiff && ~post
%                             continue;
                        else
%                             [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map);
                            if wreg
                                [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median,regions,true);
                            else
                                [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median,false,true);
                                [Att1_nt,Att1_temp_nt,m25,m75] = robust_rescale(Att1_nt,Att1_temp_nt,mask,above_median,false,true);
%                                 Att1_temp = (max(Att1_temp,0)).^thepower;
                            end
                            if temp
                                Att1_cell{sub} = Att1_temp;
                                Att1_nt_cell{sub} = Att1_temp_nt;
                            else
                                Att1_cell{sub} = Att1;
                                Att1_nt_cell{sub} = Att1_nt;
                            end
                            mask_cell{sub} = mask;
                            coord_cell{sub} = coord;
                            if plot_on_same_brain 
                                if sub == SUB
                                    if ~return_curve
                                        [frames] = VisualAtt(Att1_cell, mask_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false);
                                    else
                                        [frames,frames_brain,err,tticks] = PlotAtt(Att1_cell, mask_cell,coord_cell,region_cell,causal,plot_data_all_mni,clrmap,att_cell_c,att_cell_cpre,att_cell_a,att_cell_p,Att1_nt_cell,ref,investreg);
%                                         [frames,prob1,prob0] = VisualAtt(att_cell_cpre, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map_div,1,true,true,false,true);
                                    end
%                                     axes(haxis(sub2ind([col,row],dumm,1)));
%                                     imshow(uint8(squeeze(frames)));
%                                     axis off;
                                end
                            else
                                [frames] = VisualAtt({Att1}, {mask}, {coord}, plot_data_all_mni,color_map,1,false,false);
%                                 axes(haxis(sub2ind([col,row],dumm,sub)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
                            end
                        end
                    end
                end
            end
        end
%         if ~plot_on_same_brain
%             if dumm == 2
%                 ylabel(Subj)
%             end  
%         end
%         if plot_on_same_brain
%             if sub == SUB
%                 title(attr)
%             end
%         else
%             if sub ==1
%                 title(attr)
%             end
%         end
    end
%     if strcmp(attr,'freq_formants_hamon')
%         dumm = dumm+3;
%     else
%         dumm = dumm+1;
%     end
end

function [frames,Att1_cell] = attr_vis(attr,post,visdiff,postabs,max_att,plot_on_same_brain,above_median,temp,wreg,entire_period,avg_data,annotation)
    if ~exist('annotation','var')
        annotation=[];
    end
%     maxminp1=[];
%     maxminp2=[10000];
%     for sub =1:SUB
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/NY',subjs{sub},'_elec.mat']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/att_NY',subjs{sub},'']);
%         att = eval(attrs{a});
%         if strcmp(attrs{a},'freq_formants_hamon')
%             maxp1 = gather_att(att(ind,on:off,1,:,:),postabs,max_att);
%             maxp2 = gather_att(att(ind,on:off,2,:,:),postabs,max_att);
%             maxminp1 = [maxminp1,max(maxp1(:))];
%             maxminp2 = [maxminp2,max(maxp2(:))];
%         else
%             maxp1 = gather_att(att(ind,on:off,1,:,:),postabs,max_att);
%             maxminp1 = [maxminp1,max(maxp1(:))];
%         end
%     end
%     maxminp1 = min(maxminp1);
%     maxminp2 = min(maxminp2);
    % choose a color map and set indecies of colors for each electrode
    color_map = hot(256);%afmhot;%hot;
    color_map = color_map(end:-1:1,:);
%     [color_map]=cbrewer('seq', 'Blues', 256);
    [color_map_div]=cbrewer('div', 'RdBu', 256,'PCHIP');
    [color_map_cb]=cbrewer('seq', 'Greens', 256,'PCHIP');
%     color_map_cb = color_map_cb(end:-1:1,:);
%     color_map = color_map_cb;
%     if annotation
%         color_map = color_map_cb;
%     end
    [colorset1] = cbrewer('qual', 'Set1', 9,'PCHIP');
    color_map_div = color_map_div(end:-1:1,:);
    color_ind = value2colorind(ones(128,1),'hot',[0,1]); 
    color_ind = [color_ind,color_ind,color_ind];
    colors = {color_ind, color_map};
    clrmap = color_map;
    if strcmp(attr,'ecog')
        dsrate_org = 1;
    else
        dsrate_org = 8;
    end
    if entire_period
        on = 16/dsrate_org+1;
        off = 120/dsrate_org;
    else
        if post
            on = (16+32)/dsrate_org+1;
            off = 120/dsrate_org;
        else
            on = 16/dsrate_org+1;
            off = (16+32)/dsrate_org;
            on_post = (16+32)/dsrate_org+1;
            off_post = 120/dsrate_org;
        end
    end

    dsrate = 8/dsrate_org;
    % ind = [1:10,21:30];
    ind = [1:2];%[1:50];
    dumm = 1;
    Att1_cell = {};
    Att2_cell = {};
    diff_cell = {};
    diff2_cell = {};
    mask_cell = {};
    coord_cell = {};
    region_cell = {};
    subjs = {'717','742','749'};
%     subjs = {'749'};
    SUB = length(subjs);
    root_dir = '/Users/ranwang/Documents/writen_paper/NER2020/Visuallization_matlab';
    use_occ = false;%(strcmp(attr,'amplitudes') || strcmp(attr,'f0_hz') || strcmp(attr,'freq_formants_hamon')) && ~temp && ~visdiff && post;
    use_occentire = false;%strcmp(attr,'freq_formants_hamon') && ~temp && ~visdiff && post;
    for sub=1:SUB
        
        load(['/Users/ranwang/Documents/writen_paper/NER2020/att_NY',subjs{sub},'']);
%         if use_occentire
%             load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_cc_value_',subjs{sub},'']); 
%             ind = [1:2];
%         else
%             ind = [1:50];
%             if use_occ
%                 load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_dict_IG_value_',subjs{sub},'']);   
%     %             load(['/Users/ranwang/Documents/writen_paper/NER2020/attr_dict_IG_value_',subjs{sub},'_nonnoise_']);
%             else
%                 load(['/Users/ranwang/Documents/writen_paper/NER2020/attr_dict_IG_value_',subjs{sub},'_nonnoise_']);
%     %             load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_dict_IG_value_',subjs{sub},'']);
%             end
%         end
        load(['/Users/ranwang/Documents/writen_paper/NER2020/NY_',subjs{sub},'_elec_entiregrid.mat']);
%         onregion = {'cSTG','mSTG','parstriangularis','parsopercularis','precentral','postcentral','inferiorparietal','supramarginal'};
        onregion = {'cSTG','mSTG','parstriangularis','parsopercularis','precentral','postcentral'};
%         onregion = {'inferiorparietal','supramarginal'};
        region_cell{sub} = regions;
%         mask = isregion(regions,onregion);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_woaud/','NY',subjs{sub},'_elec_woaud.mat']);

        load(['/Users/ranwang/Documents/writen_paper/NER2020/entiregrid_ecog_',subjs{sub},'percept_imagine.mat']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_dict_IG_value_',subjs{sub},'']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_dict_IG_value_',subjs{sub},'_prearti']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_cc_value_prearti_',subjs{sub}]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_/occ_cc_temporal_8.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_/occ_cc_value.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_wpariental/occ_cc_value.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_wpariental/occ_cc_temporal_8.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/entiregrid_',subjs{sub},'_normxdim_anticausal_covnorm_NOVISWEIGHT_step1/occ_cov_temporal_8.mat',]);
        load(['/Users/ranwang/Documents/writen_paper/NER2020/entiregrid_',subjs{sub},'_normxdim_anticausal_percept_active_covnorm_NOVISWEIGHT_step1/occ_cov_temporal_8.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_normxdim_anticausal/occ_cov_value.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_normxdim/occ_cc_value.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_prearti/occ_cc_value.mat',]);
        Subj = ['NY',subjs{sub}];
        % get paths to all visualization files
        paths = get_path(root_dir,Subj);
        channel_info_all = get_channel_info(paths);
        plot_data_all_mni = get_plot_data(paths, 'mni', 'lh');
%         plot_data_all_T1 = get_plot_data(paths, 'subj', 'lh');
    %     coord = T1;
        coord = mni;
        aud = isaud(regions);
        mask_nonaud = mask & ~aud;
%         onset_diff = (onset-continuouse);
        if strcmp(attr,'freq_formants_hamon')
            att = eval(attr);
            data = att(:,on:off,1,:,:);
            if use_occ
                Att1 = gather_att2(data,postabs,max_att);
                Att1_temp = Att1;
            else
                [Att1,Att1_temp] = gather_att(data,postabs,max_att,dsrate);
            end
            sorted = sort(Att1(:));
            % get all plot frames
            if visdiff && ~post
%                 Att1_post = gather_att(att(ind,on_post:off_post,1,:,:),postabs,max_att);
%                 coef = find_normcoef(Att1_post,Att1,mask,max_anker,ankercount);
%                 diff = Att1_post*coef(1)+coef(2) - Att1;
%                 sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                 diff = diff/maxv;
%                 diff = (diff+1)/2;
%                 [frames] = VisualAtt(diff, mask, coord, plot_data_all_mni,color_map_div,1);
                [Att1_post,Att1_post_temp] = gather_att(att(:,on_post:off_post,1,:,:),postabs,max_att,dsrate);
                if wreg
                    [Att1_post,Att1_post_temp,m25,m75] = robust_rescale(Att1_post,Att1_post_temp,mask,true,regions);
                    [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,true,regions);
                else
                    [Att1_post,Att1_post_temp,m25,m75] = robust_rescale(Att1_post,Att1_post_temp,mask,true);
                    [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,true);
                end
                diff = Att1_post-Att1;
                [diff] = checkaud(diff,mask,regions);
                sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
                diff = diff/maxv;
                diff = (diff+1)/2;
                if temp
                    diff_cell{sub} = ((Att1_post_temp-Att1_temp)/maxv+1)/2;
                else
                    diff_cell{sub} = diff;
                end
                mask_cell{sub} = mask;
                coord_cell{sub} = coord;
                if plot_on_same_brain 
                    if sub == SUB
                        [frames] = VisualAtt(diff_cell, mask_cell, region_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true,annotation);
%                         axes(haxis(sub2ind([col,row],dumm,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
                    end
                else
                    [frames] = VisualAtt({diff}, {mask}, {regions},{coord}, plot_data_all_mni,color_map_div,1,false,true,annotation);
%                     axes(haxis(sub2ind([col,row],dumm,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
                end
            else
%                 [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map,sorted(end-2));
                if wreg
                    [Att1_md,Att1_md_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,true,regions);
                    [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median,regions);
                else
                    [Att1_md,Att1_md_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask);
                    [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median);
                end
                if temp
                    Att1_cell{sub} = Att1_temp;
                else
                    Att1_cell{sub} = Att1;
                end
                mask_cell{sub} = mask;
                coord_cell{sub} = coord;
                if plot_on_same_brain 
                    if sub == SUB
                        [frames] = VisualAtt(Att1_cell, mask_cell, region_cell,coord_cell, plot_data_all_mni,clrmap,1,true,false,annotation);
%                         axes(haxis(sub2ind([col,row],dumm,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
                    end
                else
                    [frames] = VisualAtt({Att1}, {mask}, {regions},{coord}, plot_data_all_mni,color_map,1,false,false,annotation);
%                     axes(haxis(sub2ind([col,row],dumm,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
                end
            end
            
            
            data = att(:,on:off,2,:,:);
            if use_occ
                Att2 = gather_att2(data,postabs,max_att);
                Att2_temp = Att2;
            else
                [Att2,Att2_temp] = gather_att(data,postabs,max_att,dsrate);
            end
%             [Att2,Att2_temp] = gather_att(data,postabs,max_att,dsrate);
            sorted = sort(Att2(:));
            if visdiff && ~post
%                 Att2_post = gather_att(att(ind,on_post:off_post,2,:,:),postabs,max_att);
%                 coef = find_normcoef(Att2_post,Att2,mask,max_anker,ankercount);
%                 diff = Att2_post*coef(1)+coef(2) - Att2;
%                 sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                 diff = diff/maxv;
%                 diff = (diff+1)/2;
%                 [frames] = VisualAtt(diff, mask, coord, plot_data_all_mni,color_map_div,1);
                [Att2_post,Att2_post_temp] = gather_att(att(:,on_post:off_post,1,:,:),postabs,max_att,dsrate);
                if wreg
                    [Att2_post,Att2_post_temp,m25,m75] = robust_rescale(Att2_post,Att2_post_temp,mask,true,regions);
                    [Att2,Att2_temp,m25,m75] = robust_rescale(Att2,Att2_temp,mask,true,regions);
                else
                    [Att2_post,Att2_post_temp,m25,m75] = robust_rescale(Att2_post,Att2_post_temp,mask);
                    [Att2,Att2_temp,m25,m75] = robust_rescale(Att2,Att2_temp,mask);
                end
                diff = Att2_post-Att2;
                [diff] = checkaud(diff,mask,regions);
                sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
                diff = diff/maxv;
                diff = (diff+1)/2;
                if temp
                    diff_cell{sub} = ((Att2_post_temp-Att2_temp)/maxv+1)/2;
                else
                    diff_cell{sub} = diff;
                end
                mask_cell{sub} = mask;
                coord_cell{sub} = coord;
                if plot_on_same_brain 
                    if sub == SUB
                        [frames2] = VisualAtt(diff_cell, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map_div,1,true,true,annotation);
%                         axes(haxis(sub2ind([col,row],dumm+1,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
                    end
                else
                    [frames2] = VisualAtt({diff}, {mask},{regions}, {coord}, plot_data_all_mni,color_map_div,1,false,true,annotation);
%                     axes(haxis(sub2ind([col,row],dumm+1,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
                end
            else
%                 [frames] = VisualAtt(Att2, mask, coord, plot_data_all_mni,color_map,sorted(end-2));
                if wreg
                    [Att2_md,Att2_md_temp,m25,m75] = robust_rescale(Att2,Att2_temp,mask,true,regions);
                    [Att2,Att2_temp,m25,m75] = robust_rescale(Att2,Att2_temp,mask,above_median,regions);
                else
                    [Att2_md,Att2_md_temp,m25,m75] = robust_rescale(Att2,Att2_temp,mask);
                    [Att2,Att2_temp,m25,m75] = robust_rescale(Att2,Att2_temp,mask,above_median);
                end
                if temp
                    Att2_cell{sub} = Att2_temp;
                else
                    Att2_cell{sub} = Att2;
                end
                Att1 = (Att1+Att2)/2;
                Att1_cell{sub} = Att1;
                mask_cell{sub} = mask;
                coord_cell{sub} = coord;
                if plot_on_same_brain 
                    if sub == SUB
                        [frames2] = VisualAtt(Att2_cell, mask_cell,region_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false,annotation);
%                         axes(haxis(sub2ind([col,row],dumm+1,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
                    end
                else
                    [frames2] = VisualAtt({Att2}, {mask},{regions}, {coord}, plot_data_all_mni,color_map,1,false,false,annotation);
%                     axes(haxis(sub2ind([col,row],dumm+1,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
                end
            end

            
%             coef = find_normcoef(Att2,Att1,mask,false,20);
%             diff = Att2*coef(1)+coef(2) - Att1;
%             [Att1,Att1_temp] = gather_att(att(ind,on:off,1,:,:),postabs,max_att,dsrate);
%             [Att2,Att2_temp] = gather_att(att(ind,on:off,2,:,:),postabs,max_att,dsrate);
            diff = Att2-Att1;
%             sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
            ind_active = find(mask==1);
            sorted = sort(diff(ind_active)); maxv = max(sorted(end-1),abs(sorted(2)));%maxv = min(sorted(end-1),abs(sorted(2)));
            diff = diff/maxv;
            diff = (diff+1)/2;
            if temp
                diff2_cell{sub} = ((Att2_md_temp-Att1_md_temp)/maxv+1)/2;
            else
                diff2_cell{sub} = diff;
            end
            mask_cell{sub} = mask;
            coord_cell{sub} = coord;
            if plot_on_same_brain 
                if sub == SUB
                    [frames3] = VisualAtt(diff2_cell, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map_div,1,true,true,annotation);
%                     axes(haxis(sub2ind([col,row],dumm+2,1)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
                end
            else
                [frames3] = VisualAtt({diff}, {mask}, {regions},{coord}, plot_data_all_mni,color_map_div,1,false,true,annotation);
%                 axes(haxis(sub2ind([col,row],dumm+2,sub)));
%                 imshow(uint8(squeeze(frames)));
%                 axis off;
            end
            if sub == SUB
                frames = cat(4,frames,frames2,frames3);
            end
            
        else
            if strcmp(attr,'freq_harm_diff')
%                 att = eval([attr,'_']);
%                 data = att(1:50,1:128,:,:,:);
% %                 data = att(ind,on:off,:,:,:);
%                 Att1_temp = data;
%                 Att1 = mean(mean(squeeze(data),1),2);
%                 maxv = min(max(Att1(:)),abs(min(Att1(:))));
%                 Att1 = Att1/maxv;
%                 Att1 = (Att1+1)/2;
                % get all plot frames
%                 [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map_div,1);
                att = eval(attr);
                [Att1,Att1_temp] = gather_att(att(:,on:off,:,:,:),postabs,max_att,dsrate);
                if wreg
                    [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,true,regions);
                else
                    [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,true);
                end
                if temp
                    Att1_cell{sub} = Att1_temp;
                else
                    Att1_cell{sub} = Att1;
                end
                mask_cell{sub} = mask;
                coord_cell{sub} = coord;
                if plot_on_same_brain 
                    if sub == SUB
                        [frames] = VisualAtt(Att1_cell, mask_cell,region_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false,annotation);
%                         axes(haxis(sub2ind([col,row],dumm,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
                    end
                else
                    [frames] = VisualAtt({Att1}, {mask},{regions}, {coord}, plot_data_all_mni,clrmap,1,false,false,annotation);
%                     axes(haxis(sub2ind([col,row],dumm,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
                end
            else
                if strcmp(attr,'onset')
                    att = onset_diff;
                else
                    if ~strcmp(attr,'avg')
                        att = eval(attr);
                    end
                end
                if ~strcmp(attr,'avg')
                    data = att(ind,on:off,1,:,:);
                    sorted = sort(data(:));
                    sorted = sorted(end-1);
                    if use_occ && (strcmp(attr,'f0_hz') || strcmp(attr,'amplitudes'))
                        Att1 = gather_att2(data,postabs,max_att);
                        Att1_temp = Att1;
                    else
                        if strcmp(attr,'ecog')
                            [Att1,Att1_temp] = gather_att3(data,postabs,max_att,dsrate);
                            
                        else
                            [Att1,Att1_temp] = gather_att(data,postabs,max_att,dsrate);
                            
                        end
                    end
                    sorted = sort(Att1(:));
    %                 sorted = sorted(end-1);
                end
                % get all plot frames
                if strcmp(attr,'loudness') || strcmp(attr,'f0_hz')
                    if visdiff && ~post
%                         Att1_post = gather_att(att(:,on_post:off_post,1,:,:),postabs,max_att);
%                         coef = find_normcoef(Att1_post,Att1,mask,max_anker,ankercount);
%                         diff = Att1_post*coef(1)+coef(2) - Att1;
%                         sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                         diff = diff/maxv;
%                         diff = (diff+1)/2;
%                         [frames] = VisualAtt(diff, mask, coord, plot_data_all_mni,color_map_div,1);
                        
                        [Att1_post,Att1_post_temp] = gather_att(att(:,on_post:off_post,1,:,:),postabs,max_att,dsrate);

                        if wreg
                            [Att1_post,Att1_post_temp,m25,m75] = robust_rescale(Att1_post,Att1_post_temp,mask,true,regions);
                            [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,true,regions);
                        else
                            [Att1_post,Att1_post_temp,m25,m75] = robust_rescale(Att1_post,Att1_post_temp,mask);
                            [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask);
                        end
                        diff = Att1_post-Att1;
                        [diff] = checkaud(diff,mask,regions);
                        sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
                        diff = diff/maxv;
                        diff = (diff+1)/2;
                        if temp
                            diff_cell{sub} = ((Att1_post_temp-Att1_temp)/maxv+1)/2;
                        else
                            diff_cell{sub} = diff;
                        end
                        mask_cell{sub} = mask;
                        coord_cell{sub} = coord;
                        if plot_on_same_brain 
                            if sub == SUB
                                [frames] = VisualAtt(diff_cell, mask_cell, region_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true,annotation);
%                                 axes(haxis(sub2ind([col,row],dumm,1)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
                            end
                        else
                            [frames] = VisualAtt({diff}, {mask}, {regions},{coord}, plot_data_all_mni,color_map_div,1,false,true,annotation);
%                             axes(haxis(sub2ind([col,row],dumm,sub)));
%                             imshow(uint8(squeeze(frames)));
%                             axis off;
                        end
                    else
%                         [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map,maxminp1);
                        if wreg
                            [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median,regions);
                        else
                            [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median);
                        end
                        if temp
                            Att1_cell{sub} = Att1_temp;
                        else
                            Att1_cell{sub} = Att1;
                        end
                        mask_cell{sub} = mask;
                        coord_cell{sub} = coord;
                        if plot_on_same_brain 
                            if sub == SUB
                                [frames] = VisualAtt(Att1_cell, mask_cell, region_cell,coord_cell, plot_data_all_mni,clrmap,1,true,false,annotation);
%                                 axes(haxis(sub2ind([col,row],dumm,1)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
                            end
                        else
                            [frames] = VisualAtt({Att1}, {mask}, {regions},{coord}, plot_data_all_mni,color_map,1,false,false,annotation);
%                             axes(haxis(sub2ind([col,row],dumm,sub)));
%                             imshow(uint8(squeeze(frames)));
%                             axis off;
                        end
                    end
                else
                    if strcmp(attr,'amplitudes') || strcmp(attr,'onset')
                        if visdiff && ~post
%                             Att1_post = gather_att(att(:,on_post:off_post,1,:,:),postabs,max_att);
%                             coef = find_normcoef(Att1_post,Att1,mask,max_anker,ankercount);
%                             diff = Att1_post*coef(1)+coef(2) - Att1;
%                             sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                             diff = diff/maxv;
%                             diff = (diff+1)/2;
%                             [frames] = VisualAtt(diff, mask, coord, plot_data_all_mni,color_map_div,1);
                            [Att1_post,Att1_post_temp] = gather_att(att(:,on_post:off_post,1,:,:),postabs,max_att,dsrate);
                            if wreg
                                [Att1_post,Att1_post_temp,m25,m75] = robust_rescale(Att1_post,Att1_post_temp,mask,true,regions);
                                [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,true,regions);
                            else
                                [Att1_post,Att1_post_temp,m25,m75] = robust_rescale(Att1_post,Att1_post_temp,mask);
                                [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask);
                            end
                            diff = Att1_post-Att1;
                            [diff] = checkaud(diff,mask,regions);
                            sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
                            diff = diff/maxv;
                            diff = (diff+1)/2;
                            if temp
                                diff_cell{sub} = ((Att1_post_temp-Att1_temp)/maxv+1)/2;
                            else
                                diff_cell{sub} = diff;
                            end
                            mask_cell{sub} = mask;
                            coord_cell{sub} = coord;
                            if plot_on_same_brain 
                                if sub == SUB
                                    [frames] = VisualAtt(diff_cell, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map_div,1,true,true,annotation);
%                                     axes(haxis(sub2ind([col,row],dumm,1)));
%                                     imshow(uint8(squeeze(frames)));
%                                     axis off;
                                end
                            else
                                [frames] = VisualAtt({diff}, {mask}, {regions},{coord}, plot_data_all_mni,color_map_div,1,false,true,annotation);
%                                 axes(haxis(sub2ind([col,row],dumm,sub)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
                            end
                        else
%                             [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map,sorted(end-1));
                            att = eval([attr,'_']);
%                             data = att(1:50,1:128,:,:,:);
                            if entire_period
                                on_ = 1;
                                off_ = 128;
                            else
                                if post
                                    on_ = (16+32)+1;
                                    off_ = 120;
                                else
                                    on_ = 16+1;
                                    off_ = (16+32);
                                    on_post_ = (16+32)+1;
                                    off_post_ = 120;
                                end
                            end
                            if size(att,2)<36
                                if temp
                                    data = att(:,on:off,:,:,:);
                                else
                                    data = att(:,on:off,:,:,:);
                                end
                            else
                                if temp
                                    data = att(:,on_:off_,:,:,:);
                                else
                                    data = att(:,on_:off_,:,:,:);
                                end
                            end
                            Att1 = mean(mean(squeeze(data),1),2);
                            sorted = sort(Att1(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                             ind_active = find(mask==1);
%                             sorted = sort(Att1(ind_active));maxv = min(sorted(end-1),abs(sorted(2)));
                            Att1 = Att1/maxv;
                            Att1 = (Att1+1)/2;
                            if temp
                                data_temp = reshape(data,size(data,1),dsrate,size(data,2)/dsrate,1,size(data,4),size(data,5));
                                Att1_temp = mean(squeeze(mean(data_temp,2)),1);
                                sizes = size(Att1_temp);
%                                 Att1_temp = reshape(Att1_temp,[sizes(1),sizes(3:end)]);
                                sorted = sort(Att1_temp(:)); maxv = min(sorted(int32(0.99*length(sorted))),abs(sorted(int32(0.01*length(sorted))+1)));
                                Att1_temp = (Att1_temp/maxv+1)/2;
                                Att1_cell{sub} = Att1_temp;
                            else
                                Att1_cell{sub} = Att1;
                            end
                            mask_cell{sub} = mask;
                            coord_cell{sub} = coord;
                            if plot_on_same_brain
                                if sub == SUB
                                    [frames] = VisualAtt(Att1_cell, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map_div,1,true,true,annotation);
%                                     axes(haxis(sub2ind([col,row],dumm,1)));
%                                     imshow(uint8(squeeze(frames)));
%                                     axis off;
                                    for subb =1:length(Att1_cell)
                                        Att1_cell{subb}=0;
                                    end
                                end
                            else
                                [frames] = VisualAtt({Att1}, {mask}, {regions},{coord}, plot_data_all_mni,color_map_div,1,false,true,annotation);
%                                 axes(haxis(sub2ind([col,row],dumm,sub)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
                            end

                        end
                    else
                        if strcmp(attr,'avg')
                            mask_cell{sub} = mask;
                            coord_cell{sub} = coord;
                            if sub == SUB
                                [frames] = VisualAtt(avg_data, mask_cell, region_cell,coord_cell, plot_data_all_mni,clrmap,1,true,false,annotation);
                            end
                        else
                            if visdiff && ~post
    %                             Att1_post = gather_att(att(:,on_post:off_post,1,:,:),postabs,max_att);
    %                             coef = find_normcoef(Att1_post,Att1,mask,max_anker,ankercount);
    %                             diff = Att1_post*coef(1)+coef(2) - Att1;
    %                             sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
    %                             diff = diff/maxv;
    %                             diff = (diff+1)/2;
    %                             [frames] = VisualAtt(diff, mask, coord, plot_data_all_mni,color_map_div,1);
                                [Att1_post,Att1_post_temp] = gather_att(att(:,on_post:off_post,1,:,:),postabs,max_att,dsrate);
                                if wreg
                                    [Att1_post,Att1_post_temp,m25,m75] = robust_rescale(Att1_post,Att1_post_temp,mask,true,regions);
                                    [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,true,regions);
                                else
                                    [Att1_post,Att1_post_temp,m25,m75] = robust_rescale(Att1_post,Att1_post_temp,mask);
                                    [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask);
                                end
                                diff = Att1_post-Att1;
                                [diff] = checkaud(diff,mask,regions);
                                sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
                                diff = diff/maxv;
                                diff = (diff+1)/2;
                                if temp
                                    diff_cell{sub} = ((Att1_post_temp-Att1_temp)/maxv+1)/2;
                                else
                                    diff_cell{sub} = diff;
                                end

                                mask_cell{sub} = mask;
                                coord_cell{sub} = coord;
                                if plot_on_same_brain
                                    if sub == SUB
                                        [frames] = VisualAtt(diff_cell, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map_div,1,true,true,annotation);
    %                                     axes(haxis(sub2ind([col,row],dumm,1)));
    %                                     imshow(uint8(squeeze(frames)));
    %                                     axis off;
                                    end
                                else
                                    [frames] = VisualAtt({diff}, {mask}, {regions}, {coord}, plot_data_all_mni,color_map_div,1,false,true,annotation);
    %                                 axes(haxis(sub2ind([col,row],dumm,sub)));
    %                                 imshow(uint8(squeeze(frames)));
    %                                 axis off;
                                end
                            else
    %                             [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map);
                                if wreg
                                    [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median,regions);
                                else
                                    [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median);
                                    Att1 = max(Att1,0).^0.43/1.25;
                                end
                                if temp
                                    Att1_cell{sub} = Att1_temp;
                                else
                                    Att1_cell{sub} = Att1;
                                end
                                mask_cell{sub} = mask;
                                coord_cell{sub} = coord;
                                if plot_on_same_brain 
                                    if sub == SUB
                                        [frames] = VisualAtt(Att1_cell, mask_cell, region_cell,coord_cell, plot_data_all_mni,clrmap,1,true,false,annotation);
    %                                     axes(haxis(sub2ind([col,row],dumm,1)));
    %                                     imshow(uint8(squeeze(frames)));
    %                                     axis off;
                                    end
                                else
                                    [frames] = VisualAtt({Att1}, {mask}, {regions}, {coord}, plot_data_all_mni,color_map,1,false,false,annotation);
    %                                 axes(haxis(sub2ind([col,row],dumm,sub)));
    %                                 imshow(uint8(squeeze(frames)));
    %                                 axis off;
                                end
                            end
                        end
                    end
                end
            end
        end
%         if ~plot_on_same_brain
%             if dumm == 2
%                 ylabel(Subj)
%             end  
%         end
%         if plot_on_same_brain
%             if sub == SUB
%                 title(attr)
%             end
%         else
%             if sub ==1
%                 title(attr)
%             end
%         end
    end
%     if strcmp(attr,'freq_formants_hamon')
%         dumm = dumm+3;
%     else
%         dumm = dumm+1;
%     end
end

function [frames1,frames2,frames,Att1_cell] = attr_vis_contrast(attr,post,visdiff,postabs,max_att,plot_on_same_brain,above_median,temp,wreg,entire_period,avg_data,annotation,bar_plot,norm_contr_elec)
%     maxminp1=[];
%     maxminp2=[10000];
%     for sub =1:SUB
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/NY',subjs{sub},'_elec.mat']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/att_NY',subjs{sub},'']);
%         att = eval(attrs{a});
%         if strcmp(attrs{a},'freq_formants_hamon')
%             maxp1 = gather_att(att(ind,on:off,1,:,:),postabs,max_att);
%             maxp2 = gather_att(att(ind,on:off,2,:,:),postabs,max_att);
%             maxminp1 = [maxminp1,max(maxp1(:))];
%             maxminp2 = [maxminp2,max(maxp2(:))];
%         else
%             maxp1 = gather_att(att(ind,on:off,1,:,:),postabs,max_att);
%             maxminp1 = [maxminp1,max(maxp1(:))];
%         end
%     end
%     maxminp1 = min(maxminp1);
%     maxminp2 = min(maxminp2);
    % choose a color map and set indecies of colors for each electrode

    color_map = hot(256);%afmhot;%hot;
    color_map = color_map(end:-1:1,:);
%     [color_map]=cbrewer('seq', 'Blues', 256);
%     [color_map_div]=cbrewer('div', 'RdBu', 256,'PCHIP');
    [color_map_div]=cbrewer('div', 'PRGn', 256,'PCHIP');
    [color_map_cb]=cbrewer('seq', 'Greens', 256,'PCHIP');
%     color_map_cb = color_map_cb(end:-1:1,:);
%     color_map = color_map_cb;
    [colorset1] = cbrewer('qual', 'Set1', 9,'PCHIP');
    % color_map_div = color_map_div(end:-1:1,:);
    color_ind = value2colorind(ones(128,1),'hot',[0,1]); 
    color_ind = [color_ind,color_ind,color_ind];
    colors = {color_ind, color_map};
    clrmap = color_map;
    if strcmp(attr,'ecog')
        dsrate_org = 1;
    else
        dsrate_org = 8;
    end
    if entire_period
        on = 16/dsrate_org+1;
        off = 120/dsrate_org;
    else
        if post
            on = (16+32)/dsrate_org+1;
            off = 120/dsrate_org;
        else
            on = 16/dsrate_org+1;
            off = (16+32)/dsrate_org;
            on_post = (16+32)/dsrate_org+1;
            off_post = 120/dsrate_org;
        end
    end

    dsrate = 8/dsrate_org;
    % ind = [1:10,21:30];
    ind = [1:2];%[1:50];
    dumm = 1;
    Att1_cell = {};
    Att2_cell = {};
    diff_cell = {};
    diff2_cell = {};
    att_cell1 = {};
    att_cell2 = {};
    mask_cell = {};
    coord_cell = {};
    region_cell = {};
    subjs = {'717','742','749'};
%     subjs = {'798'};
    cmax=1.2;
    cax=[-0.7,0.7];
%     subjs = {'749'};
    SUB = length(subjs);
    root_dir = '/Users/ranwang/Documents/writen_paper/NER2020/Visuallization_matlab';
    use_occ = false;%(strcmp(attr,'amplitudes') || strcmp(attr,'f0_hz') || strcmp(attr,'freq_formants_hamon')) && ~temp && ~visdiff && post;
    use_occentire = false;%strcmp(attr,'freq_formants_hamon') && ~temp && ~visdiff && post;
    attrs = {'spec'};
%     attrs = {'spec','loudness','f0_hz','amplitudes','freq_formants_hamon','freq_formants_noise','bandwidth_formants_noise_hz'};
%     contrast = 'sin';%'sqrt';%'minus';
    contrast = 'sqrt';
    the_power = 0.43;%0.43;%0.75;%0.33;%0.39;%0.36;%0.43;
    for sub=1:SUB
        
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/att_NY',subjs{sub},'']);
%         if use_occentire
%             load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_cc_value_',subjs{sub},'']); 
%             ind = [1:2];
%         else
%             ind = [1:50];
%             if use_occ
%                 load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_dict_IG_value_',subjs{sub},'']);   
%     %             load(['/Users/ranwang/Documents/writen_paper/NER2020/attr_dict_IG_value_',subjs{sub},'_nonnoise_']);
%             else
%                 load(['/Users/ranwang/Documents/writen_paper/NER2020/attr_dict_IG_value_',subjs{sub},'_nonnoise_']);
%     %             load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_dict_IG_value_',subjs{sub},'']);
%             end
%         end
        load(['/Users/ranwang/Documents/writen_paper/NER2020/NY_',subjs{sub},'_elec_entiregrid.mat']);
%         onregion = {'cSTG','mSTG','parstriangularis','parsopercularis','precentral','postcentral','inferiorparietal','supramarginal'};
        onregion = {'cSTG','mSTG','parstriangularis','parsopercularis','precentral','postcentral'};
%         onregion = {'inferiorparietal','supramarginal'};
%         mask = isregion(regions,onregion);
        region_cell{sub} = regions;
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_woaud/','NY',subjs{sub},'_elec_woaud.mat']);

%         load(['/Users/ranwang/Documents/writen_paper/NER2020/ecog_',subjs{sub},'prod.mat']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_dict_IG_value_',subjs{sub},'']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_dict_IG_value_',subjs{sub},'_prearti']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_cc_value_prearti_',subjs{sub}]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_/occ_cc_temporal_8.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_/occ_cc_value.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_wpariental/occ_cc_value.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_wpariental/occ_cc_temporal_8.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/entiregrid_',subjs{sub},'_normxdim_causal_covnorm_NOVISWEIGHT_step1/occ_cov_temporal_8.mat',]);
%         load(['/Users/ranwang/greene/neural_decoding/code/cnn/ALAE/training_artifacts/entiregrid_',subjs{sub},'_han5amppowerloss_ecogfinetune_alphasup3_groupnormxdim_causal/grad_cov_temporal_8.mat',]);
        load(['/Users/ranwang/greene/neural_decoding/code/cnn/ALAE/training_artifacts/entiregrid_',subjs{sub},'_han5amppowerloss_alphasup3_percept_singlebandwidth_fractive1500lim_bgnoiselearn_ecogfinetune_anticausal_universal_passive_groupnormxdim_correctbaseline/grad_cov_temporal_8.mat',]);
%         load(['/Users/ranwang/greene/neural_decoding/code/cnn/ALAE/training_artifacts/entiregrid_',subjs{sub},'_han5amppowerloss_alphasup3_percept_singlebandwidth_fractive1500lim_bgnoiselearn_ecogfinetune_anticausal_passive_groupnormxdim_correctbaseline/grad_cov_temporal_8.mat',]); % this is in green; frame2
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/entiregrid_ecog_',subjs{sub},'percept_imagine.mat']);
        load(['/Users/ranwang/Documents/writen_paper/NER2020/entiregrid_ecog_',subjs{sub},'percept_passive.mat']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/entiregrid_',subjs{sub},'_normxdim_causal_covnorm_NOVISWEIGHT_step1_grad/grad_cov_temporal_8.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_normxdim_anticausal/occ_cov_value.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_normxdim/occ_cc_value.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_prearti/occ_cc_value.mat',]);
        for a= 1:length(attrs)
            eval([attrs{a},'_a=',attrs{a},';']);
            try
                eval([attrs{a},'_a_=',attrs{a},'_;']);
            catch
                eval([attrs{a},'_a_=',attrs{a},';']);
            end
        end
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/entiregrid_',subjs{sub},'_normxdim_anticausal_covnorm_NOVISWEIGHT_step1/occ_cov_temporal_8.mat',]);
        load(['/Users/ranwang/greene/neural_decoding/code/cnn/ALAE/training_artifacts/entiregrid_',subjs{sub},'_han5amppowerloss_ecogfinetune_alphasup3_groupnormxdim_anticausal/grad_cov_temporal_8.mat',]);
        load(['/Users/ranwang/Documents/writen_paper/NER2020/entiregrid_ecog_',subjs{sub},'percept_active.mat']);
        % load(['/Users/ranwang/Documents/writen_paper/NER2020/entiregrid_ecog_',subjs{sub},'prod.mat']);
%         load(['/Users/ranwang/greene/neural_decoding/code/cnn/ALAE/training_artifacts/entiregrid_',subjs{sub},'_han5amppowerloss_alphasup3_percept_singlebandwidth_fractive1500lim_bgnoiselearn_ecogfinetune_anticausal_universal_active_groupnormxdim_correctbaseline/grad_cov_temporal_8.mat',]); % this is in purple; frame1
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/entiregrid_',subjs{sub},'_normxdim_anticausal_covnorm_NOVISWEIGHT_step1_grad/grad_cov_temporal_8.mat',]);
        for a= 1:length(attrs)
            eval([attrs{a},'=',attrs{a},';']);
            try
                eval([attrs{a},'_=',attrs{a},'_;']);
            catch
                eval([attrs{a},'_=',attrs{a},';']);
            end
        end
        Subj = ['NY',subjs{sub}];
        % get paths to all visualization files
        paths = get_path(root_dir,Subj);
        channel_info_all = get_channel_info(paths);
        plot_data_all_mni = get_plot_data(paths, 'mni', 'lh');
%         plot_data_all_T1 = get_plot_data(paths, 'subj', 'lh');
    %     coord = T1;
        coord = mni;
        aud = isaud(regions);
        mask_nonaud = mask & ~aud;
%         onset_diff = (onset-continuouse);
        if strcmp(attr,'freq_formants_hamon')
            att = eval(attr);
            data = att(:,on:off,1,:,:);
            att_a = eval([attr,'_a']);
            data_a = att_a(:,on:off,1,:,:);
            if use_occ
                Att1 = gather_att2(data,postabs,max_att);
                Att1_temp = Att1;
                Att1_a = gather_att2(data_a,postabs,max_att);
                Att1_temp_a = Att1_a;
            else
                [Att1,Att1_temp] = gather_att(data,postabs,max_att,dsrate);
                [Att1_a,Att1_temp_a] = gather_att(data_a,postabs,max_att,dsrate);
                [Att1_causal,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median);
%                 [Att1_anticausal,Att1_temp_a,m25,m75] = robust_rescale(Att1_a,Att1_temp_a,mask,above_median);
%                 Att1 = Att1_causal;
%                 Att1_a = Att1_anticausal;
                
                Att1_anticausal = (Att1_a-m25)/(m75-m25);
                Att1_temp_anticausal = (Att1_temp_a-m25)/(m75-m25); 
%                 Att1_causal_ = Att1_causal./(abs(Att1_causal)+abs(Att1_anticausal));
% %                 Att1_anticausal = Att1_anticausal./(abs(Att1_causal)+abs(Att1_anticausal));
%                 Att1_causal = Att1_causal_;
%                 Att1 = Att1_causal;
%                 Att1_a = Att1_anticausal;

                Att1_causal = (max(Att1_causal,0)).^the_power;
                Att1_anticausal = (max(Att1_anticausal,0)).^the_power;
%                 Att1 = sign(Att1-Att1_a).*(max(Att1-Att1_a,0)).^the_power;%sign(Att1-Att1_a).*sqrt(abs(Att1-Att1_a));
                Att1 = (Att1-Att1_a)./(abs(Att1)+abs(Att1_a));
                [Att1] = checkaud(Att1,mask,regions);
                sorted = sort(Att1(:)); maxv = max(sorted(end-1),abs(sorted(2)));%maxv = min(sorted(end-1),abs(sorted(2)));
%                 Att1 = Att1/maxv;
%                 Att1 = (Att1+1)/2;
            end
            diff_cell{sub} = Att1;
            att_cell1{sub} = Att1_causal;
            att_cell2{sub} = Att1_anticausal;
            mask_cell{sub} = mask;
            coord_cell{sub} = coord;
            if plot_on_same_brain 
                if sub == SUB
                    [frames11,prob1,prob0] = VisualAtt(att_cell1, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map,cmax,true,false,annotation,bar_plot);
                    [frames12,prob2,prob0] = VisualAtt(att_cell2, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map,cmax,true,false,annotation,bar_plot);
                    if bar_plot
                        [frames,prob1,prob0] = VisualAtt(diff_cell, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map_div,1,true,true,annotation,bar_plot);
                    else
                        switch contrast
                            case 'sin'
                                da = sign(prob1'-prob2').*sqrt((prob1'-prob2').^2)./sqrt((prob1'.^2+prob2'.^2)*2);
                                da(prob0<1)=0;
                                frames={};Lcrtx = load('ch2_template_lh_pial.mat'); Lcrtx.tri = Lcrtx.faces;Lcrtx.vert = Lcrtx.coords; frames{1}=tripatch(Lcrtx, 'nofigure', da); shading interp;colormap(gca,color_map_div);material dull;lighting gouraud;light;axis off;viewangle = 'l';litebrain(viewangle,.8);set(gca,'CLim',[cax(1) cax(2)]); 
                            case 'sqrt'
                                frames={};Lcrtx = load('ch2_template_lh_pial.mat'); Lcrtx.tri = Lcrtx.faces;Lcrtx.vert = Lcrtx.coords; frames{1}=tripatch(Lcrtx, 'nofigure', sign(prob1'-prob2').*(abs(prob1'-prob2')).^the_power); shading interp;colormap(gca,color_map_div);material dull;lighting gouraud;light;axis off;viewangle = 'l';litebrain(viewangle,.8);;set(gca,'CLim',[cax(1) cax(2)]); 
                            case 'minu'
                                frames={};Lcrtx = load('ch2_template_lh_pial.mat'); Lcrtx.tri = Lcrtx.faces;Lcrtx.vert = Lcrtx.coords; frames{1}=tripatch(Lcrtx, 'nofigure', prob1'-prob2'); shading interp;colormap(gca,color_map_div);material dull;lighting gouraud;light;axis off;viewangle = 'l';litebrain(viewangle,.8);;set(gca,'CLim',[cax(1) cax(2)]); 
                        end
                    end
                    
                    
%                         axes(haxis(sub2ind([col,row],dumm+1,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
                end
            else
                [frames] = VisualAtt({diff}, {mask}, {regions},{coord}, plot_data_all_mni,color_map_div,1,false,true,annotation,bar_plot);
%                     axes(haxis(sub2ind([col,row],dumm+1,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
            end
           
            
            
            data = att(:,on:off,2,:,:);
            data_a = att_a(:,on:off,2,:,:);
            if use_occ
                Att2 = gather_att2(data,postabs,max_att);
                Att2_temp = Att2;
                Att2_a = gather_att2(data_a,postabs,max_att);
                Att2_temp_a = Att2_a;
            else
                [Att2,Att2_temp] = gather_att(data,postabs,max_att,dsrate);
                [Att2_a,Att2_temp_a] = gather_att(data_a,postabs,max_att,dsrate);
                [Att2_causal,Att2_temp,m25,m75] = robust_rescale(Att2,Att2_temp,mask,above_median);
                
%                 [Att2_anticausal,Att2_temp_a,m25,m75] = robust_rescale(Att2_a,Att2_temp_a,mask,above_median);
%                 Att2 = Att2_causal;
%                 Att2_a = Att2_anticausal;
                
                Att2_anticausal = (Att2_a-m25)/(m75-m25);
                Att2_temp_anticausal = (Att2_temp_a-m25)/(m75-m25);
                
                Att2_causal = (max(Att2_causal,0)).^the_power;
                Att2_anticausal = (max(Att2_anticausal,0)).^the_power;
%                 Att2_causal_ = Att2_causal./(abs(Att2_causal)+abs(Att2_anticausal));
%                 Att2_anticausal = Att2_anticausal./(abs(Att2_causal)+abs(Att2_anticausal));
%                 Att2_causal = Att2_causal_;
%                 Att2 = Att2_causal;
%                 Att2_a = Att2_anticausal;

%                 Att2 = (Att2-Att2_a)./(abs(Att2)+abs(Att2_a));
                Att2 = sign(Att2-Att2_a).*sqrt(abs(Att2-Att2_a));
                [Att2] = checkaud(Att2,mask,regions);
                sorted = sort(Att2(:)); maxv = max(sorted(end-1),abs(sorted(2)));%maxv = min(sorted(end-1),abs(sorted(2)));
%                 Att2 = Att2/maxv;
%                 Att2 = (Att2+1)/2;
            end
            
            diff_cell{sub} = Att2;
            att_cell1{sub} = Att2_causal;
            att_cell2{sub} = Att2_anticausal;
            mask_cell{sub} = mask;
            coord_cell{sub} = coord;
            if plot_on_same_brain 
                if sub == SUB
                    [frames21,prob1,prob0] = VisualAtt(att_cell1, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map,cmax,true,false,annotation,bar_plot);
                    [frames22,prob2,prob0] = VisualAtt(att_cell2, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map,cmax,true,false,annotation,bar_plot);
                    if bar_plot
                        [frames2,prob1,prob0] = VisualAtt(diff_cell, mask_cell, region_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true,annotation,bar_plot);
                    else
                        switch contrast
                            case 'sin'
                                da = sign(prob1'-prob2').*sqrt((prob1'-prob2').^2)./sqrt((prob1'.^2+prob2'.^2)*2);
                                da(prob0<1)=0;
                                frames2={};Lcrtx = load('ch2_template_lh_pial.mat'); Lcrtx.tri = Lcrtx.faces;Lcrtx.vert = Lcrtx.coords; frames2{1}=tripatch(Lcrtx, 'nofigure', da); shading interp;colormap(gca,color_map_div);material dull;lighting gouraud;light;axis off;viewangle = 'l';litebrain(viewangle,.8);set(gca,'CLim',[cax(1) cax(2)]); 
                            case 'sqrt'
                                frames2={};Lcrtx = load('ch2_template_lh_pial.mat'); Lcrtx.tri = Lcrtx.faces;Lcrtx.vert = Lcrtx.coords; frames2{1}=tripatch(Lcrtx, 'nofigure', sign(prob1'-prob2').*(abs(prob1'-prob2')).^the_power); shading interp;colormap(gca,color_map_div);material dull;lighting gouraud;light;axis off;viewangle = 'l';litebrain(viewangle,.8);;set(gca,'CLim',[cax(1) cax(2)]); 
                            case 'minu'
                                frames2={};Lcrtx = load('ch2_template_lh_pial.mat'); Lcrtx.tri = Lcrtx.faces;Lcrtx.vert = Lcrtx.coords; frames2{1}=tripatch(Lcrtx, 'nofigure', prob1'-prob2'); shading interp;colormap(gca,color_map_div);material dull;lighting gouraud;light;axis off;viewangle = 'l';litebrain(viewangle,.8);;set(gca,'CLim',[cax(1) cax(2)]); 
                        end
                    end
%                         axes(haxis(sub2ind([col,row],dumm+1,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
                end
            else
                [frames2] = VisualAtt({diff}, {mask}, {regions}, {coord}, plot_data_all_mni,color_map_div,1,false,true,annotation,bar_plot);
%                     axes(haxis(sub2ind([col,row],dumm+1,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
            end

            
%             coef = find_normcoef(Att2,Att1,mask,false,20);
%             diff = Att2*coef(1)+coef(2) - Att1;
%             [Att1,Att1_temp] = gather_att(att(ind,on:off,1,:,:),postabs,max_att,dsrate);
%             [Att2,Att2_temp] = gather_att(att(ind,on:off,2,:,:),postabs,max_att,dsrate);
            
            if sub == SUB
                frames = cat(4,frames,frames2);
                frames1 = cat(4,frames11,frames21);
                frames2 = cat(4,frames12,frames22);
            end
            Att1_cell{sub} = (Att1+Att2)/2;
        else
            if strcmp(attr,'whatever')
                disp('freq_harm_diff not defined')
            else
                if ~strcmp(attr,'avg')
                    att = eval(attr);
                    data = att(:,on:off,1,:,:,:);
                    att_a = eval([attr,'_a']);
                    data_a = att_a(:,on:off,1,:,:,:);
                    if use_occ
                        Att1 = gather_att2(data,postabs,max_att);
                        Att1_temp = Att1;
                        Att1_a = gather_att2(data_a,postabs,max_att);
                        Att1_temp_a = Att1_a;
                    else
                        [Att1,Att1_temp] = gather_att(data,postabs,max_att,dsrate);
                        [Att1_a,Att1_temp_a] = gather_att(data_a,postabs,max_att,dsrate);
                        [Att1_causal,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median);
%                         Att1_causal = sqrt(max(Att1_causal,0));
%                         [Att1_anticausal,Att1_temp_a,m25,m75] = robust_rescale(Att1_a,Att1_temp_a,mask,above_median);
% %                         Att1_anticausal = sqrt(max(Att1_anticausal,0));
%                         Att1 = Att1_causal;
%                         Att1_a = Att1_anticausal;
%                         m25=7.1076e-06; m75=7.0691e-04;
                        Att1_anticausal = (Att1_a-m25)/(m75-m25);
                        Att1_temp_anticausal = (Att1_temp_a-m25)/(m75-m25);
                        
                        Att1_causal = (max(Att1_causal,0)).^the_power;
                        Att1_anticausal = (max(Att1_anticausal,0)).^the_power;
                
%                         Att1_anticausal_ = Att1_anticausal./(abs(Att1_causal)+abs(Att1_anticausal));
%                         Att1_causal = Att1_causal./(abs(Att1_causal)+abs(Att1_anticausal));
%                         Att1_anticausal = Att1_anticausal_;
%                         Att1 = (Att1_causal-Att1_anticausal)./(abs(Att1_causal)+abs(Att1_anticausal));
%                         Att1 = (Att1-Att1_a)./(abs(Att1)+abs(Att1_a));
                        Att1 = (Att1-Att1_a)*500000;
                        % Att1 = (Att1-Att1_a)*1;
%                         Att1 = (Att1-Att1_a);
%                         Att1 = sign(Att1-Att1_a).*sqrt(abs(Att1-Att1_a));
%                         Att1 = sign(Att1-Att1_a).*sqrt((Att1-Att1_a).^2)./sqrt(abs(Att1).^2+abs(Att1_a).^2)/sqrt(2);
%                         [Att1] = checkaud(Att1,mask,regions);
                        sorted = sort(Att1(:)); maxv = max(sorted(end-1),abs(sorted(2)));%maxv = min(sorted(end-1),abs(sorted(2)));
    %                     Att1 = Att1/maxv;
%                         Att1 = (Att1+1)/2;
                    end
                end
%                 sorted = sorted(end-1);
                % get all plot frames
                if strcmp(attr,'something') %strcmp(attr,'loudness') || strcmp(attr,'f0_hz')

%                         [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map,maxminp1);
%                     if wreg
%                         [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median,regions);
%                     else
%                         [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median);
%                     end
                    if temp
                        Att1_cell{sub} = Att1_temp;
                    else
                        Att1_cell{sub} = Att1;
                        att_cell1{sub} = Att1_causal;
                        att_cell2{sub} = Att1_anticausal;
                    end
                    mask_cell{sub} = mask;
                    coord_cell{sub} = coord;
                    if plot_on_same_brain 
                        if sub == SUB
                            [frames1,prob1,prob0] = VisualAtt(att_cell1, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map,cmax,true,false,annotation,bar_plot);
                            [frames2,prob2,prob0] = VisualAtt(att_cell2, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map,cmax,true,false,annotation,bar_plot);
                            
                            if bar_plot 
                                [frames,prob1,prob0] = VisualAtt(Att1_cell, mask_cell,region_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true,annotation,bar_plot);
                            else
                                switch contrast
                                    case 'sin'
                                        da = sign(prob1'-prob2').*sqrt((prob1'-prob2').^2)./sqrt((prob1'.^2+prob2'.^2)*2);
                                        da(prob0<1)=0;
                                        frames={};Lcrtx = load('ch2_template_lh_pial.mat'); Lcrtx.tri = Lcrtx.faces;Lcrtx.vert = Lcrtx.coords; frames{1}=tripatch(Lcrtx, 'nofigure', da); shading interp;colormap(gca,color_map_div);material dull;lighting gouraud;light;axis off;viewangle = 'l';litebrain(viewangle,.8);set(gca,'CLim',[cax(1) cax(2)]); 
                                    case 'sqrt'
                                        frames={};Lcrtx = load('ch2_template_lh_pial.mat'); Lcrtx.tri = Lcrtx.faces;Lcrtx.vert = Lcrtx.coords; frames{1}=tripatch(Lcrtx, 'nofigure', sign(prob1'-prob2').*(abs(prob1'-prob2')).^the_power); shading interp;colormap(gca,color_map_div);material dull;lighting gouraud;light;axis off;viewangle = 'l';litebrain(viewangle,.8);;set(gca,'CLim',[cax(1) cax(2)]); 
                                    case 'minu'
                                        frames={};Lcrtx = load('ch2_template_lh_pial.mat'); Lcrtx.tri = Lcrtx.faces;Lcrtx.vert = Lcrtx.coords; frames{1}=tripatch(Lcrtx, 'nofigure', prob1'-prob2'); shading interp;colormap(gca,color_map_div);material dull;lighting gouraud;light;axis off;viewangle = 'l';litebrain(viewangle,.8);;set(gca,'CLim',[cax(1) cax(2)]); 
                                end
                            end
%                                 axes(haxis(sub2ind([col,row],dumm,1)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
                        end
                    else
                        [frames] = VisualAtt({Att1}, {mask}, {regions},{coord}, plot_data_all_mni,color_map,1,false,false,annotation,bar_plot);
%                             axes(haxis(sub2ind([col,row],dumm,sub)));
%                             imshow(uint8(squeeze(frames)));
%                             axis off;
                    end
                else
                    if strcmp(attr,'empty')
%                     if strcmp(attr,'amplitudes')
% 
% %                             [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map,sorted(end-1));
%                         att = eval([attr,'_']);
%                         att_a = eval([attr,'_a_']);
% %                             data = att(1:50,1:128,:,:,:);
%                         if entire_period
%                             on_ = 1;
%                             off_ = 128;
%                         else
%                             if post
%                                 on_ = (16+32)+1;
%                                 off_ = 120;
%                             else
%                                 on_ = 16+1;
%                                 off_ = (16+32);
%                                 on_post_ = (16+32)+1;
%                                 off_post_ = 120;
%                             end
%                         end
%                         if size(att,2)<36
%                             if temp
%                                 data = att(:,on:off,:,:,:);
%                                 data_a = att_a(:,on:off,:,:,:);
%                             else
%                                 data = att(:,on:off,:,:,:);
%                                 data_a = att_a(:,on:off,:,:,:);
%                             end
%                         else
%                             if temp
%                                 data = att(:,on_:off_,:,:,:);
%                                 data_a = att_a(:,on_:off_,:,:,:);
%                             else
%                                 data = att(:,on_:off_,:,:,:);
%                                 data_a = att_a(:,on_:off_,:,:,:);
%                             end
%                         end
%                         Att1 = mean(mean(squeeze(data),1),2);
%                         Att1_a = mean(mean(squeeze(data_a),1),2);
%                         Att1_positive = Att1; Att1_a_positive = Att1_a;
%                         Att1_positive(Att1_positive<0) = 0;
%                         Att1_a_positive(Att1_a_positive<0) = 0;
%                         Att1_negative = Att1; Att1_a_negative = Att1_a;
%                         Att1_negative(Att1_negative>0) = 0;
%                         Att1_a_negative(Att1_a_negative>0) = 0;
%                         Att1_pos = (Att1_positive-Att1_a_positive)./(abs(Att1_positive)+abs(Att1_a_positive));
%                         Att1_neg = (Att1_negative-Att1_a_negative)./(abs(Att1_negative)+abs(Att1_a_negative));
%                         sorted = sort(Att1_pos(:)); maxv = max(sorted(end-1),abs(sorted(2)));
%                         sorted_neg = sort(Att1_neg(:)); maxv_neg = max(sorted_neg(end-1),abs(sorted_neg(2)));
% %                             ind_active = find(mask==1);
% %                             sorted = sort(Att1(ind_active));maxv = min(sorted(end-1),abs(sorted(2)));
%                         Att1_pos = Att1_pos/maxv;
%                         Att1_pos = (Att1_pos+1)/2;
%                         Att1_neg = Att1_neg/maxv_neg;
%                         Att1_neg = (Att1_neg+1)/2;
%                         if temp
%                             data_temp = reshape(data,size(data,1),dsrate,size(data,2)/dsrate,1,size(data,4),size(data,5));
%                             Att1_temp = mean(squeeze(mean(data_temp,2)),1);
%                             sizes = size(Att1_temp);
% %                                 Att1_temp = reshape(Att1_temp,[sizes(1),sizes(3:end)]);
%                             sorted = sort(Att1_temp(:)); maxv = min(sorted(int32(0.99*length(sorted))),abs(sorted(int32(0.01*length(sorted))+1)));
%                             Att1_temp = (Att1_temp/maxv+1)/2;
%                             Att1_cell{sub} = Att1_temp;
%                         else
%                             Att1_cell{sub} = Att1_pos;
%                             Att1_cell_neg{sub} = Att1_neg;
%                         end
%                         mask_cell{sub} = mask;
%                         coord_cell{sub} = coord;
%                         if plot_on_same_brain
%                             if sub == SUB
%                                 [frames] = VisualAtt(Att1_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                                 [frames2] = VisualAtt(Att1_cell_neg, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
% %                                     axes(haxis(sub2ind([col,row],dumm,1)));
% %                                     imshow(uint8(squeeze(frames)));
% %                                     axis off;
%                             end
%                         else
%                             [frames] = VisualAtt({Att1}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
% %                                 axes(haxis(sub2ind([col,row],dumm,sub)));
% %                                 imshow(uint8(squeeze(frames)));
% %                                 axis off;
%                         end
%                         if sub == SUB
%                             frames = cat(4,frames,frames2);
%                         end
                    else
                        if strcmp(attr,'avg')
                            mask_cell{sub} = mask;
                            coord_cell{sub} = coord;
                            if sub == SUB
                                [frames1] = VisualAtt(avg_data, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map_div,1,true,true,annotation,bar_plot);
                                [frames2] = VisualAtt(avg_data, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map_div,1,true,true,annotation,bar_plot);
                                [frames] = VisualAtt(avg_data, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map_div,1,true,true,annotation,bar_plot);
                            end
                        else

    %                             [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map);
    %                         if wreg
    %                             [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median,regions);
    %                         else
    %                             [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median);
    %                         end
                            if temp
                                Att1_cell{sub} = Att1_temp;
                            else
                                Att1_cell{sub} = Att1;
                                att_cell1{sub} = Att1_causal;
                                att_cell2{sub} = Att1_anticausal;
                                
                            end
                            mask_cell{sub} = mask;
                            coord_cell{sub} = coord;
                            if plot_on_same_brain 
                                if sub == SUB
                                    [frames1,prob1,prob0] = VisualAtt(att_cell1, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map,cmax,true,false,annotation,bar_plot);
                                    [frames2,prob2,prob0] = VisualAtt(att_cell2, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map,cmax,true,false,annotation,bar_plot);
                                    if bar_plot
                                        [frames,prob1,prob0] = VisualAtt(Att1_cell, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map_div,1,true,true,annotation,bar_plot);
%                                         [frames,prob1,prob0] = VisualAtt4barplot({att_cell1,att_cell2}, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map_div,1,true,true,annotation,bar_plot);
                                    else
                                        if norm_contr_elec
                                           [frames,prob1,prob0] = VisualAtt4normcontr(Att1_cell,att_cell1,att_cell2, mask_cell, region_cell,coord_cell, plot_data_all_mni,color_map_div,1,false,false,annotation);
                                        else
                                            switch contrast
                                                case 'sin'
                                                    da = sign(prob1'-prob2').*sqrt((prob1'-prob2').^2)./sqrt((prob1'.^2+prob2'.^2)*2);
                                                    da(prob0<1)=0;
                                                    frames={};Lcrtx = load('ch2_template_lh_pial.mat'); Lcrtx.tri = Lcrtx.faces;Lcrtx.vert = Lcrtx.coords; frames{1}=tripatch(Lcrtx, 'nofigure', da); shading interp;colormap(gca,color_map_div);material dull;lighting gouraud;light;axis off;viewangle = 'l';litebrain(viewangle,.8); set(gca,'CLim',[cax(1) cax(2)]); 
                                                case 'sqrt'
                                                    frames={};Lcrtx = load('ch2_template_lh_pial.mat'); Lcrtx.tri = Lcrtx.faces;Lcrtx.vert = Lcrtx.coords; frames{1}=tripatch(Lcrtx, 'nofigure', sign(prob1'-prob2').*(abs(prob1'-prob2')).^the_power); shading interp;colormap(gca,color_map_div);material dull;lighting gouraud;light;axis off;viewangle = 'l';litebrain(viewangle,.8);;set(gca,'CLim',[cax(1) cax(2)]); 
                                                case 'minu'
                                                    frames={};Lcrtx = load('ch2_template_lh_pial.mat'); Lcrtx.tri = Lcrtx.faces;Lcrtx.vert = Lcrtx.coords; frames{1}=tripatch(Lcrtx, 'nofigure', prob1'-prob2'); shading interp;colormap(gca,color_map_div);material dull;lighting gouraud;light;axis off;viewangle = 'l';litebrain(viewangle,.8);set(gca,'CLim',[cax(1) cax(2)]); 
                                            end
                                        end
                                    end
    %                                     axes(haxis(sub2ind([col,row],dumm,1)));
    %                                     imshow(uint8(squeeze(frames)));
    %                                     axis off;
                                end
                            else
                                [frames] = VisualAtt({Att1}, {mask}, {regions}, {coord}, plot_data_all_mni,color_map,1,false,false,annotation,bar_plot);
    %                                 axes(haxis(sub2ind([col,row],dumm,sub)));
    %                                 imshow(uint8(squeeze(frames)));
    %                                 axis off;
                            end
                        end

                    end
                end
            end
        end
%         if ~plot_on_same_brain
%             if dumm == 2
%                 ylabel(Subj)
%             end  
%         end
%         if plot_on_same_brain
%             if sub == SUB
%                 title(attr)
%             end
%         else
%             if sub ==1
%                 title(attr)
%             end
%         end
    end
%     if strcmp(attr,'freq_formants_hamon')
%         dumm = dumm+3;
%     else
%         dumm = dumm+1;
%     end
end

function [frames,Att1_cell] = attr_vis_occbased(attr,post,visdiff,postabs,max_att,plot_on_same_brain,above_median,temp,wreg,avg_data)
%     maxminp1=[];
%     maxminp2=[10000];
%     for sub =1:SUB
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/NY',subjs{sub},'_elec.mat']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/att_NY',subjs{sub},'']);
%         att = eval(attrs{a});
%         if strcmp(attrs{a},'freq_formants_hamon')
%             maxp1 = gather_att(att(:,on:off,1,:,:),postabs,max_att);
%             maxp2 = gather_att(att(:,on:off,2,:,:),postabs,max_att);
%             maxminp1 = [maxminp1,max(maxp1(:))];
%             maxminp2 = [maxminp2,max(maxp2(:))];
%         else
%             maxp1 = gather_att(att(:,on:off,1,:,:),postabs,max_att);
%             maxminp1 = [maxminp1,max(maxp1(:))];
%         end
%     end
%     maxminp1 = min(maxminp1);
%     maxminp2 = min(maxminp2);
    % choose a color map and set indecies of colors for each electrode
    color_map = hot(256);%afmhot;%hot;
    color_map = color_map(end:-1:1,:);
    % [color_map]=cbrewer('seq', 'Blues', 256);
    [color_map_div]=cbrewer('div', 'RdBu', 256,'PCHIP');
    [colorset1] = cbrewer('qual', 'Set1', 9,'PCHIP');
    color_map_div = color_map_div(end:-1:1,:);
    color_ind = value2colorind(ones(128,1),'hot',[0,1]); 
    color_ind = [color_ind,color_ind,color_ind];
    colors = {color_ind, color_map};
    clrmap = color_map;
    if post
        on = 32+1;
        off = 120-16;
    else
        on = 17;
        off = 16+16;
        on_post = 16+32+1;
        off_post = 120;
    end
    dsrate = 8;
    % ind = [1:10,21:30];
    ind = [1:50];%[1:50];
    dumm = 1;
    Att1_cell = {};
    Att2_cell = {};
    diff_cell = {};
    diff2_cell = {};
    mask_cell = {};
    coord_cell = {};
    subjs = {'717','742','749'};
%     subjs = {'749'};
    SUB = length(subjs);
    root_dir = '/Users/ranwang/Documents/writen_paper/NER2020/Visuallization_matlab';
    use_occ = false;%(strcmp(attr,'amplitudes') || strcmp(attr,'f0_hz') || strcmp(attr,'freq_formants_hamon')) && ~temp && ~visdiff && post;
    use_occentire = false;%strcmp(attr,'freq_formants_hamon') && ~temp && ~visdiff && post;
    for sub=1:SUB
        load(['/Users/ranwang/Documents/writen_paper/NER2020/NY_',subjs{sub},'_elec.mat']);
        load(['/Users/ranwang/Documents/writen_paper/NER2020/att_NY',subjs{sub},'']);
%         if use_occentire
%             load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_cc_value_',subjs{sub},'']); 
%             ind = [1:2];
%         else
%             ind = [1:50];
%             if use_occ
%                 load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_dict_IG_value_',subjs{sub},'']);   
%     %             load(['/Users/ranwang/Documents/writen_paper/NER2020/attr_dict_IG_value_',subjs{sub},'_nonnoise_']);
%             else
%                 load(['/Users/ranwang/Documents/writen_paper/NER2020/attr_dict_IG_value_',subjs{sub},'_nonnoise_']);
%     %             load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_dict_IG_value_',subjs{sub},'']);
%             end
%         end
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_dict_IG_value_',subjs{sub},'']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_dict_IG_value_',subjs{sub},'_prearti']);
        
        Subj = ['NY',subjs{sub}];
        % get paths to all visualization files
        paths = get_path(root_dir,Subj);
        channel_info_all = get_channel_info(paths);
        plot_data_all_mni = get_plot_data(paths, 'mni', 'lh');
        plot_data_all_T1 = get_plot_data(paths, 'subj', 'lh');
    %     coord = T1;
        coord = mni;
        aud = isaud(regions);
        mask_nonaud = mask & ~aud;
        if ~strcmp(attr,'avg')
            load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'/occ_value.mat']);
            eval([attr,'_occ=','eval(attr);']);
            load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'/org_value.mat']);
            eval([attr,'_org=','eval(attr);']);
            load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'/select_all_value.mat']);
            eval([attr,'_all=','eval(attr);']);
            att = eval([attr,'_occ']);
            if strcmp(attr,'spec')
                att = squeeze(att);
            end
            att_org = eval([attr,'_org']);
            att_all = eval([attr,'_all']);
        end
        if strcmp(attr,'freq_formants_hamon')
            if use_occ
                Att1 = gather_att2(data,postabs,max_att);
                Att1_temp = Att1;
            else
                [Att1,Att1_temp] = occ2att(att_all(:,on:off,1),att_org(:,on:off,1),att(:,on:off,1,:,:),true);
            end
            sorted = sort(Att1(:));
            % get all plot frames
            if visdiff && ~post
%                 Att1_post = gather_att(att(ind,on_post:off_post,1,:,:),postabs,max_att);
%                 coef = find_normcoef(Att1_post,Att1,mask,max_anker,ankercount);
%                 diff = Att1_post*coef(1)+coef(2) - Att1;
%                 sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                 diff = diff/maxv;
%                 diff = (diff+1)/2;
%                 [frames] = VisualAtt(diff, mask, coord, plot_data_all_mni,color_map_div,1);
                [Att1_post,Att1_post_temp] = gather_att(att(:,on_post:off_post,1,:,:),postabs,max_att,dsrate);
                if wreg
                    [Att1_post,Att1_post_temp,m25,m75] = robust_rescale(Att1_post,Att1_post_temp,mask,true,regions);
                    [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,true,regions);
                else
                    [Att1_post,Att1_post_temp,m25,m75] = robust_rescale(Att1_post,Att1_post_temp,mask,true);
                    [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,true);
                end
                diff = Att1_post-Att1;
                [diff] = checkaud(diff,mask,regions);
                sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
                diff = diff/maxv;
                diff = (diff+1)/2;
                if temp
                    diff_cell{sub} = ((Att1_post_temp-Att1_temp)/maxv+1)/2;
                else
                    diff_cell{sub} = diff;
                end
                mask_cell{sub} = mask;
                coord_cell{sub} = coord;
                if plot_on_same_brain 
                    if sub == SUB
                        [frames] = VisualAtt(diff_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                         axes(haxis(sub2ind([col,row],dumm,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
                    end
                else
                    [frames] = VisualAtt({diff}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                     axes(haxis(sub2ind([col,row],dumm,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
                end
            else
%                 [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map,sorted(end-2));
                if wreg
                    [Att1_md,Att1_md_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,true,regions);
                    [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median,regions);
                else
                    [Att1_md,Att1_md_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask);
                    [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median);
                end
                if temp
                    Att1_cell{sub} = Att1_temp;
                else
                    Att1_cell{sub} = Att1;
                end
                mask_cell{sub} = mask;
                coord_cell{sub} = coord;
                if plot_on_same_brain 
                    if sub == SUB
                        [frames] = VisualAtt(Att1_cell, mask_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false);
%                         axes(haxis(sub2ind([col,row],dumm,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
                    end
                else
                    [frames] = VisualAtt({Att1}, {mask}, {coord}, plot_data_all_mni,color_map,1,false,false);
%                     axes(haxis(sub2ind([col,row],dumm,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
                end
            end
            
            if use_occ
                Att2 = gather_att2(data,postabs,max_att);
                Att2_temp = Att2;
            else
                [Att2,Att2_temp] = occ2att(att_all(:,on:off,2),att_org(:,on:off,2),att(:,on:off,2,:,:),true);
            end
%             [Att2,Att2_temp] = gather_att(data,postabs,max_att,dsrate);
            sorted = sort(Att2(:));
            if visdiff && ~post
%                 Att2_post = gather_att(att(ind,on_post:off_post,2,:,:),postabs,max_att);
%                 coef = find_normcoef(Att2_post,Att2,mask,max_anker,ankercount);
%                 diff = Att2_post*coef(1)+coef(2) - Att2;
%                 sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                 diff = diff/maxv;
%                 diff = (diff+1)/2;
%                 [frames] = VisualAtt(diff, mask, coord, plot_data_all_mni,color_map_div,1);
                [Att2_post,Att2_post_temp] = gather_att(att(:,on_post:off_post,1,:,:),postabs,max_att,dsrate);
                if wreg
                    [Att2_post,Att2_post_temp,m25,m75] = robust_rescale(Att2_post,Att2_post_temp,mask,true,regions);
                    [Att2,Att2_temp,m25,m75] = robust_rescale(Att2,Att2_temp,mask,true,regions);
                else
                    [Att2_post,Att2_post_temp,m25,m75] = robust_rescale(Att2_post,Att2_post_temp,mask);
                    [Att2,Att2_temp,m25,m75] = robust_rescale(Att2,Att2_temp,mask);
                end
                diff = Att2_post-Att2;
                [diff] = checkaud(diff,mask,regions);
                sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
                diff = diff/maxv;
                diff = (diff+1)/2;
                if temp
                    diff_cell{sub} = ((Att2_post_temp-Att2_temp)/maxv+1)/2;
                else
                    diff_cell{sub} = diff;
                end
                mask_cell{sub} = mask;
                coord_cell{sub} = coord;
                if plot_on_same_brain 
                    if sub == SUB
                        [frames2] = VisualAtt(diff_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                         axes(haxis(sub2ind([col,row],dumm+1,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
                    end
                else
                    [frames2] = VisualAtt({diff}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                     axes(haxis(sub2ind([col,row],dumm+1,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
                end
            else
%                 [frames] = VisualAtt(Att2, mask, coord, plot_data_all_mni,color_map,sorted(end-2));
                if wreg
                    [Att2_md,Att2_md_temp,m25,m75] = robust_rescale(Att2,Att2_temp,mask,true,regions);
                    [Att2,Att2_temp,m25,m75] = robust_rescale(Att2,Att2_temp,mask,above_median,regions);
                else
                    [Att2_md,Att2_md_temp,m25,m75] = robust_rescale(Att2,Att2_temp,mask);
                    [Att2,Att2_temp,m25,m75] = robust_rescale(Att2,Att2_temp,mask,above_median);
                end
                if temp
                    Att2_cell{sub} = Att2_temp;
                else
                    Att2_cell{sub} = Att2;
                end
                mask_cell{sub} = mask;
                coord_cell{sub} = coord;
                if plot_on_same_brain 
                    if sub == SUB
                        [frames2] = VisualAtt(Att2_cell, mask_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false);
%                         axes(haxis(sub2ind([col,row],dumm+1,1)));
%                         imshow(uint8(squeeze(frames)));
%                         axis off;
                    end
                else
                    [frames2] = VisualAtt({Att2}, {mask}, {coord}, plot_data_all_mni,color_map,1,false,false);
%                     axes(haxis(sub2ind([col,row],dumm+1,sub)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
                end
            end

            
%             coef = find_normcoef(Att2,Att1,mask,false,20);
%             diff = Att2*coef(1)+coef(2) - Att1;
            
            diff = Att2-Att1;
            sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
            diff = diff/maxv;
            diff = (diff+1)/2;
            if temp
                diff2_cell{sub} = ((Att2_md_temp-Att1_md_temp)/maxv+1)/2;
            else
                diff2_cell{sub} = diff;
            end
            mask_cell{sub} = mask;
            coord_cell{sub} = coord;
            if plot_on_same_brain 
                if sub == SUB
                    [frames3] = VisualAtt(diff2_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                     axes(haxis(sub2ind([col,row],dumm+2,1)));
%                     imshow(uint8(squeeze(frames)));
%                     axis off;
                end
            else
                [frames3] = VisualAtt({diff}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                 axes(haxis(sub2ind([col,row],dumm+2,sub)));
%                 imshow(uint8(squeeze(frames)));
%                 axis off;
            end
            if sub == SUB
                frames = cat(4,frames,frames2,frames3);
            end
            
        else
            if strcmp(attr,'freq_harm_diff')
              disp(' ')
            else
%                 att = eval(attr);
                if ~strcmp(attr,'avg')
                    data = att(:,on:off,1,:,:);
                    sorted = sort(data(:));
                    sorted = sorted(end-1);
                    if use_occ && (strcmp(attr,'f0_hz') || strcmp(attr,'amplitudes'))
                        Att1 = gather_att2(data,postabs,max_att);
                        Att1_temp = Att1;
                    else
                        [Att1,Att1_temp] = occ2att(att_all(:,on:off,:),att_org(:,on:off,:),att(:,on:off,:,:,:),true);
                    end
                    sorted = sort(Att1(:));
    %                 sorted = sorted(end-1);
                end
                % get all plot frames
                if strcmp(attr,'loudness') || strcmp(attr,'f0_hz')
                    if visdiff && ~post
%                         Att1_post = gather_att(att(:,on_post:off_post,1,:,:),postabs,max_att);
%                         coef = find_normcoef(Att1_post,Att1,mask,max_anker,ankercount);
%                         diff = Att1_post*coef(1)+coef(2) - Att1;
%                         sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                         diff = diff/maxv;
%                         diff = (diff+1)/2;
%                         [frames] = VisualAtt(diff, mask, coord, plot_data_all_mni,color_map_div,1);
                        
                        [Att1_post,Att1_post_temp] = gather_att(att(:,on_post:off_post,1,:,:),postabs,max_att,dsrate);

                        if wreg
                            [Att1_post,Att1_post_temp,m25,m75] = robust_rescale(Att1_post,Att1_post_temp,mask,true,regions);
                            [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,true,regions);
                        else
                            [Att1_post,Att1_post_temp,m25,m75] = robust_rescale(Att1_post,Att1_post_temp,mask);
                            [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask);
                        end
                        diff = Att1_post-Att1;
                        [diff] = checkaud(diff,mask,regions);
                        sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
                        diff = diff/maxv;
                        diff = (diff+1)/2;
                        if temp
                            diff_cell{sub} = ((Att1_post_temp-Att1_temp)/maxv+1)/2;
                        else
                            diff_cell{sub} = diff;
                        end
                        mask_cell{sub} = mask;
                        coord_cell{sub} = coord;
                        if plot_on_same_brain 
                            if sub == SUB
                                [frames] = VisualAtt(diff_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                                 axes(haxis(sub2ind([col,row],dumm,1)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
                            end
                        else
                            [frames] = VisualAtt({diff}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                             axes(haxis(sub2ind([col,row],dumm,sub)));
%                             imshow(uint8(squeeze(frames)));
%                             axis off;
                        end
                    else
%                         [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map,maxminp1);
                        if wreg
                            [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median,regions);
                        else
                            [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median);
                        end
                        if temp
                            Att1_cell{sub} = Att1_temp;
                        else
                            Att1_cell{sub} = Att1;
                        end
                        mask_cell{sub} = mask;
                        coord_cell{sub} = coord;
                        if plot_on_same_brain 
                            if sub == SUB
                                [frames] = VisualAtt(Att1_cell, mask_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false);
%                                 axes(haxis(sub2ind([col,row],dumm,1)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
                            end
                        else
                            [frames] = VisualAtt({Att1}, {mask}, {coord}, plot_data_all_mni,color_map,1,false,false);
%                             axes(haxis(sub2ind([col,row],dumm,sub)));
%                             imshow(uint8(squeeze(frames)));
%                             axis off;
                        end
                    end
                else
                    if strcmp(attr,'amplitudes')
                        if visdiff && ~post
%                             Att1_post = gather_att(att(:,on_post:off_post,1,:,:),postabs,max_att);
%                             coef = find_normcoef(Att1_post,Att1,mask,max_anker,ankercount);
%                             diff = Att1_post*coef(1)+coef(2) - Att1;
%                             sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
%                             diff = diff/maxv;
%                             diff = (diff+1)/2;
%                             [frames] = VisualAtt(diff, mask, coord, plot_data_all_mni,color_map_div,1);
                            [Att1_post,Att1_post_temp] = gather_att(att(:,on_post:off_post,1,:,:),postabs,max_att,dsrate);
                            if wreg
                                [Att1_post,Att1_post_temp,m25,m75] = robust_rescale(Att1_post,Att1_post_temp,mask,true,regions);
                                [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,true,regions);
                            else
                                [Att1_post,Att1_post_temp,m25,m75] = robust_rescale(Att1_post,Att1_post_temp,mask);
                                [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask);
                            end
                            diff = Att1_post-Att1;
                            [diff] = checkaud(diff,mask,regions);
                            sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
                            diff = diff/maxv;
                            diff = (diff+1)/2;
                            if temp
                                diff_cell{sub} = ((Att1_post_temp-Att1_temp)/maxv+1)/2;
                            else
                                diff_cell{sub} = diff;
                            end
                            mask_cell{sub} = mask;
                            coord_cell{sub} = coord;
                            if plot_on_same_brain 
                                if sub == SUB
                                    [frames] = VisualAtt(diff_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                                     axes(haxis(sub2ind([col,row],dumm,1)));
%                                     imshow(uint8(squeeze(frames)));
%                                     axis off;
                                end
                            else
                                [frames] = VisualAtt({diff}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                                 axes(haxis(sub2ind([col,row],dumm,sub)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
                            end
                        else
%                             [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map,sorted(end-1));
                            Att1 = mean(mean(squeeze(att_org(:,on:off,:,:,:)-att(:,on:off,:,:,:)),1),2);
                            sorted = sort(Att1(:)); maxv = min(sorted(end-1),abs(sorted(2)));
                            Att1 = Att1/maxv;
                            Att1 = (Att1+1)/2;
                            if temp
                                data_temp = reshape(data,size(data,1),dsrate,size(data,2)/dsrate,1,size(data,4),size(data,5));
                                Att1_temp = mean(mean(squeeze(data_temp),1),2);
                                sizes = size(Att1_temp);
                                Att1_temp = reshape(Att1_temp,[sizes(1),sizes(3:end)]);
                                Att1_temp = (Att1_temp/maxv+1)/2;
                                Att1_cell{sub} = Att1_temp;
                            else
                                Att1_cell{sub} = Att1;
                            end
                            mask_cell{sub} = mask;
                            coord_cell{sub} = coord;
                            if plot_on_same_brain
                                if sub == SUB
                                    [frames] = VisualAtt(Att1_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
%                                     axes(haxis(sub2ind([col,row],dumm,1)));
%                                     imshow(uint8(squeeze(frames)));
%                                     axis off;
                                end
                            else
                                [frames] = VisualAtt({Att1}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
%                                 axes(haxis(sub2ind([col,row],dumm,sub)));
%                                 imshow(uint8(squeeze(frames)));
%                                 axis off;
                            end

                        end
                    else
                        if ~strcmp(attr,'avg')
%                             [frames] = VisualAtt(Att1_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
                            if visdiff && ~post
    %                             Att1_post = gather_att(att(:,on_post:off_post,1,:,:),postabs,max_att);
    %                             coef = find_normcoef(Att1_post,Att1,mask,max_anker,ankercount);
    %                             diff = Att1_post*coef(1)+coef(2) - Att1;
    %                             sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
    %                             diff = diff/maxv;
    %                             diff = (diff+1)/2;
    %                             [frames] = VisualAtt(diff, mask, coord, plot_data_all_mni,color_map_div,1);
                                [Att1_post,Att1_post_temp] = gather_att(att(:,on_post:off_post,1,:,:),postabs,max_att,dsrate);
                                if wreg
                                    [Att1_post,Att1_post_temp,m25,m75] = robust_rescale(Att1_post,Att1_post_temp,mask,true,regions);
                                    [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,true,regions);
                                else
                                    [Att1_post,Att1_post_temp,m25,m75] = robust_rescale(Att1_post,Att1_post_temp,mask);
                                    [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask);
                                end
                                diff = Att1_post-Att1;
                                [diff] = checkaud(diff,mask,regions);
                                sorted = sort(diff(:)); maxv = min(sorted(end-1),abs(sorted(2)));
                                diff = diff/maxv;
                                diff = (diff+1)/2;
                                if temp
                                    diff_cell{sub} = ((Att1_post_temp-Att1_temp)/maxv+1)/2;
                                else
                                    diff_cell{sub} = diff;
                                end

                                mask_cell{sub} = mask;
                                coord_cell{sub} = coord;
                                if plot_on_same_brain
                                    if sub == SUB
                                        [frames] = VisualAtt(diff_cell, mask_cell, coord_cell, plot_data_all_mni,color_map_div,1,true,true);
    %                                     axes(haxis(sub2ind([col,row],dumm,1)));
    %                                     imshow(uint8(squeeze(frames)));
    %                                     axis off;
                                    end
                                else
                                    [frames] = VisualAtt({diff}, {mask}, {coord}, plot_data_all_mni,color_map_div,1,false,true);
    %                                 axes(haxis(sub2ind([col,row],dumm,sub)));
    %                                 imshow(uint8(squeeze(frames)));
    %                                 axis off;
                                end
                            else
    %                             [frames] = VisualAtt(Att1, mask, coord, plot_data_all_mni,color_map);
                                if wreg
                                    [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median,regions);
                                else
                                    [Att1,Att1_temp,m25,m75] = robust_rescale(Att1,Att1_temp,mask,above_median);
                                end
                                if temp
                                    Att1_cell{sub} = Att1_temp;
                                else
                                    Att1_cell{sub} = Att1;
                                end
                                mask_cell{sub} = mask;
                                coord_cell{sub} = coord;
                                if plot_on_same_brain 
                                    if sub == SUB
                                        [frames] = VisualAtt(Att1_cell, mask_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false);
    %                                     axes(haxis(sub2ind([col,row],dumm,1)));
    %                                     imshow(uint8(squeeze(frames)));
    %                                     axis off;
                                    end
                                else
                                    [frames] = VisualAtt({Att1}, {mask}, {coord}, plot_data_all_mni,color_map,1,false,false);
    %                                 axes(haxis(sub2ind([col,row],dumm,sub)));
    %                                 imshow(uint8(squeeze(frames)));
    %                                 axis off;
                                end
                            end
                        else
                            mask_cell{sub} = mask;
                            coord_cell{sub} = coord;
                            if sub == SUB
                                [frames] = VisualAtt(avg_data, mask_cell, coord_cell, plot_data_all_mni,clrmap,1,true,false);
                            end
                        end
                    end
                end
            end
        end
%         if ~plot_on_same_brain
%             if dumm == 2
%                 ylabel(Subj)
%             end  
%         end
%         if plot_on_same_brain
%             if sub == SUB
%                 title(attr)
%             end
%         else
%             if sub ==1
%                 title(attr)
%             end
%         end
    end
%     if strcmp(attr,'freq_formants_hamon')
%         dumm = dumm+3;
%     else
%         dumm = dumm+1;
%     end
end

function [frames,Att1_cell] = region_contrast()
    if ~exist('annotation','var')
        annotation=[];
    end
    color_map = hot(256);%afmhot;%hot;
    color_map = color_map(end:-1:1,:);
    [color_map_div]=cbrewer('div', 'RdBu', 256,'PCHIP');
    [color_map_cb]=cbrewer('seq', 'Greens', 256,'PCHIP');
%     color_map_cb = color_map_cb(end:-1:1,:);
%     color_map = color_map_cb;
%     if annotation
%         color_map = color_map_cb;
%     end
    [colorset1] = cbrewer('qual', 'Set1', 9,'PCHIP');
    color_map_div = color_map_div(end:-1:1,:);
    color_ind = value2colorind(ones(128,1),'hot',[0,1]); 
    color_ind = [color_ind,color_ind,color_ind];
    colors = {color_ind, color_map};
    clrmap = color_map;
    if 0
        dsrate_org = 1;
    else
        dsrate_org = 8;
    end
    if 1
        on = 16/dsrate_org+1;
        off = 120/dsrate_org;
    else
        if post
            on = (16+32)/dsrate_org+1;
            off = 120/dsrate_org;
        else
            on = 16/dsrate_org+1;
            off = (16+32)/dsrate_org;
            on_post = (16+32)/dsrate_org+1;
            off_post = 120/dsrate_org;
        end
    end

    dsrate = 8/dsrate_org;
    Att1_cell = [];
    mask_cell = [];
    coord_cell = [];
    region_cell = [];
    subjs = {'717','742','749'};
%     subjs = {'749'};
    SUB = length(subjs);
    root_dir = '/Users/ranwang/Documents/writen_paper/NER2020/Visuallization_matlab';
    attrs = {'loudness','f0_hz','amplitudes','freq_formants_hamon1','freq_formants_hamon2','freq_formants_noise','bandwidth_formants_noise_hz'};
    colors = linspace(0,1,length(attrs)+1); colors = colors(1:end-1);
    colors = [colors',0.48*ones(length(colors),1),0.5*ones(length(colors),1)];
    colors = hsl2rgb(colors);
    colormaps = [];
    for a=1:length(attrs)
        colormaps = [colormaps;[linspace(1,colors(a,1),128)',linspace(1,colors(a,2),128)',linspace(1,colors(a,3),128)']];
    end
    
    for sub=1:SUB
        
        load(['/Users/ranwang/Documents/writen_paper/NER2020/NY_',subjs{sub},'_elec_entiregrid.mat']);
%         mask = isregion(regions,onregion);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_woaud/','NY',subjs{sub},'_elec_woaud.mat']);

        load(['/Users/ranwang/Documents/writen_paper/NER2020/entiregrid_ecog_',subjs{sub},'percept_active.mat']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_dict_IG_value_',subjs{sub},'']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_dict_IG_value_',subjs{sub},'_prearti']);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/occ_cc_value_prearti_',subjs{sub}]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_/occ_cc_temporal_8.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_/occ_cc_value.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_wpariental/occ_cc_value.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_wpariental/occ_cc_temporal_8.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/entiregrid_',subjs{sub},'_normxdim_anticausal_covnorm_NOVISWEIGHT_step1/occ_cov_temporal_8.mat',]);
        load(['/Users/ranwang/Documents/writen_paper/NER2020/entiregrid_',subjs{sub},'_normxdim_causal_percept_active_covnorm_NOVISWEIGHT_step1/occ_cov_temporal_8.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_normxdim_anticausal/occ_cov_value.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_normxdim/occ_cc_value.mat',]);
%         load(['/Users/ranwang/Documents/writen_paper/NER2020/',subjs{sub},'_prearti/occ_cc_value.mat',]);
        Subj = ['NY',subjs{sub}];
        % get paths to all visualization files
        paths = get_path(root_dir,Subj);
        channel_info_all = get_channel_info(paths);
        plot_data_all_mni = get_plot_data(paths, 'mni', 'lh');
%         plot_data_all_T1 = get_plot_data(paths, 'subj', 'lh');
    %     coord = T1;
        coord = mni;
        region_cell = cat(1,region_cell,regions);
        coord_cell = cat(2,coord_cell,coord);
        mask_cell = cat(1,mask_cell,mask);
        comps = zeros(1,1,15,15,length(attrs));
        for a=1:length(attrs)
            switch attrs{a}
                case 'freq_formants_hamon1'
                    att = eval('freq_formants_hamon');
                    data = att(:,on:off,1,:,:);
                case 'freq_formants_hamon2'
                    att = eval('freq_formants_hamon');
                    data = att(:,on:off,2,:,:);
                otherwise
                    att = eval(attrs{a});
                    data = att(:,on:off,1,:,:);
            end
            [Att1,Att1_temp] = gather_att(data,true,false,dsrate);
            comps(:,:,:,:,a) = Att1;
        end
        Att1_cell = cat(3,Att1_cell,comps);
    end
    Att1_cell = (Att1_cell-mean(Att1_cell,5))./mean(Att1_cell,5);
    [Att1_cell,polar_cell] = max(Att1_cell,[],5);
    Att1_cell = Att1_cell/max(Att1_cell(:));
    Att1_cell = Att1_cell/length(attrs)+(polar_cell-1)/length(attrs);
    [frames] = VisualAtt({Att1_cell}, {mask_cell}, {region_cell}, {coord_cell}, plot_data_all_mni,colormaps,1,true,false,annotation,false,false);
end