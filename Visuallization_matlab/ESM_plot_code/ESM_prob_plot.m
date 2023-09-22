% plotting
load('MNI.mat')
load('ESM_languagehits_weight.mat');
load('S.mat')
Lcrtx = load('ch2_template_lh_pial.mat');
Rcrtx = load('ch2_template_rh_pial.mat');

chans = MNI;
Rcrtx.tri = Rcrtx.faces;Rcrtx.vert = Rcrtx.coords;
Lcrtx.tri = Lcrtx.faces;Lcrtx.vert = Lcrtx.coords;
% get xyz limits for plotting purposes
perim=1; % how many millimeters away from brain/electrodes boundaries to set the colorcoded plane perimeter, recommend >0 to avoid skimming brain surface (default 1mm)
axl(:,:,1)=[min([Rcrtx.coords; Lcrtx.coords]); max([Rcrtx.coords; Lcrtx.coords])]; %min and max of MESH VERTICES' x y z coordinates (2x3)
axl(:,:,2)=[min(chans); max(chans)]; %%min and max of ELECTRODES' x y z coordinates (2x3) and stack them (2x3x2)
axl=[min(axl(1,:,:),[],3); max(axl(2,:,:),[],3)]; % Get the minima and maxima of both
axislim=reshape(axl,1,6)+[-1 1 -1 1 -1 1]*perim; clear axl %Use the, to define the axis boundaries, and add additional perimeter (perim)
%language
titles = {'speech arrest','visual naming','auditory naming','comprehension'};
S.gsp=18;
[~,prob_all]=ctmr_gauss_plot_edited(Lcrtx,MNI,ones(size(MNI,1),1),S.cax,0,S.cm,S.gsp,[]); 
%%
close all

for i = 1:numel(titles)

set(figure(i),'Position',[1 5 1000 1000]);
weightsCrt = weights(:,i);
% calculate current hit prob. 
S.cax = [0 0.3];
[hh,prob]=ctmr_gauss_plot_edited(Lcrtx,chans,weightsCrt,S.cax,0,S.cm,S.gsp,prob_all); 
alpha(hh,1); % Adjust opacity specified for that row
viewangle = 'l';   
pltshowplanes = 0;
cameratoolbar('setmode',''); 
litebrain(viewangle,.8); 
axis(axislim); 
hold on; colormap(gca,S.cm); 
c = colorbar;
c.Color = 'w';
set(gcf,'Clipping','off')
set(gcf,'InvertHardcopy','off');
set(gcf,'Color','k')
%%
title(titles{i},'Color','w','FontSize',20)
print(gcf,['ESM_lang_LH_3mm' num2str(i) '_prob.jpg'],'-djpeg','-r300');
end