load('active.mat');
load('passive.mat');
load('imagine.mat');
load('coordinate.mat');
load('mask.mat');
load('region.mat');
invest_regions = {'cSTG','rSTG','mMTG','cMTG','rMTG',...
                  'inferiorprecentral','superiorprecentral','postcentral','supramarginal',...
                  'parsopercularis','parstriangularis','parsorbitalis','rostralmiddlefrontal','caudalmiddlefrontal'};
act=[];
pass=[];
imgn=[];
msk=[];
cod=[];
reg = [];
for sub = 1:length(active)
   act = cat(3,act,active{sub});
   pass = cat(3,pass,passive{sub});
   imgn = cat(3,imgn,imagine{sub});
   msk = cat(1,msk,mask{sub});
   cod = cat(2,cod,coordinate{sub});
   reg = cat(1,reg,region{sub});
end
act=squeeze(act);
pass=squeeze(pass);
imgn=squeeze(imgn);
region_value_active = cell(length(invest_regions),1);
region_value_passive = cell(length(invest_regions),1);
region_value_imagine = cell(length(invest_regions),1);
for m=1:size(msk,1)
   for n=1:size(msk,2)
      region = reshape(reg(m,n,:),1,[]);
      [found,ind]=ismember(region(find(region~=0)),invest_regions);
      [ismstg,~]=ismember(region(find(region~=0)),{'mSTG'});
      if ismstg
         if cod(2,m,n)>=0.553486529*cod(3,m,n)-2.527049117
            [~,ind]=ismember('rSTG',invest_regions);
         else
            [~,ind]=ismember('cSTG',invest_regions);
         end
      end
      if found && msk(m,n)
         region_value_active{ind} = [region_value_active{ind},act(m,n)];
         region_value_imagine{ind} = [region_value_imagine{ind},imgn(m,n)];
         region_value_passive{ind} = [region_value_passive{ind},pass(m,n)];
      end
   end
end
for r=1:length(region_value_active)
   p=zeros(1000,1);z=zeros(1000,1);
   n=50;
   for iii=1:1000
%       data=[randn(size(region_value_active{r}',2)*n,1)*std(region_value_active{r})'+mean(region_value_active{r})', ...
%             randn(size(region_value_imagine{r}',2)*n,1)*std(region_value_imagine{r})'+mean(region_value_imagine{r})', ...
%             randn(size(region_value_passive{r}',2)*n,1)*std(region_value_passive{r})'+mean(region_value_passive{r})'];
        data=[randn(size(region_value_imagine{r}',2)*n,1)*std(region_value_imagine{r})'+mean(region_value_imagine{r})', ...
            randn(size(region_value_passive{r}',2)*n,1)*std(region_value_passive{r})'+mean(region_value_passive{r})'];
      [p(iii),~,stats] = anova1(data,[1,2],'off');
   end
   disp(mean(p));
   
%    multcompare(stats);
end