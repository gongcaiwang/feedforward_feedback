load('curve_f0.mat');
f0 = frames;
load('curve_fdiff.mat');
fdiff = frames;
load('curve_f0_err.mat');
f0_err = frames;
load('curve_fdiff_err.mat');
fdiff_err = frames;
name4plot = {'cSTG','rSTG','mMTG','cMTG','rMTG','vPreCG','dPreCG','PostCG','supramarginal','pOp','pTri','pOb','rMFG','cMFG'};
for i=1:14
    figure;
    plot(squeeze(f0(:,i))); hold on;
    plot(squeeze(fdiff(:,i)));
    title(name4plot{i})
end