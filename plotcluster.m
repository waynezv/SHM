function plotcluster(feat, feat2, i)
% PLOTFIT TAKES FEAT2 AND PLOTS RESULTS OF CLUSTERING 
% i INDICATES THE ith FIGURE TO BE PLOT

x=feat2';
[a, model, L] = vbgm_wz_1(x, 4); % TUNE NO. OF CLUSTERS
label_name=unique(a);
noCluster=length(label_name); % NO. OF CLUSTERS
spread_wz_1(flipud(x), a, i); % PLOT SCATTERS OF CLUSTERING RESULT
% spread(flipud(x), a)

end