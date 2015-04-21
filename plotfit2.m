function plotfit2(feat, feat2, i)
% PLOTFIT TAKES FEAT2 AND PLOTS RESULTS OF CLUSTERING AND POLYFITTING
% i INDICATES THE ith FIGURE TO PLOT

x=feat2';
[a, model, L] = vbgm_wz_1(x, 4); % TUNE NO. OF CLUSTERS
label_name=unique(a);
noCluster=length(label_name); % NO. OF CLUSTERS
% subplot(3,3,i);
for k=1:noCluster
    [P(k,:)]= polyfit(feat(find(a==k),2),feat(find(a==k),1),2);
    t1=polyval(P(k,:),feat(find(a==k),2));
    errorpolyfit(k)=mean((feat(find(a==k),1)-t1).^2);% Square Mean Error
    plot(feat(find(a==k),2), t1);
    hold on;
end
end