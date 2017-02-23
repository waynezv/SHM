function plotweight(feat, feat2, i)
% PLOTFIT TAKES FEAT2 AND PLOTS RESULTS OF WEIGHTING AFTER CLUSTERING 
% AND POLYFITTING. 
% i INDICATES THE ith FIGURE TO PLOT

x=feat2';
[a, model, L] = vbgm_wz_1(x, 4); % TUNE NO. OF CLUSTERS
label_name=unique(a);
noCluster=length(label_name); % NO. OF CLUSTERS

figure(1)
spread_wz_1(flipud(x), a, i); % PLOT SCATTERS OF CLUSTERING RESULT

figure(2)
subplot(3,3,i);
for k=1:noCluster
    [P(k,:)]= polyfit(feat(find(a==k),2),feat(find(a==k),1),2);
    t1=polyval(P(k,:),feat(find(a==k),2));
    errorpolyfit(k)=mean((feat(find(a==k),1)-t1).^2);% Square Mean Error
    plot(feat(find(a==k),2), t1);
    hold on;
end

figure(3)
subplot(3,3,i);
weight = ones(500,991);
penalty2 = 0.5;
for k=1: noCluster
       if(errorpolyfit(k)<50)
           for i=1:991 % Q
               if(i>min(feat(find(a==k),2)))
                  t1 = P(k,1)*i^2 + P(k,2)*i +P(k,3);
                  t2 = round(t1);
                  %t3=max(t2-1,1);
                  %t4=min(t2+1,size(weight,1));
                  %weight(t3:t4,i)=weight(t3:t4,i)*penalty2; 
                    if(((t2-1)>=1)&&((t2+1)<size(weight,1)))
                        weight(t2-1:t2+1,i)=penalty2;
                    end
              end
           end
       end
end
contour(weight);

end