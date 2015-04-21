clear all, close all, clc

addpath ./data_mat
% LOAD DATA
load Vn_1m_li.mat;
% INITIALIZE
orimap=abs(Vn)>0.02;
feat=[];
nClus = 35; % set # of clusters
m = 500; % # of measurements
Q = 991; % # of frequencies
weight = ones(m, Q); % init weight
penalty2 = 0.5;
% EXTRACT FEATURE
for i=1:size(orimap,1)
    for j=1:size(orimap,2)
        if(orimap(i,j)>0)
            feat=[feat;[i,j]];% DIM (i * j) * 2
        end
    end
end
feat2=[feat(:,1) feat(:,2)-feat(:,1)]; % same DIM with feat
x=feat2'; % DIM 2 * ( i * j )
% CLUSTER
% ADJUST NO. OF CLUSTERS TO 4 FOR 1MHz CASE
[a, model, L] = vbgm_wz_1(x, nClus); % DIM ( i* j ) * 1 
label_name=unique(a); % Get the label name without repetition
kmean_K=length(label_name);% Get the number of labels           
% POLYNOMIAL FIT
% fit for cluster #i

% COUNT NUMBER OF DOTS IN EACH CLUSTER
ClusRed = []; % init reduced # of clusters
for k = 1 : kmean_K
    nDot = length(find(a==k));
    if  nDot > 15 % assume a cluster size threshold
        ClusRed = [ClusRed k]; % record cluster label that larger than threshold size
    end
end
p = zeros(length(ClusRed), 3); % init space for fit coefficients
for k = ClusRed
    xx = feat(find(a==k),2);
    yy = feat(find(a==k),1);
    [p(k,:),~,mu] = polyfit(xx,yy,2);
    valuePred = polyval(p(k,:),xx,[],mu);
%     x1 = linspace(min(xx), max(xx)); % use a finer axis
%     valuePred2 = polyval(p,x1); % predicted values with fitting coefficients
%     errorFit(k,:) = abs(yy - valuePred);
%     erropolyfit(k) = mean((yy-valuePred).^2);

    % PLOT DATA AND POLYNOMIAL FIT
%     subplot(round(sqrt(nClus))+1,round(sqrt(nClus))+1,k) % auto adjust # of subplots, use SQRT for a slower descent than linerly bi-part 
    plot(xx,yy,'o')
    hold on
    plot(xx, valuePred, 'c-')
end

% ADJUST WEIGHT
for k = ClusRed
    xx = feat(find(a==k),2);
    yy = feat(find(a==k),1);
    [p(k,:),~,mu] = polyfit(xx,yy,2);
    for i = 1:Q
        if i>min(xx) && i<max(xx)
             valuePred = polyval(p(k,:),i,[],mu);  
             weight(round(valuePred), i) = penalty2;
        end
    end
end
             
    
    
% ind = []; % init space for indices of penalty
% for k = ClusRed
%     xx = feat(find(a==k),2);
%     yy = feat(find(a==k),1);
%     [p(k,:),~,mu] = polyfit(xx,yy,2);
%     valuePred = polyval(p(k,:),xx,[],mu);
%     weight(f
%     ind = feat(find(a==k),:);
%     weight(ind(:,1), ind(:,2)) = penalty2;
% %     ind = [ind; feat(find(a==k),:)]; % store indices 
% end
% for i = ind(:,1)
%     for j = ind(:,2)
% %         if(i>size(weight,1) && j<size(weight,2))
%             weight(i,j) = penalty2;
% %         end   
%     end
%     
% end

% for i = 1:Q
% %     if(i>min(xx) && i<max(xx))
%         for k = ClusRed 
%             valuePred = polyval(p(k,:),i,[],mu);
%             valuePred2 = round(valuePred);
%             if((valuePred2-1)>=1 && (valuePred2+1)<m)
%                 weight(valuePred2-1:valuePred2+1,i) = penalty2;
%             end
%         end
% %     end       
% end

figure
contour(weight)
figure
X = [1:1:Q];
Y = [1:1:m];
mesh(X,Y,-weight)

% for k=ClusRed
% %    if(errorpolyfit(k)<50)
%        for i=1:Q
%            if(i>min(feat(find(a==k),2)))
%               t1 = P(k,1)*i^2 + P(k,2)*i +P(k,3);
%               t2 = round(t1);
%               %t3=max(t2-1,1);
%               %t4=min(t2+1,size(weight,1));
%               %weight(t3:t4,i)=weight(t3:t4,i)*penalty2; 
%                 if(((t2-1)>=1)&&((t2+1)<size(weight,1))) % # of measurements
%                     weight(t2-1:t2+1,i)=penalty2;
%                 end
%           end
%        end
% %    end
% end
% % ADJUST WEIGHT
% for i=1:Q
%       for k=ClusRed
%           if(i>min(feat(find(a==k),2)))
%                 t1= P(k,1)*i^2 + P(k,2)*i +P(k,3);
%                 t2=round(t1);
%                 weight(t2:t2+1,i)=weight(t2:t2+1,i)*penalty2; 
%                 
%                 % OR
%                 if(((t2-1)>=1)&&((t2+1)<size(weight,1)))
%                 weight(t2-1:t2+1,i)=penalty2;
%                 end
%           end
%           
%       end
% end


% % T = table(xx,yy,f,yy-f,'VariableNames',{'X','Y','Fit','FitError'})
% 
% % hold on
% % plot(xx,f,'r--')
% % plot(x1,f1,'c+')
% % legend('y','f','f1')
% % for k=1:kmean_K
% %     P(k,:) = polyfit(feat(find(a==k),2),feat(find(a==k),1),2);
% %     t1=polyval(P(k,:),feat(find(a==k),2));
% %     errorpolyfit(k)=mean((feat(find(a==k),1)-t1).^2);% Square Mean Error
% % end
% % 
% 
% 
% 