clear all, close all, clc;
load feat_1m_li, load feat2_1m_li;

% orimap=abs(Vn)>0.02;
% 
% feat=[];
% for i=1:size(orimap,1)
%     for j=1:size(orimap,2)
%         if(orimap(i,j)>0)
%             feat=[feat;[i,j]];
%         end
%     end
% end
% 
% feat2=[feat(:,1) feat(:,2)-feat(:,1)];

% SHOT 1:9
% figure(1)
% for i = 1:9
% plotcluster(feat, feat2, i);
% end
% 
% figure(2)
% for i = 1:9
% plotfit(feat, feat2, i);
% end
% 
% figure(3)
for i = 1:9
plotweight(feat, feat2, i);
end
