clear all, close all

% paths and data
add_all_path();
% load('Vn2_w4.0-0.9_1M_s0.02.mat');
% load V_bp_2M.mat;
% Vn = V;
load Vn_bp_2M.mat;
% pre-process
orimap=abs(Vn)>0.02;
[fx, fy] = find(orimap > 0);
feat = [fx fy];

%% add Hough 
% % convert each nonzero pixel in orimap to hough space (r,theta)
% houghSpace = 3;
% dTheta = pi/(2*houghSpace);
% %theta = [0:dTheta:pi/2];
% theta = pi/12;
% [M,N] = size(orimap);
% centerX = M/2;
% centerY = N/2;
% r = zeros(length(theta), length(fx));
% for i = 1:length(fx)
% %r(:,i) = ((fx(i)-centerX).*cos(theta)+(fy(i)-centerY).*sin(theta))';   
%     r(:,i) = (fx(i).*cos(theta)+fy(i).*sin(theta))'; 
% end
% % form feature          
% % data_dev2 = r';
% alpha = 5; % weight
% data_dev2 = [feat alpha .* r'];
%============TO Do======================
% use hough and find lines above threshold
% directional & smoothing
%=======================================

%% line projection
projBase = [1 1]; % 45 degree
r = feat*projBase';
% form feature       
alpha = 1.5; % weight
data_dev2 = [feat alpha*r];

% normalize
% for i=1:size(data_dev2,2)
%     data_dev2(:,i)=data_dev2(:,i)./norm(data_dev2(:,i),2);
% end

% might:pre-select
K = 1;
% idx_kms = kmeans(data_dev2,K,'distance','cosine','start','cluster');
% idx_sub = cell(K,1);
% 
% for k = 1:K
%     idx_sub{k,1} = idx_kms==k;
% end

% dist
dist_mat = pdist2(data_dev2, data_dev2);

%% cluster
method = 'TD'
% method = 'VBGMM'
switch lower(method)
    case 'td'
% TD-cluster
    nClus = 50;
    display('Transitive Distance');
    trans_dist_mat = MST_Data(dist_mat);
    opts = statset('Display', 'iter');
    [label, Center, sumDinC, D] = kmeans(trans_dist_mat, nClus, 'distance', 'sqeuclidean', 'emptyaction', 'drop', 'replicates', 3, 'options', opts);
    unique(label)
    length(unique(label))
    case 'vbgmm'
% VB-GMM cluster
nClus = 35;
display('Transitive Distance');
trans_dist_mat = MST_Data(dist_mat);
[a, ~, ~] = vbgm_wz_1(trans_dist_mat, nClus);
unique(a)
label=unique(a);
    otherwise disp('Cluster method not defined!\n')
end

%% visualize
color = 'brgmcyk';
m = length(color);
for i = 1:length(label)
%     plot(feat(a==i,1),feat(a==i,2), 'b.')
    scatter(feat(label==i,1),feat(label==i,2),36,color(mod(i-1,m)+1));
    hold on
end
