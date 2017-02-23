function feature_set = Texton_Hist_MS_KMS(img_rgb, label)
%% Calculate Filter Response
[filterBank] = createFilterBank2();
[filter_response_col] = extractFilterResponses(img_rgb, filterBank);

% mu = mean(filter_response_col, 1);
% sigma = sqrt(sum((filter_response_col - repmat(mu, [size(filter_response_col, 1), 1])).^2, 1)./(size(filter_response_col, 1)-1));
% filter_response_col = (filter_response_col-repmat(mu, [size(filter_response_col, 1), 1]))./repmat(sigma, [size(filter_response_col, 1), 1]);

%% Calculating Posterior Probability
tic
display('Calculating dictionary using mean shift')
height = size(img_rgb, 1);
width = size(img_rgb, 2);
dim = size(filter_response_col, 2);
sub_sample_idx = [];
ratio = 3;
for i = 1:ratio:height
    for j = 1:ratio:width
        sub_sample_idx = [sub_sample_idx (i-1)*width+j];
    end
end
[clustCent,IDX,cluster2dataCell] = MeanShiftCluster(filter_response_col(sub_sample_idx, :)', 35, 0);

K = size(clustCent, 2);
K
threshold = 40;
threshold2 = 4;

if(K>3&&K<=threshold)
    display('Calculating dictionary using k-means')
    [kms_idx, dictionary] = kmeans(filter_response_col, [], 'EmptyAction','singleton', 'start', clustCent');
elseif(K<=threshold2)
    display('Calculating dictionary using k-means')
    [kms_idx, dictionary] = kmeans(filter_response_col, threshold2, 'EmptyAction','singleton', 'start', 'cluster');
    K = threshold2;
elseif(K>threshold)
    display('Calculating dictionary using k-means')
    cluster_num = threshold;
    cluster_size = zeros(K, 1);
    for j=1:K
        idx = IDX == j;
        cluster_size(j, 1) = sum(idx);
    end
    [sorted_size IX] = sort(cluster_size, 'descend');
    Centers = zeros(cluster_num, dim);
    tmp = clustCent';
    for j = 1:cluster_num
        Centers(j, :) = tmp(IX(j), :);
    end
    [kms_idx, dictionary] = kmeans(filter_response_col, [], 'EmptyAction','singleton', 'start', Centers);
    K = threshold;
end

P = pdist2(filter_response_col, dictionary, 'euclidean');
P = (1./P.^2);
P = P./repmat(sum(P, 2), [1 K]);
sup_num = max(label(:));
feature_set = zeros(sup_num, K);
for i = 1:sup_num
    idx = label == i;
    feature_set(i, :) = mean(P(idx, :), 1);
end
toc