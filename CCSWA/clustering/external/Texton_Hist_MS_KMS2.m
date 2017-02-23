function feature_set = Texton_Hist_MS_KMS2(img_rgb, label)
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
ratio = 2;
for i = 1:ratio:height
    for j = 1:ratio:width
        sub_sample_idx = [sub_sample_idx (i-1)*width+j];
    end
end

threshold = 40;
threshold2 = 8;

bandwidth = 35;
[clustCent,IDX,cluster2dataCell] = MeanShiftCluster(filter_response_col(sub_sample_idx, :)', bandwidth, 0);

K = size(clustCent, 2);
K

while 1
    if(K<threshold2)
        display('Bandwidth too large, recalculating dictionary using mean shift')
        bandwidth = bandwidth*0.9;
        [clustCent,IDX,cluster2dataCell] = MeanShiftCluster(filter_response_col(sub_sample_idx, :)', bandwidth, 0);
        K = size(clustCent, 2);
        K
    elseif(K>threshold)
        display('Bandwidth too small, recalculating dictionary using mean shift')
        bandwidth = bandwidth/0.82;
        [clustCent,IDX,cluster2dataCell] = MeanShiftCluster(filter_response_col(sub_sample_idx, :)', bandwidth, 0);
        K = size(clustCent, 2);
        K
    end
    if(K>=threshold2 && K<=threshold)
        break
    end
end

display('Calculating dictionary using k-means')
[kms_idx, dictionary] = kmeans(filter_response_col, [], 'EmptyAction','singleton', 'start', clustCent');

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