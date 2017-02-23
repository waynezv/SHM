function feature_set = Texton_Hist_MS_EM(img_rgb, label)
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

bandwidth = 28;
[clustCent,IDX,cluster2dataCell] = MeanShiftCluster(filter_response_col(sub_sample_idx, :)', bandwidth, 0);

K = size(clustCent, 2);
K
threshold = 30;
threshold2 = 50;

while 1
    if(K<threshold)
        display('Bandwidth too large, recalculating dictionary using mean shift')
        bandwidth = bandwidth*0.9;
        [clustCent,IDX,cluster2dataCell] = MeanShiftCluster(filter_response_col(sub_sample_idx, :)', bandwidth, 0);
        K = length(unique(IDX));
        K
    elseif(K>threshold2)
        display('Bandwidth too small, recalculating dictionary using mean shift')
        bandwidth = bandwidth/0.86;
        [clustCent,IDX,cluster2dataCell] = MeanShiftCluster(filter_response_col(sub_sample_idx, :)', bandwidth, 0);
        K = length(unique(IDX));
        K
    end
    if(K>=threshold&&K<=threshold2)
        break
    end
end

IDX_Uniq = unique(IDX);
for i = 1:K
    tmp_idx = IDX == IDX_Uniq(i);
    IDX(tmp_idx) = i;
end

display('Calculating dictionary using EM')
obj = gmdistribution.fit(filter_response_col(sub_sample_idx, :),K, 'Regularize', 0.02, 'start', IDX);
[em_idx,nlogl,P] = cluster(obj, filter_response_col);

sup_num = max(label(:));
feature_set = zeros(sup_num, K);
for i = 1:sup_num
    idx = label == i;
    feature_set(i, :) = mean(P(idx, :), 1);
end
toc