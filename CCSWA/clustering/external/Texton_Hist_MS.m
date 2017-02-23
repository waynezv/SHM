function feature_set = Texton_Hist_MS(img_rgb, label)
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
ratio = 1;
for i = 1:ratio:height
    for j = 1:ratio:width
        sub_sample_idx = [sub_sample_idx (i-1)*width+j];
    end
end

bandwidth = 35;
[clustCent,IDX,cluster2dataCell] = MeanShiftCluster(filter_response_col(sub_sample_idx, :)', bandwidth, 0);

K = size(clustCent, 2);
K
threshold = 40;
threshold2 = 5;

while 1
    if(K<threshold2)
        display('Bandwidth too large, recalculating dictionary using mean shift')
        bandwidth = bandwidth*0.9;
        [clustCent,IDX,cluster2dataCell] = MeanShiftCluster(filter_response_col(sub_sample_idx, :)', bandwidth, 0);
        K = size(clustCent, 2);
        K
        if(K>=threshold2)
            break
        end
    elseif(K>threshold)
        display('Bandwidth too small, recalculating dictionary using mean shift')
        bandwidth = bandwidth/0.9;
        [clustCent,IDX,cluster2dataCell] = MeanShiftCluster(filter_response_col(sub_sample_idx, :)', bandwidth, 0);
        K = size(clustCent, 2);
        K
        if(K<=threshold2)
            break
        end
    else
        break
    end
end

sample_num = length(IDX);
IDX2 = zeros(1, sample_num);
unique_IDX = unique(IDX);
dict_num = length(unique_IDX);
iter = 0;
P = zeros(sample_num, dict_num);
for i=1:dict_num
    iter = iter+1;
    tmp_idx = IDX == unique_IDX(i);
    IDX2(tmp_idx) = iter;
    P(tmp_idx, iter) = 1;
end

sup_num = max(label(:));
feature_set = zeros(sup_num, K);
for i = 1:sup_num
    idx = label == i;
    feature_set(i, :) = mean(P(idx, :), 1);
end
toc