function dictionary = learn_dict_ms_kms_hard(img_rgb)
%% Calculate Filter Response
[filterBank] = createFilterBank2();
[filter_response_col] = extractFilterResponses(img_rgb, filterBank);

%% Calculating Posterior Probability
display('Calculating dictionary using mean shift')
height = size(img_rgb, 1);
width = size(img_rgb, 2);
sub_sample_idx = [];
ratio = 2;
for i = 1:ratio:height
    for j = 1:ratio:width
        sub_sample_idx = [sub_sample_idx (i-1)*width+j];
    end
end

threshold = 40;
threshold2 = 20;

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
        bandwidth = bandwidth/0.85;
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