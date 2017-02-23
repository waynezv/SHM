function feature_set = Texton_Hist_Soft_Overlap(img_rgb, dictionary, label, win_size)
[height width channel] = size(img_rgb);
[X Y] = meshgrid(1:width, 1:height);
[filterBank] = createFilterBank2();
[filter_response_col] = extractFilterResponses(img_rgb, filterBank);
K = size(dictionary, 1);

P = pdist2(filter_response_col, dictionary, 'euclidean');
P = (1./P.^2);
P = P./repmat(sum(P, 2), [1 K]);
sup_num = max(label(:));
feature_set = zeros(sup_num, K);
for i = 1:sup_num
    idx = label == i;
    mean_x = round(mean(X(idx)));
    mean_y = round(mean(Y(idx)));
    min_x = max(mean_x-win_size, 1);
    max_x = min(mean_x+win_size, width);
    min_y = max(mean_y-win_size, 1);
    max_y = min(mean_y+win_size, height);
    tmp_idx = logical(zeros(height, width));
    tmp_idx(min_y:max_y,min_x:max_x)=1;
    feature_set(i, :) = mean(P(tmp_idx(:), :), 1);
end