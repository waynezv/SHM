function feature_set = Texton_Hist_Soft(img_rgb, dictionary, label)
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
    feature_set(i, :) = mean(P(idx, :), 1);
end