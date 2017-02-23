function [feature_set] = Filter_Bank_Feature(img_rgb, label)
%% Calculate Filter Response
[filterBank] = createFilterBank2();
[filter_response_col] = extractFilterResponses(img_rgb, filterBank);

% mu = mean(filter_response_col, 1);
% sigma = sqrt(sum((filter_response_col - repmat(mu, [size(filter_response_col, 1), 1])).^2, 1)./(size(filter_response_col, 1)-1));
% filter_response_col = (filter_response_col-repmat(mu, [size(filter_response_col, 1), 1]))./repmat(sigma, [size(filter_response_col, 1), 1]);
sup_num = max(label(:));
feature_set = zeros(sup_num, size(filter_response_col, 2));
for i = 1:sup_num
    idx = label == i;
    feature_set(i, :) = mean(filter_response_col(idx, :), 1);
end