function [feature_set] = Texton_Hist(img_rgb, label, K)
%% Calculate Filter Response
[filterBank] = createFilterBank2();
[filter_response_col] = extractFilterResponses(img_rgb, filterBank);

% mu = mean(filter_response_col, 1);
% sigma = sqrt(sum((filter_response_col - repmat(mu, [size(filter_response_col, 1), 1])).^2, 1)./(size(filter_response_col, 1)-1));
% filter_response_col = (filter_response_col-repmat(mu, [size(filter_response_col, 1), 1]))./repmat(sigma, [size(filter_response_col, 1), 1]);

%% Calculating Posterior Probability
display('Calculating dictionary using kmeans')
tic
[kms_idx, dictionary] = kmeans(filter_response_col, K, 'EmptyAction','singleton', 'start', 'cluster');
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