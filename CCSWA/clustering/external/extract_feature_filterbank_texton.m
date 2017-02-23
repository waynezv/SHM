function [feature_set coordinate_set] = extract_feature_filterbank_texton(img, label)
[width height channel] = size(img);
sup_num = max(max(label));
[X Y] = meshgrid(1:width, 1:height);
filterBank = createFilterBank();
filter_response = extractFilterResponses(img, filterBank);
feature_set = zeros(sup_num, size(filter_response, 2));
coordinate_set = zeros(sup_num, 2);
for i = 1:sup_num
    label_idx = label == i;
    feature_set(i, :) = mean(filter_response(label_idx(:), :), 1);
    coordinate_set(i, 1) = mean(X(label_idx(:)));
    coordinate_set(i, 2) = mean(Y(label_idx(:)));
end