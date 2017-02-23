function [IDX] = kmeans4(D, k1, k2)
IDX_Set = cell(k2,1);
MSE = zeros(k2,1);
for i = 1:k2
    [IDX_Set{i,1}, ~, sumd] = kmeans(D, k1, 'emptyaction', 'singleton', 'distance', 'cosine');
    MSE(i,1) = sum(sumd);
end
[~, min_idx] = min(MSE);
IDX = IDX_Set{min_idx,1};