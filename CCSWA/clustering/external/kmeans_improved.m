function [IDX] = kmeans_improved(data, k)
cluster_num = 2;
[IDX_tmp C_Orig] = kmeans(data,cluster_num, 'start', 'cluster', 'emptyaction', 'singleton');
while(cluster_num<k)
    split_gain = zeros(cluster_num, 1);
    C_split_set = cell(cluster_num, 1);
    for j = 1:cluster_num
        cluster_idx = IDX_tmp == j;
        if(sum(cluster_idx)>1)
            cluster_data = data(cluster_idx, :);
            MSE_Orig = sum(sum((cluster_data - repmat(C_Orig(j, :), [size(cluster_data, 1) 1])).^2));
            [IDX_split C] = kmeans(cluster_data,2, 'start', 'uniform', 'emptyaction', 'singleton');
            cluster_data1 = cluster_data(IDX_split==1, :);
            cluster_data2 = cluster_data(IDX_split==2, :);
            MSE_Split = sum(sum((cluster_data1 - repmat(C(1, :), [size(cluster_data1, 1) 1])).^2)) + sum(sum((cluster_data2 - repmat(C(2, :), [size(cluster_data2, 1) 1])).^2));
            split_gain(j, 1) = MSE_Orig-MSE_Split;
            C_split_set{j, 1} = C;
        else
            split_gain(j, 1) = -1000000000000;
        end
    end
    [max_val, max_idx] = max(split_gain);
    C_Nxt = C_Orig;
    C_Nxt(max_idx, :) = C_split_set{max_idx, 1}(1, :);
    C_Nxt(cluster_num+1, :) = C_split_set{max_idx, 1}(2, :);
    cluster_num = cluster_num+1;
    [IDX_tmp C_Orig] = kmeans(data,[], 'start', C_Nxt, 'emptyaction', 'singleton');
end
IDX = IDX_tmp;