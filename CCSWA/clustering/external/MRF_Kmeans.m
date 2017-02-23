function [IDX] = MRF_Kmeans(data, neigh, K, alpha)
[IDX Center] = kmeans(data, K, 'start', 'cluster', 'emptyaction', 'singleton');
Center_Pre = Center;
IDX_Pre = IDX;
sup_num = size(data, 1);
epsilon = 0.00001;
error = 10000;
iter = 0;
while(iter<20)
    iter = iter+1;
    iter
    for i = 1:sup_num
        dist_K = zeros(1, K);
        dist_K_Orig = zeros(1, K);
        neigh_size = size(neigh{i, 1}, 2);
        nei_label_check = zeros(1, neigh_size);
        added_penalty = zeros(1, K);
        for k = 1:K
            dist_K(1, k) = norm(data(i, :) - Center(k,:));
            dist_K_Orig(1, k) = dist_K(1, k);
            for l = 1:neigh_size
                if(k~=IDX_Pre(neigh{i, 1}(1, l)))
                    dist_K(1, k) = dist_K(1, k) + alpha*neigh{i, 2}(1, l);
                    added_penalty(1, k) = added_penalty(1, k)+1;%alpha*neigh{i, 2}(1, l);
                end
                nei_label_check(1, l) = IDX_Pre(neigh{i, 1}(1, l));
            end
        end
        [dist_val label] = min(dist_K);
        IDX(i) = label;
    end
    check = find(IDX_Pre ~= IDX);
    for k = 1:K
        idx = IDX == k;
        Center(k, :) = mean(data(idx, :), 1);
    end
    
    
    IDX_Pre = IDX;

    error = sqrt(sum(sum((Center - Center_Pre).^2)));
    Center_Pre = Center;
end