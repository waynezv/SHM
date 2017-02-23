function [IDX] = agglomerative_clustering(D, K, thre)
display('Performing Agglomerative Clustering')
sup_num = size(D, 1);
IDX_tmp = zeros(sup_num, 1);
for i=1:sup_num
    IDX_tmp(i, 1) = i;
end
cluster_num = sup_num;
label_set = unique(IDX_tmp);
outlier_set = [];
cost_set = zeros(sup_num, 1);

count_merge = 0;
while(cluster_num>K)
    count_merge = count_merge+1;
    min_cost_diff = inf;
    min_cost = inf;
    min_label1 = 0;
    min_member_set2 = [];
    min_outlier = [];
    for i=1:cluster_num-1
        label1 = label_set(i);
        member_set1 = find(IDX_tmp==label1);
        for j=i+1:cluster_num
            label2 = label_set(j);
            member_set2 = find(IDX_tmp==label2);
            member_set = [member_set1; member_set2];
            dist_tmp = sum(D(member_set, member_set), 2)./length(member_set);
            idx_outlier = dist_tmp > thre;
            idx_inlier = 1-idx_outlier;
            tmp_outlier = member_set(idx_outlier);
            curr_cost = sum(dist_tmp(logical(idx_inlier))) + length([outlier_set; tmp_outlier])*thre;
            cost_diff = curr_cost - cost_set(label1) - cost_set(label2);
            if(cost_diff<min_cost_diff)
                min_label1 = label1;
                min_member_set2 = member_set2;
                min_cost_diff = cost_diff;
                min_cost = curr_cost;
                min_outlier = tmp_outlier;
            end
        end
    end
    IDX_tmp(min_member_set2) = min_label1;
    label_set = unique(IDX_tmp(IDX_tmp>0));
    cost_set(label1) = min_cost;
    cluster_num = length(label_set);
    cluster_num
    IDX_tmp(min_outlier) = 0;
    outlier_set = [outlier_set; min_outlier];
end
IDX = IDX_tmp;
for i=1:cluster_num
    label = label_set(i);
    idx = IDX_tmp==label;
    IDX(idx) = i;
end
for i=1:sup_num
    if(IDX(i)==0)
        min_dist = inf;
        merge_label = 0;
        for j=1:sup_num
            if(IDX(j)~=0)
                if(D(i, j)<min_dist)
                    min_dist = D(i, j);
                    merge_label = IDX(j);
                end
            end
        end
        IDX(i)=merge_label;
    end
end