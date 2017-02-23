function [label_out] = Region_Merging(image, label_in, thre, thre_centroid)
label_out = label_in;
[height width channel] = size(image);
label_num = max(label_in(:));
centroid = zeros(label_num, channel+1);
Adjacency_Matrix = zeros(label_num, label_num);
for i = 1:height
    for j = 1:width
        centroid(label_in(i,j), 1:channel) = centroid(label_in(i,j), 1:channel)+reshape(image(i, j, :), [1 channel]);
        centroid(label_in(i,j), channel+1) = centroid(label_in(i,j), channel+1)+1;
        for k = -1:1
            for l = -1:1
                if i+k>=1 && i+k<=height && j+l>=1 && j+l<=width
                    if label_in(i+k, j+l)~=label_in(i,j)
                        Adjacency_Matrix(label_in(i+k, j+l), label_in(i,j))=1;
                    end
                end
            end
        end
    end
end
centroid(:, 1:channel) = centroid(:, 1:channel)./repmat(centroid(:, channel+1), [1 3]);
neigh_set = cell(label_num, 1);
dist_set = cell(label_num, 1);
for i = 1:label_num
    neigh_set{i, 1}=[];
    dist_set{i, 1}=[];
end
for i = 1:label_num-1
    for j = i:label_num
        if Adjacency_Matrix(i, j)==1
            neigh_set{i, 1} = [neigh_set{i, 1} j];
            dist_set{i, 1} = [dist_set{i, 1} norm(centroid(i, 1:channel)-centroid(j, 1:channel))];
            neigh_set{j, 1} = [neigh_set{j, 1} i];
            dist_set{j, 1} = [dist_set{j, 1} norm(centroid(i, 1:channel)-centroid(j, 1:channel))];
        end
    end
end

iterate = 0;
[min_size min_size_idx] = min(centroid(:, 4));
if min_size<thre
    region_idx = min_size_idx;
    iterate = 1;
else
    search_label=0;
    search_neigh=0;
    tmp_dist=inf;
    for i=1:label_num
        for j=1:length(dist_set{i, 1})
            if dist_set{i, 1}(j)<tmp_dist
                tmp_dist = dist_set{i, 1}(j);
                search_label = i;
                search_neigh = j;
            end
        end
    end
    min_dist = dist_set{search_label, 1}(search_neigh);
    if min_dist<thre_centroid
        iterate = 1;
        region_idx = search_label;
    end
end
while(iterate)
    [edge_weight neigh_idx1] = min(dist_set{region_idx, 1}); %#ok<ASGLU>

    merge_node = neigh_set{region_idx, 1}(neigh_idx1);
    neigh_idx2 = find(neigh_set{merge_node, 1}==region_idx);

    centroid(merge_node, 1:channel) = (centroid(merge_node, 1:channel).*centroid(merge_node, channel+1) + ...
        centroid(region_idx, 1:channel).*centroid(region_idx, channel+1)) ./ (centroid(merge_node, channel+1)+centroid(region_idx, channel+1));
    centroid(merge_node, channel+1) = (centroid(merge_node, channel+1)+centroid(region_idx, channel+1));
    centroid(region_idx, :) = inf;
    label_idx = label_out==region_idx;
    label_out(label_idx)=merge_node;

    neigh_set{region_idx, 1}(neigh_idx1)=[];
    dist_set{region_idx, 1}(neigh_idx1)=[];
    neigh_set{merge_node, 1}(neigh_idx2)=[];
    dist_set{merge_node, 1}(neigh_idx2)=[];
    
    common_node = intersect(neigh_set{region_idx, 1}, neigh_set{merge_node, 1});
    for i = 1:length(common_node)
        idx_del = neigh_set{common_node(i), 1}==region_idx;
        neigh_set{common_node(i), 1}(idx_del) = [];
        dist_set{common_node(i), 1}(idx_del) = [];
        idx_change1 = neigh_set{common_node(i), 1}==merge_node;
        idx_change2 = neigh_set{merge_node, 1}==common_node(i);
        dist_set{common_node(i), 1}(idx_change1) = norm(centroid(common_node(i), 1:channel) - centroid(merge_node, 1:channel));
        dist_set{merge_node, 1}(idx_change2) = dist_set{common_node(i), 1}(idx_change1);
    end

    [diff_node,idx_diff] = setdiff(neigh_set{region_idx, 1}, neigh_set{merge_node, 1});
    for i = 1:length(diff_node)
        d_new = norm(centroid(diff_node(i), 1:channel) - centroid(merge_node, 1:channel));
        neigh_set{merge_node, 1} = [neigh_set{merge_node, 1} diff_node(i)];
        dist_set{merge_node, 1} = [dist_set{merge_node, 1} d_new];
        idx_change = neigh_set{diff_node(i), 1}==region_idx;
        neigh_set{diff_node(i), 1}(idx_change) = merge_node;
        dist_set{diff_node(i), 1}(idx_change) = d_new;
    end
    
    neigh_set{region_idx, 1}=[];
    dist_set{region_idx, 1}=[];

    iterate = 0;
    [min_size min_size_idx] = min(centroid(:, 4));
    if min_size<thre
        iterate = 1;
        region_idx = min_size_idx;
    else
        search_label=0;
        search_neigh=0;
        tmp_dist=inf;
        for i=1:label_num
            for j=1:length(dist_set{i, 1})
                if dist_set{i, 1}(j)<tmp_dist
                    tmp_dist = dist_set{i, 1}(j);
                    search_label = i;
                    search_neigh = j;
                end
            end
        end
        min_dist = dist_set{search_label, 1}(search_neigh);
        if min_dist<thre_centroid
            iterate = 1;
            region_idx = search_label;
        end
    end
end
label_unique = unique(label_out);
for i=1:length(label_unique)
    label_out(label_out==label_unique(i)) = i;
end