% This function eliminates and merges small connected components from the 
% input label, given an input edge map. The dissimilarity between any two
% neighboring regions is defined as mean edge strength (from the edge map)
% over the set of boundary pixels.

function [label_out] = Region_Merging_Edge_Map(image, label_in, bound_pxl_mask, edge_map, thre)

par_k = 240;
label_out = label_in;
[height width channel] = size(image);

[X Y] = meshgrid(1:width, 1:height);

label_num = max(label_in(:));
sup_size = zeros(label_num, 1);
Adjacency_Matrix = zeros(label_num, label_num);
for i = 1:height
    for j = 1:width
        sup_size(label_in(i,j), 1) = sup_size(label_in(i,j), 1)+1;
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
neigh_set = cell(label_num, 1);
dist_set = cell(label_num, 1);
for i = 1:label_num
    neigh_set{i, 1}=[];
    dist_set{i, 1}=[];
end
bandwidth = 4;
for i = 1:label_num-1
    for j = i:label_num
        if Adjacency_Matrix(i, j)==1
            
            label_sup_i = label_out==i;
            label_sup_j = label_out==j;
            size_i = sum(label_sup_i(:));
            size_j = sum(label_sup_j(:));
            label_i = (label_sup_i & bound_pxl_mask==1);
            label_j = (label_sup_j & bound_pxl_mask==1);
            set_xy1 = [X(label_i) Y(label_i)];
            set_xy2 = [X(label_j) Y(label_j)];
            pdist_xy = double(pdist2(set_xy1, set_xy2, 'euclidean')<=sqrt(8));
            normalize_weight = sum(pdist_xy(:));
            edge_strength1 = edge_map(label_i);
            edge_strength2 = edge_map(label_j);
            edge_strength = repmat(edge_strength1,[1, length(edge_strength2)])+repmat(edge_strength2',[length(edge_strength1), 1]);
            dist2 = sum(sum(edge_strength.*pdist_xy))/normalize_weight;
            min_size = min(size_i, size_j);
            dist = dist2-(par_k/min_size);
            
            neigh_set{i, 1} = [neigh_set{i, 1} j];
            dist_set{i, 1} = [dist_set{i, 1} dist];
            neigh_set{j, 1} = [neigh_set{j, 1} i];
            dist_set{j, 1} = [dist_set{j, 1} dist];
        end
    end
end

iter = 0;
%iterate = 1;
while(1)
    iter = iter+1;
    iter
    %iterate = 0;
    
    search_label=0;
    search_neigh=0;
    tmp_dist=inf;
    for i=1:label_num
        if(sup_size(i, 1)==0)
            if(~isempty(dist_set{i, 1}))
                display('There is Error!')
                return;
            end
        end
        for j=1:length(dist_set{i, 1})
            if dist_set{i, 1}(j)<tmp_dist
                tmp_dist = dist_set{i, 1}(j);
                search_label = i;
                search_neigh = j;
            end
        end
    end
    
    
    min_dist = dist_set{search_label, 1}(search_neigh);
    if min_dist<thre
        %iterate = 1;
        region_idx = search_label;
    else
        break;
    end
    
    [~, neigh_idx1] = min(dist_set{region_idx, 1});

    merge_node = neigh_set{region_idx, 1}(neigh_idx1);
    neigh_idx2 = find(neigh_set{merge_node, 1}==region_idx);

    sup_size(merge_node, 1) = sup_size(merge_node, 1) + sup_size(region_idx, 1);
    sup_size(region_idx, 1) = 0;
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

        label_sup_i = label_out==merge_node;
        label_sup_j = label_out==common_node(i);
        size_i = sum(label_sup_i(:));
        size_j = sum(label_sup_j(:));
        label_i = (label_sup_i & bound_pxl_mask==1);
        label_j = (label_sup_j & bound_pxl_mask==1);
        set_xy1 = [X(label_i) Y(label_i)];
        set_xy2 = [X(label_j) Y(label_j)];
        pdist_xy = double(pdist2(set_xy1, set_xy2, 'euclidean')<=sqrt(8));
        normalize_weight = sum(pdist_xy(:));
        edge_strength1 = edge_map(label_i);
        edge_strength2 = edge_map(label_j);
        edge_strength = repmat(edge_strength1,[1, length(edge_strength2)])+repmat(edge_strength2',[length(edge_strength1), 1]);
        dist2 = sum(sum(edge_strength.*pdist_xy))/normalize_weight;
        min_size = min(size_i, size_j);
        dist = dist2-(par_k/min_size);

        dist_set{common_node(i), 1}(idx_change1) = dist;
        dist_set{merge_node, 1}(idx_change2) = dist;
    end

    diff_node1 = setdiff(neigh_set{merge_node, 1}, neigh_set{region_idx, 1});
    diff_node2 = setdiff(neigh_set{region_idx, 1}, neigh_set{merge_node, 1});
    
    
    for i = 1:length(diff_node1)
        idx_change1 = neigh_set{diff_node1(i), 1}==merge_node;
        idx_change2 = neigh_set{merge_node, 1}==diff_node1(i);

        label_sup_i = label_out==merge_node;
        label_sup_j = label_out==diff_node1(i);
        size_i = sum(label_sup_i(:));
        size_j = sum(label_sup_j(:));
        label_i = (label_sup_i & bound_pxl_mask==1);
        label_j = (label_sup_j & bound_pxl_mask==1);
        set_xy1 = [X(label_i) Y(label_i)];
        set_xy2 = [X(label_j) Y(label_j)];
        pdist_xy = double(pdist2(set_xy1, set_xy2, 'euclidean')<=sqrt(8));
        normalize_weight = sum(pdist_xy(:));
        edge_strength1 = edge_map(label_i);
        edge_strength2 = edge_map(label_j);
        edge_strength = repmat(edge_strength1,[1, length(edge_strength2)])+repmat(edge_strength2',[length(edge_strength1), 1]);
        dist2 = sum(sum(edge_strength.*pdist_xy))/normalize_weight;
        min_size = min(size_i, size_j);
        dist = dist2-(par_k/min_size);

        dist_set{diff_node1(i), 1}(idx_change1) = dist;
        dist_set{merge_node, 1}(idx_change2) = dist;
    end
    
    for i = 1:length(diff_node2)
        label_sup_i = label_out==merge_node;
        label_sup_j = label_out==diff_node2(i);
        size_i = sum(label_sup_i(:));
        size_j = sum(label_sup_j(:));
        
        label_i = (label_sup_i & bound_pxl_mask==1);
        label_j = (label_sup_j & bound_pxl_mask==1);
        set_xy1 = [X(label_i) Y(label_i)];
        set_xy2 = [X(label_j) Y(label_j)];
        pdist_xy = double(pdist2(set_xy1, set_xy2, 'euclidean')<=sqrt(8));
        normalize_weight = sum(pdist_xy(:));
        edge_strength1 = edge_map(label_i);
        edge_strength2 = edge_map(label_j);
        edge_strength = repmat(edge_strength1,[1, length(edge_strength2)])+repmat(edge_strength2',[length(edge_strength1), 1]);
        dist2 = sum(sum(edge_strength.*pdist_xy))/normalize_weight;
        min_size = min(size_i, size_j);
        dist = dist2-(par_k/min_size);
        
        neigh_set{merge_node, 1} = [neigh_set{merge_node, 1} diff_node2(i)];
        dist_set{merge_node, 1} = [dist_set{merge_node, 1} dist];
        idx_change = neigh_set{diff_node2(i), 1}==region_idx;
        neigh_set{diff_node2(i), 1}(idx_change) = merge_node;
        dist_set{diff_node2(i), 1}(idx_change) = dist;
    end
    
    neigh_set{region_idx, 1}=[];
    dist_set{region_idx, 1}=[];
end
label_unique = unique(label_out);
for i=1:length(label_unique)
    label_out(label_out==label_unique(i)) = i;
end