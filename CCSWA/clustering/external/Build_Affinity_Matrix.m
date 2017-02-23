function [Affinity_Matrix] = Build_Affinity_Matrix(label, vertex, bound_pxl_mask, edge_map, bandwidth)
[height width] = size(label);
[X Y] = meshgrid(1:width, 1:height);

sup_num = max(label(:));
neigh_info = neigh_finding(label);
Affinity_Matrix = zeros(sup_num,sup_num);
count = 0;
for i = 1:sup_num
    for j = 1:length(neigh_info{i, 1})
        if neigh_info{i, 1}(1, j)~= 0
            count = count+1;
            node1 = i;
            node2 = neigh_info{i, 1}(1, j);
            % Compute Edge Strength Distance
            label_i = (label==node1 & bound_pxl_mask==1);
            label_j = (label==node2 & bound_pxl_mask==1);
            set_xy1 = [X(label_i) Y(label_i)];
            set_xy2 = [X(label_j) Y(label_j)];
            edge_strength1 = edge_map(label_i);
            edge_strength2 = edge_map(label_j);
            edge_strength = repmat(edge_strength1,[1, length(edge_strength2)])+repmat(edge_strength2',[length(edge_strength1), 1]);
            pdist_xy = double(pdist2(set_xy1, set_xy2, 'euclidean')<=sqrt(8));
            normalize_weight = sum(pdist_xy(:));
            dist_edge = sum(sum(edge_strength.*pdist_xy))/normalize_weight;
            % Compute Chi2 Distance
            dist_chi = pdist3(vertex(node1, :), vertex(node2, :), 'chisq');
            % Update Edge and Neighbor Info
            Affinity_Matrix(node1, node2) = exp(-(dist_edge+0.5*dist_chi)/(2*bandwidth^2));
            Affinity_Matrix(node2, node1) = Affinity_Matrix(node1, node2);
            idx_tmp = neigh_info{neigh_info{i, 1}(1, j), 1} == i;
            neigh_info{neigh_info{i, 1}(1, j), 1}(1, idx_tmp) = 0;
        end
    end
end
Affinity_Matrix = Affinity_Matrix+diag(ones(sup_num,1));