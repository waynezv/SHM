function [Graph_Edge] = Build_RAG_Img_Hybrid(label, vertex, image, bound_pxl_mask)
[height width channel] = size(image);
R = image(:,:,1);
G = image(:,:,2);
B = image(:,:,3);
[L,a,b] = RGB2Lab(double(R), double(G), double(B));
[X Y] = meshgrid(1:width, 1:height);

sup_num = max(label(:));
neigh_info = neigh_finding(label);
bandwidth = 4;
Graph_Edge = [];
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
            set_rgb1 = [L(label_i) a(label_i) b(label_i)];
            set_xy1 = [X(label_i) Y(label_i)];
            set_rgb2 = [L(label_j) a(label_j) b(label_j)];
            set_xy2 = [X(label_j) Y(label_j)];
            pdist_rgb = pdist2(set_rgb1, set_rgb2, 'euclidean');
            pdist_xy = pdist2(set_xy1, set_xy2, 'euclidean');
            pdist_xy = exp(-(pdist_xy.^2)./(2*bandwidth^2));
            normalize_weight = sum(pdist_xy(:));
            dist_edge = sum(sum(pdist_rgb.*pdist_xy))/normalize_weight;
            % Compute Chi2 Distance
            dist_chi = pdist3(vertex(node1, :), vertex(node2, :), 'chisq');            
            % Update Edge and Neighbor Info
            Graph_Edge(count, 1) = node1;
            Graph_Edge(count, 2) = node2;
            Graph_Edge(count, 3) = dist_edge;
            idx_tmp = neigh_info{neigh_info{i, 1}(1, j), 1} == i;
            neigh_info{neigh_info{i, 1}(1, j), 1}(1, idx_tmp) = 0;
        end
    end
end