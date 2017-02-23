function Trans_Dist_Mat = MSF(Graph_Edge, tree_num)
[~, IDX] = sort(Graph_Edge(:, 3)); %col1 & col2 stand for edge vertices; col3 stands for edge weight;
Graph_Edge = Graph_Edge(IDX, :);
edge_num = size(Graph_Edge, 1);
vertex_num = max(max(Graph_Edge(:,1)), max(Graph_Edge(:,2)));
Trans_Dist_Mat_Set = zeros(vertex_num, vertex_num, tree_num);

%% Construct MST
connection = cell(tree_num,1);
Edge_Forest = cell(tree_num,1);
for i = 1:tree_num
    connection{i,1} = (1:vertex_num)';
    Edge_Forest{i,1} = cell(edge_num,2);
end

tree_iter = 1;
for i = 1 : edge_num
    tmp = 0;
    while tmp<tree_num
        if connection{tree_iter,1}(Graph_Edge(i, 1), 1) < connection{tree_iter,1}(Graph_Edge(i, 2), 1)
            min_point = Graph_Edge(i, 1);
            max_point = Graph_Edge(i, 2);
            label_min = connection{tree_iter,1}(min_point, 1);
            label_max = connection{tree_iter,1}(max_point, 1);
            Graph_Edge(i, 4) = tree_iter;
            Edge_Forest{tree_iter,1}{i, 1} = find(connection{tree_iter,1}==label_min);
            Edge_Forest{tree_iter,1}{i, 2} = find(connection{tree_iter,1}==label_max);
            connection{tree_iter,1}(Edge_Forest{tree_iter,1}{i, 2}, 1) = label_min;
            tree_iter = mod(tree_iter,tree_num)+1;
            break
        elseif connection{tree_iter,1}(Graph_Edge(i, 1), 1) > connection{tree_iter,1}(Graph_Edge(i, 2), 1)
            min_point = Graph_Edge(i, 2);
            max_point = Graph_Edge(i, 1);
            label_min = connection{tree_iter,1}(min_point, 1);
            label_max = connection{tree_iter,1}(max_point, 1);
            Graph_Edge(i, 4) = tree_iter;
            Edge_Forest{tree_iter,1}{i, 2} = find(connection{tree_iter,1}==label_min);
            Edge_Forest{tree_iter,1}{i, 1} = find(connection{tree_iter,1}==label_max);
            connection{tree_iter,1}(Edge_Forest{tree_iter,1}{i, 1}, 1) = label_min;
            tree_iter = mod(tree_iter,tree_num)+1;
            break
        else
            tmp = tmp + 1;
            tree_iter = mod(tree_iter,tree_num)+1;
        end
        if tmp >= tree_num
            Graph_Edge(i, 4) = -1;
        end
    end
    continue_flag = 0;
    for j = 1:tree_num
        if length(unique(connection{j,1})) ~= 1
            continue_flag = 1;
            break
        end
    end
    if continue_flag == 0
        break
    end
end

idx_SST = cell(tree_num,1);
MST_Edge = cell(tree_num,1);
for tree_iter = 1:tree_num
    idx_SST{tree_iter,1} = Graph_Edge(:,4) == tree_iter;
    MST_Edge{tree_iter,1} = Graph_Edge(idx_SST{tree_iter,1}, 1:3);
    Edge_Forest{tree_iter,1} = Edge_Forest{tree_iter,1}(idx_SST{tree_iter,1},:);
end

% Trans_Norm = zeros(vertex_num,vertex_num);
for tree_iter = 1:tree_num
    mst_edge_num = size(MST_Edge{tree_iter,1}, 1);
    for i = 1:mst_edge_num
        Trans_Dist_Mat_Set(Edge_Forest{tree_iter,1}{i, 1}, Edge_Forest{tree_iter,1}{i, 2}, tree_iter) = MST_Edge{tree_iter,1}(i, 3);
        Trans_Dist_Mat_Set(Edge_Forest{tree_iter,1}{i, 2}, Edge_Forest{tree_iter,1}{i, 1}, tree_iter) = MST_Edge{tree_iter,1}(i, 3);
%         Trans_Norm(Edge_Forest{tree_iter,1}{i, 1}, Edge_Forest{tree_iter,1}{i, 2}) = Trans_Norm(Edge_Forest{tree_iter,1}{i, 1}, Edge_Forest{tree_iter,1}{i, 2})+1;
%         Trans_Norm(Edge_Forest{tree_iter,1}{i, 2}, Edge_Forest{tree_iter,1}{i, 1}) = Trans_Norm(Edge_Forest{tree_iter,1}{i, 2}, Edge_Forest{tree_iter,1}{i, 1})+1;
%         for j = 1:length(Edge_Forest{tree_iter,1}{i, 1})
%             for k = 1:length(Edge_Forest{tree_iter,1}{i, 2})
%                 Trans_Dist_Mat_Set(tree_iter, Edge_Forest{tree_iter,1}{i, 1}(j, 1), Edge_Forest{tree_iter,1}{i, 2}(k, 1)) = MST_Edge(i, 3);
%                 Trans_Dist_Mat_Set(tree_iter, Edge_Forest{tree_iter,1}{i, 2}(k, 1), Edge_Forest{tree_iter,1}{i, 1}(j, 1)) = MST_Edge(i, 3);
%             end
%         end
    end
end
Trans_Dist_Mat = squeeze(max(Trans_Dist_Mat_Set,[],3));
% Trans_Dist_Mat = sum(Trans_Dist_Mat_Set,3)./(Trans_Norm+eps);