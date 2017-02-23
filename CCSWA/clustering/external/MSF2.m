function Trans_Dist_Mat = MSF2(Graph_Edge, tree_num)
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

idx_SST = cell(tree_num,1);
for i = 1:tree_num
    idx_SST{i,1} = zeros(edge_num,1);
end

for i = 1 : edge_num
    tree_iter = 1;
    while tree_iter<=tree_num
        if connection{tree_iter,1}(Graph_Edge(i, 1), 1) < connection{tree_iter,1}(Graph_Edge(i, 2), 1)
            accept_flag = 0;
            if tree_iter == 1
                accept_flag = 1;
            else
                if rand(1)<=0.8
                    accept_flag = 1;
                end
            end
            if accept_flag == 1
                min_point = Graph_Edge(i, 1);
                max_point = Graph_Edge(i, 2);
                label_min = connection{tree_iter,1}(min_point, 1);
                label_max = connection{tree_iter,1}(max_point, 1);
                Edge_Forest{tree_iter,1}{i, 1} = find(connection{tree_iter,1}==label_min);
                Edge_Forest{tree_iter,1}{i, 2} = find(connection{tree_iter,1}==label_max);
                connection{tree_iter,1}(Edge_Forest{tree_iter,1}{i, 2}, 1) = label_min;
                idx_SST{tree_iter,1}(i,1)=1;
            end
        elseif connection{tree_iter,1}(Graph_Edge(i, 1), 1) > connection{tree_iter,1}(Graph_Edge(i, 2), 1)
            accept_flag = 0;
            if tree_iter == 1
                accept_flag = 1;
            else
                if rand(1)<=0.8
                    accept_flag = 1;
                end
            end
            if accept_flag == 1
                min_point = Graph_Edge(i, 2);
                max_point = Graph_Edge(i, 1);
                label_min = connection{tree_iter,1}(min_point, 1);
                label_max = connection{tree_iter,1}(max_point, 1);
                Edge_Forest{tree_iter,1}{i, 2} = find(connection{tree_iter,1}==label_min);
                Edge_Forest{tree_iter,1}{i, 1} = find(connection{tree_iter,1}==label_max);
                connection{tree_iter,1}(Edge_Forest{tree_iter,1}{i, 1}, 1) = label_min;
                idx_SST{tree_iter,1}(i,1)=1;
            end
        end
        tree_iter = tree_iter+1;
    end
%     continue_flag = 0;
%     for j = 1:tree_num
%         if length(unique(connection{j,1})) ~= 1
%             continue_flag = 1;
%             break
%         end
%     end
%     if continue_flag == 0
%         break
%     end
end

MST_Edge = cell(tree_num,1);
for tree_iter = 1:tree_num
    idx_SST{tree_iter,1} = find(idx_SST{tree_iter,1});
    MST_Edge{tree_iter,1} = Graph_Edge(idx_SST{tree_iter,1}, 1:3);
    Edge_Forest{tree_iter,1} = Edge_Forest{tree_iter,1}(idx_SST{tree_iter,1},:);
end

Trans_Norm = zeros(vertex_num,vertex_num);
for tree_iter = 1:tree_num
    mst_edge_num = size(MST_Edge{tree_iter,1}, 1);
    for i = 1:mst_edge_num
        Trans_Dist_Mat_Set(Edge_Forest{tree_iter,1}{i, 1}, Edge_Forest{tree_iter,1}{i, 2}, tree_iter) = MST_Edge{tree_iter,1}(i, 3);
        Trans_Dist_Mat_Set(Edge_Forest{tree_iter,1}{i, 2}, Edge_Forest{tree_iter,1}{i, 1}, tree_iter) = MST_Edge{tree_iter,1}(i, 3);
        Trans_Norm(Edge_Forest{tree_iter,1}{i, 1}, Edge_Forest{tree_iter,1}{i, 2}) = Trans_Norm(Edge_Forest{tree_iter,1}{i, 1}, Edge_Forest{tree_iter,1}{i, 2})+1;
        Trans_Norm(Edge_Forest{tree_iter,1}{i, 2}, Edge_Forest{tree_iter,1}{i, 1}) = Trans_Norm(Edge_Forest{tree_iter,1}{i, 2}, Edge_Forest{tree_iter,1}{i, 1})+1;
    end
end
% Trans_Dist_Mat = squeeze(max(Trans_Dist_Mat_Set,[],3));
Trans_Dist_Mat = sum(Trans_Dist_Mat_Set,3)./(Trans_Norm+eps);