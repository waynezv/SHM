function [MST_Edge Trans_Dist_Mat Edge_Forest] = MST_Per(Graph_Edge, Graph_Edge_Orig)
[~, IDX] = sort(Graph_Edge(:, 3)); %col1 & col2 stand for edge vertices; col3 stands for edge weight;
Graph_Edge = Graph_Edge(IDX, :);
Graph_Edge_Orig = Graph_Edge_Orig(IDX, :);
edge_num = size(Graph_Edge, 1);
vertex_num = max(max(Graph_Edge(:,1)), max(Graph_Edge(:,2)));
Trans_Dist_Mat = zeros(vertex_num, vertex_num);

%% Construct MST
connection = (1:vertex_num)';
Edge_Forest = cell(edge_num, 2);
% for i = 1:edge_num
%     Edge_Forest{i,1} = Graph_Edge(i,1);
%     Edge_Forest{i,2} = Graph_Edge(i,2);
% end

for i = 1 : edge_num
    if connection(Graph_Edge(i, 1), 1) < connection(Graph_Edge(i, 2), 1)
        min_point = Graph_Edge(i, 1);
        max_point = Graph_Edge(i, 2);
        label_min = connection(min_point, 1);
        label_max = connection(max_point, 1);
        Graph_Edge(i, 4) = 2;
        Edge_Forest{i, 1} = find(connection==label_min);
        Edge_Forest{i, 2} = find(connection==label_max);
        connection(Edge_Forest{i, 2}, 1) = label_min;
    elseif connection(Graph_Edge(i, 1), 1) > connection(Graph_Edge(i, 2), 1)
        min_point = Graph_Edge(i, 2);
        max_point = Graph_Edge(i, 1);
        label_min = connection(min_point, 1);
        label_max = connection(max_point, 1);
        Graph_Edge(i, 4) = 2;
        Edge_Forest{i, 2} = find(connection==label_min);
        Edge_Forest{i, 1} = find(connection==label_max);
        connection(Edge_Forest{i, 1}, 1) = label_min;
    else
        Graph_Edge(i, 4) = 1;
    end
    if length(unique(connection)) == 1
        break
    end
end

idx_SST = Graph_Edge(:, 4) ==2;
MST_Edge = Graph_Edge_Orig(idx_SST, 1:3);
Edge_Forest = Edge_Forest(idx_SST, 1:2);

mst_edge_num = size(MST_Edge, 1);
for i = 1:mst_edge_num
    for j = 1:length(Edge_Forest{i, 1})
        for k = 1:length(Edge_Forest{i, 2})
            Trans_Dist_Mat(Edge_Forest{i, 1}(j, 1), Edge_Forest{i, 2}(k, 1)) = MST_Edge(i, 3);
            Trans_Dist_Mat(Edge_Forest{i, 2}(k, 1), Edge_Forest{i, 1}(j, 1)) = MST_Edge(i, 3);
        end
    end
end