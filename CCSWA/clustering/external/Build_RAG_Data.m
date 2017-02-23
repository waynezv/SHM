function [Graph_Edge] = Build_RAG_Data(data)
[point_num, dim] = size(data); %#ok<NASGU>
edge_num = 0;
for i = 1 : point_num-1
    for j = i+1 : point_num
        edge_num = edge_num + 1;
    end
end
Graph_Edge = zeros(edge_num, 3);
count_edge = 0;
for i = 1 : point_num-1
    i %#ok<NOPRT>
    for j = i+1 : point_num
        count_edge = count_edge + 1;
        Graph_Edge(count_edge, 1) = i;
        Graph_Edge(count_edge, 2) = j;
        Graph_Edge(count_edge, 3) = norm(data(i, :) - data(j, :), 2);
    end
end