function [Graph_Edge] = Build_RAG_Img(label, vertex)
sup_num = max(label(:));
neigh_info = neigh_finding(label);
Graph_Edge = [];
count = 0;
for i = 1:sup_num
    for j = 1:length(neigh_info{i, 1})
        if neigh_info{i, 1}(1, j)~= 0
            count = count+1;
            Graph_Edge(count, 1) = i;
            Graph_Edge(count, 2) = neigh_info{i, 1}(1, j);
            Graph_Edge(count, 3) = norm(vertex(i, :)-vertex(neigh_info{i, 1}(1, j), :));
            idx_tmp = neigh_info{neigh_info{i, 1}(1, j), 1} == i;
            neigh_info{neigh_info{i, 1}(1, j), 1}(1, idx_tmp) = 0;
        end
    end
end