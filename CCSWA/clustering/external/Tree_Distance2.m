function [Node_Statistics Node_Connection Dist_Mat] = Tree_Distance2(SST_Edge, hr)
edge_num = size(SST_Edge, 1);
node_num = edge_num + 1;
Node_Statistics = cell(node_num, 1);
Node_Connection = cell(node_num, 2);
for i = 1 : node_num
    Node_Connection{i, 2} = 0;
end
for i = 1 : edge_num
    node1 = SST_Edge(i, 1);
    node2 = SST_Edge(i, 2);
    Node_Connection{node1, 2} = Node_Connection{node1, 2} + 1;
    Node_Connection{node2, 2} = Node_Connection{node2, 2} + 1;
    Node_Connection{node1, 1}(Node_Connection{node1, 2}, 1) = node2;
    Node_Connection{node1, 1}(Node_Connection{node1, 2}, 2) = i;
    Node_Connection{node1, 1}(Node_Connection{node1, 2}, 3) = SST_Edge(i, 3);
    Node_Connection{node2, 1}(Node_Connection{node2, 2}, 1) = node1;
    Node_Connection{node2, 1}(Node_Connection{node2, 2}, 2) = i;
    Node_Connection{node2, 1}(Node_Connection{node2, 2}, 3) = SST_Edge(i, 3);
end
for i = 1 : node_num
    Node_Statistics{i, 1} = cell(Node_Connection{i, 2}, 1);
    for j = 1 : Node_Connection{i, 2}
        if Node_Connection{i, 1}(j, 3) <= hr
            count_node = 1;
            node_curr = Node_Connection{i, 1}(j, 1);
            node_next = [];
            Node_Statistics{i, 1}{j, 1}(1, 1:2) = [node_curr Node_Connection{i, 1}(j, 3)];
            count_next = 0;
            for k = 1 : size(Node_Connection{node_curr, 1}, 1)
                if Node_Connection{node_curr, 1}(k, 1) ~= i && Node_Connection{node_curr, 1}(k, 3) + Node_Connection{i, 1}(j, 3) <= hr
                    count_next = count_next + 1;
                    node_next(count_next, 1:3) = [Node_Connection{node_curr, 1}(k, 1) Node_Connection{node_curr, 1}(k, 3) + Node_Connection{i, 1}(j, 3) node_curr]; %#ok<AGROW>
                end
            end
            while size(node_next, 1) > 0
                Node_Statistics{i, 1}{j, 1}(count_node+1:count_node+size(node_next, 1), :) = node_next(:, 1:2);
                count_node = count_node + size(node_next, 1);
                node_curr = node_next;
                node_next = [];
                count_next = 0;
                for k = 1 : size(node_curr, 1)
                    single_curr_node = node_curr(k, 1);
                    for l = 1 : size(Node_Connection{single_curr_node, 1}, 1)
                        if Node_Connection{single_curr_node, 1}(l, 1) ~= node_curr(k, 3) && Node_Connection{single_curr_node, 1}(l, 3)+node_curr(k, 2) <= hr
                            count_next = count_next + 1;
                            node_next(count_next, 1:3) = [Node_Connection{single_curr_node, 1}(l, 1) Node_Connection{single_curr_node, 1}(l, 3)+node_curr(k, 2) single_curr_node]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end
end

Node_Statistics2 = cell(node_num, 1);
for i = 1 : node_num
    count_store_node = 0;
    temp_store_node = [];
    for j = 1 : size(Node_Statistics{i, 1}, 1)
        for k = 1 : size(Node_Statistics{i, 1}{j, 1}, 1)
            count_store_node = count_store_node + 1;
            temp_store_node(count_store_node, 1:3) = [Node_Statistics{i, 1}{j, 1}(k, :) j];
        end
    end
    Node_Statistics2{i, 1} = temp_store_node;
    Node_Statistics2{i, 1}(count_store_node+1, :) = [i, 0, 0];
end

Dist_Mat = zeros(node_num, node_num);
for i = 1:node_num
    for j = 1:node_num
        Dist_Mat(i, Node_Statistics2{i, 1}(j, 1)) = Node_Statistics2{i, 1}(j, 2);
    end
end