function neigh_info = neigh_finding(sup_label)
[height width] = size(sup_label);
sup_num = max(sup_label(:));
connection_matrix = zeros(sup_num, sup_num);
for i = 1:height
    for j =1:width
        for k=-1:1
            for l=-1:1
                if i+k>=1 && i+k<=height && j+l>=1 && j+l<=width && (k~=0||l~=0)
                    if sup_label(i+k,j+l)~=sup_label(i,j)
                        connection_matrix(sup_label(i+k,j+l), sup_label(i,j)) = 1;
                        connection_matrix(sup_label(i,j), sup_label(i+k,j+l)) = 1;
                    end
                end
            end
        end
    end
end
neigh_info = cell(sup_num, 1);
for i = 1:sup_num
    neigh_info{i, 1} = find(connection_matrix(i, :)==1);
end