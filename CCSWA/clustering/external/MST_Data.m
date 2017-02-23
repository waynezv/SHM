function trans_dist_mat = MST_Data(W)
N = size(W, 1);
trans_dist_mat = zeros(N,N);
tic
[X, Y] = meshgrid(1:N,1:N);
idx_tril = logical(1-tril(ones(N,N)));
Edgelist = [X(idx_tril) Y(idx_tril) W(idx_tril)];
clear X Y;
display(['Establishing edge list takes ' num2str(toc) ' seconds']);

tic
[~,IDX_Sort] = sort(Edgelist(:,3),'ascend');
Edgelist = Edgelist(IDX_Sort,:);
clear IDX_Sort;
display(['Sorting all edges takes ' num2str(toc) ' seconds']);

label = 1:N;
count = 0;
clust_num = N;
tic
while(clust_num>1)
    count = count+1;
    if(label(Edgelist(count, 1))~=label(Edgelist(count, 2)))
        idx1 = label==label(Edgelist(count, 1));
        idx2 = label==label(Edgelist(count, 2));
        label(idx1) = label(Edgelist(count, 2));
        trans_dist_mat(idx1,idx2) = Edgelist(count, 3);
        trans_dist_mat(idx2,idx1) = Edgelist(count, 3);
        clust_num = clust_num-1;
    end
end
display(['Getting mst edges takes ' num2str(toc) ' seconds']);