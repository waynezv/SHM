function label = AHC(W, desired_clust_num)
N = size(W, 1);
tic
[X, Y] = meshgrid(1:N,1:N);
idx_tril = logical(1-tril(ones(N,N)));
Edgelist = [X(idx_tril) Y(idx_tril) W(idx_tril)];
clear X Y;
toc
display(['Establishing edge list takes ' num2str(toc) ' seconds']);

tic
[~,IDX_Sort] = sort(Edgelist(:,3),'ascend');
Edgelist = Edgelist(IDX_Sort,:);
clear IDX_Sort;
toc
display(['Sorting all edges takes ' num2str(toc) ' seconds']);

label = 1:N;
count = 0;
clust_num = N;
tic
while(clust_num>desired_clust_num)
    count = count+1;
    if(label(Edgelist(count, 1))~=label(Edgelist(count, 2)))
        idx_tmp = label==label(Edgelist(count, 2));
        label(idx_tmp) = label(Edgelist(count, 1));
        clust_num = clust_num-1;
    end
end
toc
display(['Getting mst edges takes ' num2str(toc) ' seconds']);
label = label';