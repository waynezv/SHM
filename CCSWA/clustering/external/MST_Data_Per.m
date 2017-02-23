function trans_dist_mat = MST_Data_Per(W,pert)
N = size(W, 1);
trans_dist_mat = zeros(N,N);
% tic
[X, Y] = meshgrid(1:N,1:N);
idx_tril = logical(1-tril(ones(N,N)));
Edgelist_Orig = [X(idx_tril) Y(idx_tril) W(idx_tril)];
Edgelist_Pert(:,3) = Edgelist_Orig(:,3) + pert.*rand(size(Edgelist_Orig,1),1);
clear X Y;
% display(['Establishing edge list takes ' num2str(toc) ' seconds']);

% tic
[~,IDX_Sort] = sort(Edgelist_Pert(:,3),'ascend');
Edgelist_Orig = Edgelist_Orig(IDX_Sort,:);
clear IDX_Sort;
% display(['Sorting all edges takes ' num2str(toc) ' seconds']);

label = 1:N;
count = 0;
clust_num = N;
% tic
while(clust_num>1)
    count = count+1;
    if(label(Edgelist_Orig(count, 1))~=label(Edgelist_Orig(count, 2)))
        idx1 = label==label(Edgelist_Orig(count, 1));
        idx2 = label==label(Edgelist_Orig(count, 2));
        label(idx1) = label(Edgelist_Orig(count, 2));
        trans_dist_mat(idx1,idx2) = Edgelist_Orig(count, 3);
        trans_dist_mat(idx2,idx1) = Edgelist_Orig(count, 3);
        clust_num = clust_num-1;
    end
end
% display(['Getting mst edges takes ' num2str(toc) ' seconds']);