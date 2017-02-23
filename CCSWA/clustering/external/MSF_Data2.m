function trans_dist_mat = MSF_Data2(W, tree_num, pooling)
N = size(W, 1);
trans_dist_mat_set = zeros(N,N,tree_num);
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
edge_num = size(Edgelist,1);

label = cell(tree_num,1);
cluster_num = zeros(tree_num,1);
for i = 1:tree_num
    label{i,1} = 1:N;
    cluster_num(i,1) = N;
end

tic

tree_iter = 1;
for i = 1:edge_num
    tmp = 0;
    while tmp<tree_num
        if(label{tree_iter,1}(Edgelist(i, 1))~=label{tree_iter,1}(Edgelist(i, 2)))
            idx1 = label{tree_iter,1}==label{tree_iter,1}(Edgelist(i, 1));
            idx2 = label{tree_iter,1}==label{tree_iter,1}(Edgelist(i, 2));
            label{tree_iter,1}(idx1) = label{tree_iter,1}(Edgelist(i, 2));
            trans_dist_mat_set(idx1,idx2,tree_iter) = Edgelist(i, 3);
            trans_dist_mat_set(idx2,idx1,tree_iter) = Edgelist(i, 3);
            cluster_num(tree_iter,1) = cluster_num(tree_iter,1) - 1;
            tree_iter = mod(tree_iter,tree_num)+1;
            break
        else
            tmp = tmp+1;
            tree_iter = mod(tree_iter,tree_num)+1;
        end
    end
    if(sum(cluster_num~=1)==0)
        break
    end
end
if strcmp(pooling, 'max')
    trans_dist_mat = squeeze(max(trans_dist_mat_set,[],3));
elseif strcmp(pooling, 'mean')
    trans_dist_mat = mean(trans_dist_mat_set,3);
else
    display('Wrong pooling method indication.')
end
display(['Getting mst edges takes ' num2str(toc) ' seconds']);