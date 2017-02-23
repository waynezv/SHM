function dist_mat_out = max_var_anal2(dist_mat_in, max_dim)
[U,S,~] = svd(dist_mat_in,0);
max_svd_dim = size(S,1);
dist_mat_out = U(:,1:min(max_dim,max_svd_dim));