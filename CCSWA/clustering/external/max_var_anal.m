function [output_mat] = max_var_anal(input_mat, max_dim)
sigma = std(input_mat, 1);
[~, idx] = sort(sigma, 'descend');
idx2 = idx(1:max_dim);
output_mat = input_mat(:,idx2);