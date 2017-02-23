function [D_out,mu_out,std_out] = zscore_normalize(D_in)
[D, N] = size(D_in);
mu_out = mean(D_in,2);
D_out = D_in-repmat(mu_out, [1 N]);
std_out = std(D_out,0,2);
D_out = D_out./repmat(std_out, [1 N]);