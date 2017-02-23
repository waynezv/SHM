function IDX = ncut_clust(W, k)
%% Normalized Cuts
[NcutDiscrete,NcutEigenvectors,NcutEigenvalues] = ncutW(W,k);
N = size(NcutDiscrete,1);
IDX = zeros(N,1);
for j=1:k,
    id = find(NcutDiscrete(:,j));
    IDX(id) = j;
end