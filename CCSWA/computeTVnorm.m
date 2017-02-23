close all, clear all
% load files
addpath ./TVnorm
addpath ../data_mat_t
load V_1M_s0.02_w1.5-0.5-0.25_t0.45.mat
V_cc = abs(V);
load V_bp_1M_ori.mat
V_bp = abs(V);

% min ||F1||tv+||F2||1 s.t. F=F1+F2
lamda = 1;

[V1,V2]=TVL1(V_cc,lamda); % V1 is TV, V2 is sparse, if lambda is large, strengthens sparse
% 
% figure
% imagesc(V_cc, [0 0.5])
% 
% figure
% imagesc(V1, [0 0.5]), title('V1 is TV');
% 
% figure
% imagesc(V2, [0 0.5]), title('V2 is sparse');

% compute TV-norm
[Gmag1, Gdir1] = imgradient(V1);
% [Gmag2, Gdir2] = imgradient(V2);
[grad_x, grad_y]= imgradientxy(V_cc);

%anisotropic TV
TV_ani1 = sum(abs(Gdir1(:)))
% TV_ani2 = sum(abs(Gdir2(:))) 
TVnorm = TVnorm(V_cc)

%isotropic TV
TV_iso = sum(sqrt(grad_x(:).^2 + grad_y(:).^2))
% TV_iso_2 = sum(sqrt(Gdir(:).^2))