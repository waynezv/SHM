function [filterBank] = createFilterBank2()
%% Parameter Setting
k = 0.45;
h_siz = 5;

%% Filter Bank
%% Gaussians
h1 = fspecial2('gaussian', h_siz, k);
h2 = fspecial2('gaussian', h_siz, 2*k);
h3 = fspecial2('gaussian', h_siz, 4*k);
%% Laplacians of Gaussians
h4 = fspecial2('log', h_siz, k);
h5 = fspecial2('log', h_siz, 2*k);
h6 = fspecial2('log', h_siz, 4*k);
h7 = fspecial2('log', h_siz, 8*k);
%% Derivatives of Gaussians
h8 = fspecial2('DoGx', h_siz, 2*k);
h9 = fspecial2('DoGx', h_siz, 4*k);
h10 = fspecial2('DoGy', h_siz, 2*k);
h11 = fspecial2('DoGy', h_siz, 4*k);
filterBank = {h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11};