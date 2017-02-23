function [label] = superpixel_slic(img, reg_size, regularizer)
[height width channel] = size(img);
imglab = single(vl_xyz2lab(vl_rgb2xyz(img)));
label_pre = vl_slic(imglab, reg_size, regularizer);
label_unique = unique(label_pre);
label = zeros(height, width);
for i = 1:length(label_unique)
    idx = label_pre == label_unique(i);
    label(idx) = i;
end
tic
label = eight_connectivity(label);
toc
% tic
% label = Region_Merging(double(img), label, 300, 0);
% toc