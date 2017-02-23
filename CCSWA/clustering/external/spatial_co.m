function [coordinate] = spatial_co(label)
sup_num = max(label(:));
[height width] = size(label);
coordinate = zeros(sup_num, 2);
[X Y] = meshgrid(1:width,1:height);
for i = 1:sup_num
    idx = label==i;
    coordinate(i, 1) = mean(X(idx));
    coordinate(i, 2) = mean(Y(idx));
end