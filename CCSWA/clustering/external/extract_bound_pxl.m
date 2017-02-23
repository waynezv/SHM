function [bound_pxl_mask] = extract_bound_pxl(seg_label, bandwidth)
[height width] = size(seg_label);
bound_pxl_mask = zeros(height, width);
% h = fspecial('average', bandwidth);
% seg_label2 = imfilter(seg_label, h, 'replicate');
% for i=1:height
%     for j=1:width
%         if (abs(seg_label2(i,j)-seg_label(i,j)) > 0.1/(bandwidth*bandwidth))
%             bound_pxl_mask(i, j) = 1;
%         end
%     end
% end
[X Y] = meshgrid(1:width, 1:height);
for i=1:height
    for j=1:width
        min_height = max(1, i-bandwidth);
        max_height = min(height, i+bandwidth);
        min_width = max(1, j-bandwidth);
        max_width = min(width, j+bandwidth);
        
        tmp = seg_label(min_height:max_height, min_width:max_width)~=seg_label(i, j);
        tmp2 = sqrt((Y(min_height:max_height, min_width:max_width)-i).^2 + (X(min_height:max_height, min_width:max_width)-j).^2) <= bandwidth;
        tmp3 = tmp.*tmp2;
        if(sum(tmp3(:)) ~= 0)
            bound_pxl_mask(i, j) = 1;
        end
    end
end