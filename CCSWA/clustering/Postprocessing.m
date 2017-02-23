filelist = dir('*.bmp');
for i = 1 : length(filelist)
% for i = 11 : 11
    img = imread(filelist(i).name, 'bmp');
    width = size(img, 2);
    height = size(img, 1);
    bound = zeros(height, width);
    for j = 1 : height
        for k = 1 : width
            if img(j, k, 1)==0 && img(j, k, 2)==0 && img(j, k, 3)==0
                bound(j, k) = 1;
            end
        end
    end
    SE = strel('rectangle',[4 4]);
    bound2 = imdilate(bound,SE);
    for j = 1 : height
        for k = 1 : width
            if bound2(j, k) == 1
                img(j, k, :) = 0;
            end
        end
    end
    imwrite(img, [filelist(i).name(1:end-4) '_post' '.bmp'], 'bmp');
end