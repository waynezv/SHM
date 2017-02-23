function [feature] = extract_hist(img, label)
[height width channel] = size(img);
R = img(:,:,1);
G = img(:,:,2);
B = img(:,:,3);
sup_num = max(label(:));
hist_set = cell(1, channel);
for i = 1:channel
    hist_set{1, i} = zeros(sup_num, 8);
end

for rr = 1 : height
    for cc = 1 : width
        for i = 1 : channel
            bin_num = ceil((img(rr, cc, i)+1)/32);
            hist_set{1, i}(label(rr, cc), bin_num) = hist_set{1, i}(label(rr, cc), bin_num) + 1;
        end
    end
end
for i = 1:channel
    for j = 1:sup_num
        hist_set{1, i}(j, :) = hist_set{1, i}(j, :)./sum(hist_set{1, i}(j, :));
    end
end

feature = [];
for i = 1:channel
    feature = [feature hist_set{1, i}];
end